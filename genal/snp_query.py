from __future__ import annotations

import aiohttp
import asyncio
import nest_asyncio
import random
import re
import threading
import time
from collections import Counter, deque
from decimal import Decimal, InvalidOperation
from typing import Any
from urllib.parse import urlencode

from tqdm import tqdm

GWAS_V2_ASSOCIATIONS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
GWAS_V2_SNP_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/{rs_id}"
# Pages of 100 associations: larger pages time out for heavily studied SNPs
PAGE_SIZE = 100
MAX_RETRIES = 3
CLIENT_RATE_LIMIT_QPS = 14
MAX_CONCURRENT_REQUESTS = 12
RETRYABLE_STATUS_CODES = {408, 429, 500, 502, 503, 504}
SERVER_ERROR_STATUS_CODES = {500, 502, 503, 504}
STATUS_ORDER = [
    "ok",
    "no_associations",
    "invalid_format",
    "not_found",
    "timeout",
    "rate_limited",
    "server_error",
    "client_error",
]
_RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


def _ensure_nested_event_loop_support():
    """Patch a running event loop only when needed, typically in notebooks."""
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        return
    if loop.is_running():
        nest_asyncio.apply()


def _normalize_rsid(snp: Any) -> str | None:
    text = str(snp).strip()
    if not text:
        return None
    if not _RSID_PATTERN.fullmatch(text):
        return None
    return f"rs{text[2:]}"


def _build_url(base_url: str, params: dict[str, Any]) -> str:
    filtered = {key: value for key, value in params.items() if value is not None}
    return f"{base_url}?{urlencode(filtered)}"


def _coerce_p_value(value: Any) -> float | None:
    try:
        p_value = float(value)
    except (TypeError, ValueError):
        return None
    if p_value != p_value:
        return None
    return p_value


def _exact_p_value(association: dict[str, Any]) -> Decimal | None:
    """
    Exact p-value of an association, from its mantissa and exponent when available
    (the float p_value field underflows to 0 for p-values below ~1e-308).
    """
    mantissa = association.get("pvalue_mantissa")
    exponent = association.get("pvalue_exponent")
    if mantissa is not None and exponent is not None:
        try:
            p_exact = Decimal(f"{mantissa}e{int(exponent)}")
            if p_exact.is_finite() and p_exact >= 0:
                return p_exact
        except (InvalidOperation, TypeError, ValueError):
            pass
    p_value = _coerce_p_value(association.get("p_value"))
    return None if p_value is None else Decimal(repr(p_value))


def _make_outcome(
    query_snp: str | None,
    associations: list[dict[str, Any]] | None,
    status: str,
    error: str | None = None,
) -> dict[str, Any]:
    return {
        "query_snp": query_snp,
        "associations": associations or [],
        "status": status,
        "error": error,
    }


def _clone_outcome(outcome: dict[str, Any]) -> dict[str, Any]:
    return {
        "query_snp": outcome["query_snp"],
        "associations": [dict(record) for record in outcome["associations"]],
        "status": outcome["status"],
        "error": outcome["error"],
    }


def _format_http_error(status_code: int, payload: Any) -> str:
    if isinstance(payload, dict):
        message = payload.get("errorMessage") or payload.get("message") or payload.get("error")
        if message:
            return f"HTTP {status_code}: {message}"
    return f"HTTP {status_code}"


def _classify_transport_failure(
    status_code: int | None,
    failure_type: str | None,
) -> str:
    if failure_type == "timeout" or status_code == 408:
        return "timeout"
    if status_code == 429:
        return "rate_limited"
    if failure_type == "server_error" or status_code in SERVER_ERROR_STATUS_CODES:
        return "server_error"
    return "client_error"


def _retry_delay_seconds(attempt_index: int, retry_after: str | None = None) -> float:
    if retry_after is not None:
        try:
            return max(float(retry_after), 0.0)
        except ValueError:
            pass
    base_delay = min(0.5 * (2**attempt_index), 4.0)
    return base_delay + random.uniform(0.0, 0.25)


class _RateLimiter:
    def __init__(self, rate: int, period_s: float = 1.0) -> None:
        self.rate = rate
        self.period_s = period_s
        self._timestamps: deque[float] = deque()
        self._lock = threading.Lock()

    async def acquire(self) -> None:
        while True:
            with self._lock:
                now = time.monotonic()
                while self._timestamps and now - self._timestamps[0] >= self.period_s:
                    self._timestamps.popleft()
                if len(self._timestamps) < self.rate:
                    self._timestamps.append(now)
                    return
                wait_s = max(self.period_s - (now - self._timestamps[0]), 0.0)
            await asyncio.sleep(wait_s)


_SHARED_RATE_LIMITER = _RateLimiter(CLIENT_RATE_LIMIT_QPS)


async def _read_json_response(
    session: aiohttp.ClientSession,
    url: str,
    timeout: float | None,
) -> tuple[int, dict[str, str], Any]:
    async def _request() -> tuple[int, dict[str, str], Any]:
        response = await session.get(url)
        async with response:
            try:
                payload = await response.json()
            except (aiohttp.ContentTypeError, ValueError):
                # Error responses can have a non-JSON body (e.g. an HTML error page): they are classified by status code
                if response.status == 200:
                    raise
                payload = None
            return response.status, dict(getattr(response, "headers", {})), payload

    if timeout is None:
        return await _request()
    return await asyncio.wait_for(_request(), timeout=timeout)


async def _fetch_json_with_retries(
    session: aiohttp.ClientSession,
    url: str,
    timeout: float | None,
    semaphore: asyncio.Semaphore,
    rate_limiter: _RateLimiter,
    retry_counts: Counter,
) -> dict[str, Any]:
    last_status = None
    last_error = None
    last_failure_type = None

    for attempt in range(MAX_RETRIES + 1):
        try:
            await rate_limiter.acquire()
            async with semaphore:
                status_code, headers, payload = await _read_json_response(session, url, timeout)

            if status_code == 200:
                return {
                    "ok": True,
                    "status_code": status_code,
                    "payload": payload,
                    "error": None,
                    "failure_type": None,
                }

            last_status = status_code
            last_error = _format_http_error(status_code, payload)
            last_failure_type = (
                "timeout"
                if status_code == 408
                else "rate_limited"
                if status_code == 429
                else "server_error"
                if status_code in SERVER_ERROR_STATUS_CODES
                else "client_error"
            )

            if status_code in RETRYABLE_STATUS_CODES and attempt < MAX_RETRIES:
                retry_counts[_classify_transport_failure(status_code, last_failure_type)] += 1
                await asyncio.sleep(_retry_delay_seconds(attempt, headers.get("Retry-After")))
                continue

            return {
                "ok": False,
                "status_code": status_code,
                "payload": payload,
                "error": last_error,
                "failure_type": last_failure_type,
            }
        except asyncio.TimeoutError:
            last_status = None
            last_error = "Request timed out"
            last_failure_type = "timeout"
            if attempt < MAX_RETRIES:
                retry_counts["timeout"] += 1
                await asyncio.sleep(_retry_delay_seconds(attempt))
                continue
            break
        except aiohttp.ClientError as exc:
            last_status = None
            last_error = str(exc) or exc.__class__.__name__
            last_failure_type = "client_error"
            if attempt < MAX_RETRIES:
                retry_counts["client_error"] += 1
                await asyncio.sleep(_retry_delay_seconds(attempt))
                continue
            break
        except ValueError as exc:
            last_status = None
            last_error = f"Invalid JSON response: {exc}"
            last_failure_type = "server_error"
            if attempt < MAX_RETRIES:
                retry_counts["server_error"] += 1
                await asyncio.sleep(_retry_delay_seconds(attempt))
                continue
            break

    return {
        "ok": False,
        "status_code": last_status,
        "payload": None,
        "error": last_error,
        "failure_type": last_failure_type,
    }


def _association_sort_key(entry: tuple[Decimal, dict[str, Any]]) -> tuple[Decimal, str, str]:
    p_exact, record = entry
    return (
        p_exact,
        record["trait"],
        record["study_accession"] or "",
    )


def _update_best_association(
    best_associations: dict[tuple[str, str | None], tuple[Decimal, dict[str, Any]]],
    trait: str,
    p_exact: Decimal,
    p_value: float,
    study_accession: str | None,
) -> None:
    """Keep, for each (trait, study) pair, the association with the smallest exact p-value."""
    key = (trait, study_accession)
    current = best_associations.get(key)
    if current is None or p_exact < current[0]:
        best_associations[key] = (
            p_exact,
            {
                "trait": trait,
                "p_value": p_value,
                "study_accession": study_accession,
            },
        )


def _finalize_associations(
    best_associations: dict[tuple[str, str | None], tuple[Decimal, dict[str, Any]]],
    max_associations: int | None,
) -> list[dict[str, Any]]:
    """Sort the associations by exact p-value (then trait and study) and keep the first max_associations."""
    ordered = [record for _, record in sorted(best_associations.values(), key=_association_sort_key)]
    if max_associations is not None:
        ordered = ordered[:max_associations]
    return ordered


async def _query_single_rsid(
    session: aiohttp.ClientSession,
    rs_id: str,
    p_threshold: float,
    max_associations: int | None,
    timeout: float | None,
    semaphore: asyncio.Semaphore,
    rate_limiter: _RateLimiter,
    retry_counts: Counter,
) -> dict[str, Any]:
    # Pages are requested in the default API order (association ID, unique), which keeps the paging deterministic:
    # sorting by p-value server-side has ties, and associations can then be duplicated or lost at page boundaries.
    # All pages are read; filtering, deduplication and ranking use the exact p-value (mantissa and exponent).
    def page_url(page: int) -> str:
        return _build_url(GWAS_V2_ASSOCIATIONS_URL, {"rs_id": rs_id, "page": page, "size": PAGE_SIZE})

    async def fetch(url: str) -> dict[str, Any]:
        return await _fetch_json_with_retries(session, url, timeout, semaphore, rate_limiter, retry_counts)

    # Once the number of pages is known from the first page, the other pages are requested concurrently
    # (through the shared semaphore and rate limiter) and processed in page order.
    responses = [await fetch(page_url(0))]
    first_payload = (responses[0]["payload"] or {}) if responses[0]["ok"] else {}
    total_pages = first_payload.get("page", {}).get("totalPages") if isinstance(first_payload, dict) else None
    if isinstance(total_pages, int) and total_pages > 1:
        responses += await asyncio.gather(*(fetch(page_url(page)) for page in range(1, total_pages)))
    elif responses[0]["ok"] and total_pages is None:
        # Without page information, the next links are followed
        next_url = first_payload.get("_links", {}).get("next", {}).get("href")
        while next_url:
            responses.append(await fetch(next_url))
            if not responses[-1]["ok"]:
                break
            next_url = (responses[-1]["payload"] or {}).get("_links", {}).get("next", {}).get("href")

    threshold = Decimal(repr(float(p_threshold)))
    best_associations: dict[tuple[str, str | None], tuple[Decimal, dict[str, Any]]] = {}
    saw_association_rows = False

    for response in responses:
        if not response["ok"]:
            status = _classify_transport_failure(response["status_code"], response["failure_type"])
            return _make_outcome(rs_id, [], status, response["error"])

        payload = response["payload"] or {}
        page_associations = payload.get("_embedded", {}).get("associations", [])
        if page_associations:
            saw_association_rows = True

        for association in page_associations:
            p_exact = _exact_p_value(association)
            if p_exact is None or p_exact > threshold:
                continue
            p_value = _coerce_p_value(association.get("p_value"))
            if p_value is None:
                p_value = float(p_exact)

            study_accession = association.get("accession_id")
            for trait in association.get("efo_traits", []):
                trait_name = trait.get("efo_trait")
                if not trait_name:
                    continue
                _update_best_association(
                    best_associations,
                    trait_name,
                    p_exact,
                    p_value,
                    study_accession,
                )

    if best_associations:
        return _make_outcome(rs_id, _finalize_associations(best_associations, max_associations), "ok")

    if saw_association_rows:
        return _make_outcome(rs_id, [], "no_associations")

    existence_url = GWAS_V2_SNP_URL.format(rs_id=rs_id)
    existence = await _fetch_json_with_retries(
        session,
        existence_url,
        timeout,
        semaphore,
        rate_limiter,
        retry_counts,
    )
    if existence["ok"]:
        return _make_outcome(rs_id, [], "no_associations")
    if existence["status_code"] == 404:
        return _make_outcome(rs_id, [], "not_found", existence["error"])

    status = _classify_transport_failure(existence["status_code"], existence["failure_type"])
    return _make_outcome(rs_id, [], status, existence["error"])


def _build_summary(
    outcomes: list[dict[str, Any]],
    requested_count: int,
    unique_query_count: int,
    retry_counts: Counter,
) -> dict[str, Any]:
    status_counts = Counter(outcome["status"] for outcome in outcomes)
    return {
        "requested_count": requested_count,
        "valid_count": sum(outcome["status"] != "invalid_format" for outcome in outcomes),
        "unique_query_count": unique_query_count,
        "invalid_count": status_counts.get("invalid_format", 0),
        "retry_counts": {status: retry_counts.get(status, 0) for status in STATUS_ORDER if retry_counts.get(status, 0)},
        "retry_total": sum(retry_counts.values()),
        "status_counts": {status: status_counts.get(status, 0) for status in STATUS_ORDER},
    }


def async_query_gwas_catalog(
    snps,
    p_threshold=5e-8,
    max_associations=None,
    timeout=100,
):
    _ensure_nested_event_loop_support()
    try:
        loop = asyncio.get_event_loop()
    except RuntimeError:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
    return loop.run_until_complete(
        query_gwas_catalog_coroutine(
            snps,
            p_threshold=p_threshold,
            max_associations=max_associations,
            timeout=timeout,
        )
    )


async def query_gwas_catalog_coroutine(
    snps,
    p_threshold=5e-8,
    max_associations=None,
    timeout=100,
):
    """
    Query the GWAS Catalog v2 API for SNP associations.

    Parameters:
        snps (list): List of SNPs to query.
        p_threshold (float): Inclusive p-value threshold for filtering associations (compared with the exact
            p-value given by the mantissa and exponent of each association).
        max_associations (int): Maximum number of qualifying trait-study associations to return for each SNP
            (those with the smallest exact p-values).
        timeout (int | None): Timeout for each HTTP request in seconds.

    Returns:
        tuple[list[dict], dict]: Per-input outcomes and query summary metadata.
    """
    prepared_inputs = []
    unique_rsids = []
    seen_rsids = set()
    for snp in snps:
        normalized = _normalize_rsid(snp)
        prepared_inputs.append({"input_snp": snp, "query_snp": normalized})
        if normalized and normalized not in seen_rsids:
            seen_rsids.add(normalized)
            unique_rsids.append(normalized)

    retry_counts: Counter = Counter()
    cached_outcomes: dict[str, dict[str, Any]] = {}

    if unique_rsids:
        client_timeout = aiohttp.ClientTimeout(total=timeout)
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_REQUESTS)
        rate_limiter = _SHARED_RATE_LIMITER
        connector = aiohttp.TCPConnector(
            limit=MAX_CONCURRENT_REQUESTS,
            limit_per_host=MAX_CONCURRENT_REQUESTS,
        )

        async with aiohttp.ClientSession(timeout=client_timeout, connector=connector) as session:
            async def _query_and_tag(rs_id: str) -> tuple[str, dict[str, Any]]:
                return (
                    rs_id,
                    await _query_single_rsid(
                        session,
                        rs_id,
                        p_threshold,
                        max_associations,
                        timeout,
                        semaphore,
                        rate_limiter,
                        retry_counts,
                    ),
                )

            tasks = [asyncio.create_task(_query_and_tag(rs_id)) for rs_id in unique_rsids]
            with tqdm(total=len(tasks), desc="Processing GWAS Catalog rsIDs") as progress_bar:
                for task in asyncio.as_completed(tasks):
                    rs_id, outcome = await task
                    cached_outcomes[rs_id] = outcome
                    progress_bar.update(1)

    outcomes = []
    for item in prepared_inputs:
        normalized = item["query_snp"]
        if normalized is None:
            outcomes.append(
                _make_outcome(
                    None,
                    [],
                    "invalid_format",
                    "Expected a canonical rsID matching rs + digits.",
                )
            )
            continue
        outcomes.append(_clone_outcome(cached_outcomes[normalized]))

    summary = _build_summary(outcomes, len(snps), len(unique_rsids), retry_counts)
    return outcomes, summary
