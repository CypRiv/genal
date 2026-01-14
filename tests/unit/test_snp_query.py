from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Any

import pytest


@dataclass
class _FakeResponse:
    status: int
    payload: dict[str, Any] | None = None
    delay_s: float = 0.0

    async def __aenter__(self) -> "_FakeResponse":
        return self

    async def __aexit__(self, exc_type, exc, tb) -> bool:
        return False

    async def json(self) -> dict[str, Any]:
        return self.payload or {}


class _FakeClientSession:
    def __init__(self, responses: dict[str, _FakeResponse]) -> None:
        self._responses = responses

    async def __aenter__(self) -> "_FakeClientSession":
        return self

    async def __aexit__(self, exc_type, exc, tb) -> bool:
        return False

    async def get(self, url: str) -> _FakeResponse:
        resp = self._responses[url]
        if resp.delay_s:
            await asyncio.sleep(resp.delay_s)
        return resp


class _NoopTqdm:
    def __init__(self, total: int, desc: str) -> None:
        self.total = total
        self.desc = desc

    def __enter__(self) -> "_NoopTqdm":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        return False

    def update(self, n: int = 1) -> None:
        return


def test_async_query_gwas_catalog_handles_empty_snp_list(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession({}))

    results, errors, timeouts = async_query_gwas_catalog([])
    assert results == {}
    assert errors == []
    assert timeouts == []


def test_async_query_gwas_catalog_returns_correct_structure(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    snp = "rsTEST"
    base_url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations"
        "?projection=associationBySnp"
    )
    payload = {
        "_embedded": {
            "associations": [
                {"pvalue": 1e-9, "efoTraits": [{"trait": "TraitA"}]},
                {"pvalue": 1e-9, "efoTraits": [{"trait": "TraitA"}]},  # duplicated trait
            ]
        }
    }
    responses = {base_url: _FakeResponse(status=200, payload=payload)}

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession(responses))

    results, errors, timeouts = async_query_gwas_catalog([snp], p_threshold=5e-8)
    assert set(results.keys()) == {snp}
    assert set(results[snp]) == {"TraitA"}
    assert errors == []
    assert timeouts == []


def test_async_query_gwas_catalog_p_threshold_filtering(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    snp = "rsP"
    base_url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations"
        "?projection=associationBySnp"
    )
    payload = {
        "_embedded": {
            "associations": [
                {"pvalue": 1e-9, "efoTraits": [{"trait": "Keep"}]},
                {"pvalue": 1e-3, "efoTraits": [{"trait": "Drop"}]},
            ]
        }
    }
    responses = {base_url: _FakeResponse(status=200, payload=payload)}

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession(responses))

    results, errors, timeouts = async_query_gwas_catalog([snp], p_threshold=5e-8)
    assert results[snp] == ["Keep"]
    assert errors == []
    assert timeouts == []


def test_async_query_gwas_catalog_max_associations_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    snp = "rsLIMIT"
    base_url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations"
        "?projection=associationBySnp"
    )
    payload = {
        "_embedded": {
            "associations": [
                {"pvalue": 1e-9, "efoTraits": [{"trait": "T1"}]},
                {"pvalue": 1e-9, "efoTraits": [{"trait": "T2"}]},
                {"pvalue": 1e-9, "efoTraits": [{"trait": "T3"}]},
            ]
        }
    }
    responses = {base_url: _FakeResponse(status=200, payload=payload)}

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession(responses))

    results, errors, timeouts = async_query_gwas_catalog(
        [snp], p_threshold=5e-8, max_associations=2
    )
    assert len(results[snp]) <= 2
    assert errors == []
    assert timeouts == []


def test_async_query_gwas_catalog_handles_timeout(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    snp = "rsTIMEOUT"
    base_url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations"
        "?projection=associationBySnp"
    )
    responses = {base_url: _FakeResponse(status=200, payload={"_embedded": {"associations": []}}, delay_s=0.01)}

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession(responses))

    results, errors, timeouts = async_query_gwas_catalog([snp], timeout=0.001)
    assert results == {}
    assert errors == []
    assert timeouts == [snp]


def test_async_query_gwas_catalog_invalid_snp_is_reported_as_error(monkeypatch: pytest.MonkeyPatch) -> None:
    from genal.snp_query import async_query_gwas_catalog

    snp = "rsINVALID"
    base_url = (
        f"https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{snp}/associations"
        "?projection=associationBySnp"
    )
    responses = {base_url: _FakeResponse(status=404, payload=None)}

    monkeypatch.setattr("genal.snp_query.tqdm", _NoopTqdm)
    monkeypatch.setattr("genal.snp_query.aiohttp.ClientSession", lambda: _FakeClientSession(responses))

    results, errors, timeouts = async_query_gwas_catalog([snp])
    assert results == {}
    assert errors == [snp]
    assert timeouts == []


@pytest.mark.network
def test_async_query_gwas_catalog_real_snp() -> None:
    from genal.snp_query import async_query_gwas_catalog

    results, errors, timeouts = async_query_gwas_catalog(["rs7412"], p_threshold=1, timeout=30)
    assert "rs7412" in results or "rs7412" in errors or "rs7412" in timeouts


@pytest.mark.network
def test_async_query_gwas_catalog_multiple_snps() -> None:
    from genal.snp_query import async_query_gwas_catalog

    snps = ["rs7412", "rs429358"]
    results, errors, timeouts = async_query_gwas_catalog(snps, p_threshold=1, timeout=30)
    for snp in snps:
        assert snp in results or snp in errors or snp in timeouts


@pytest.mark.network
def test_async_query_gwas_catalog_handles_timeout_network() -> None:
    from genal.snp_query import async_query_gwas_catalog

    # Use an extremely low timeout; this is expected to fail via timeout or error.
    _, errors, timeouts = async_query_gwas_catalog(["rs7412"], p_threshold=1, timeout=0.001)
    assert errors or timeouts


@pytest.mark.network
def test_async_query_gwas_catalog_invalid_snp_network() -> None:
    from genal.snp_query import async_query_gwas_catalog

    invalid = "rsTHIS_SNP_DOES_NOT_EXIST_123"
    results, errors, timeouts = async_query_gwas_catalog([invalid], p_threshold=1, timeout=30)
    assert invalid not in results
    assert invalid in errors or invalid in timeouts

