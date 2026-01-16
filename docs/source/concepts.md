# Core concepts

## The `Geno` object model

{py:class}`genal.Geno` is the central class. It wraps a SNP-level table and accumulates intermediate results during workflows.

Key attributes (you don't need to manipulate these directly):

- `G.data`: main SNP table (`pandas.DataFrame`)
- `G.phenotype`: set by {py:meth}`genal.Geno.set_phenotype` (phenotype DataFrame + metadata)
- `G.MR_data`: set by {py:meth}`genal.Geno.query_outcome` (exposure/outcome tables used by MR)
- `G.MR_results`: set by {py:meth}`genal.Geno.MR` (results table + harmonized SNP table; used by plotting)
- `G.MR_loo_results`: set by {py:meth}`genal.Geno.MR_loo` (leave-one-out results tuple; used by `MR_loo_plot`)
- `G.MRpresso_results` / `G.MRpresso_subset_data`: set by {py:meth}`genal.Geno.MRpresso`

## Standard columns

Most `genal` workflows assume (a subset of) the following “standard” columns:

| Column | Meaning |
|---|---|
| `CHR` | Chromosome (integer; X is typically encoded as 23) |
| `POS` | Base-pair position |
| `SNP` | Variant identifier (rsID or a `CHR:POS:...` fallback) |
| `EA` | Effect allele |
| `NEA` | Non-effect allele |
| `BETA` | Effect estimate (beta; odds ratios can be log-transformed during preprocessing) |
| `SE` | Standard error of `BETA` |
| `P` | P-value |
| `EAF` | Effect allele frequency (aligned to `EA` when possible) |

### What is required where (rule of thumb)

This is a *practical guide*, not an exhaustive contract. When a method can work with either an rsID or genomic coordinates, this is written as `SNP (or CHR+POS)`.

| Method | Minimal required columns | Notes / recommended inputs |
|---|---|---|
| {py:meth}`genal.Geno.preprocess_data` | partial inputs | Can fill/validate columns using a build-only reference. Filling rsIDs requires `CHR+POS`; filling coordinates requires `SNP`. |
| {py:meth}`genal.Geno.clump` | `SNP`, `P` | LD clumping via PLINK; returns a new `Geno` (or `None` if nothing passes). |
| {py:meth}`genal.Geno.prs` | `EA`, `BETA`, plus `SNP (or CHR+POS)` | If `CHR+POS` are available, genal will prefer position-based matching to your genotype dataset to reduce ID-mismatch losses. |
| {py:meth}`genal.Geno.query_outcome` | `SNP`, `EA`, `NEA`, `BETA`, `SE` (exposure and outcome) | Outcome querying is rsID-based; proxy search is optional. If you plan to use `action=2` later, `EAF` in both datasets is strongly recommended. |
| {py:meth}`genal.Geno.MR` / {py:meth}`genal.Geno.MR_loo` / {py:meth}`genal.Geno.MRpresso` | `MR_data` | All consume `MR_data` produced by `query_outcome()`. |
| {py:meth}`genal.Geno.colocalize` | `BETA`, `SE`, plus `CHR+POS` (preferred) **or** `SNP` (in both datasets) | If `EA/NEA` are present in both datasets, effects are allele-aligned; otherwise results assume both GWAS use the same reference allele. For quantitative traits, provide `sdY` or (`EAF` + `n`) to avoid the default `sdY=1` assumption. |
| {py:meth}`genal.Geno.update_eaf` | `EA`, plus `CHR+POS` **or** `SNP` | Uses PLINK to compute allele frequencies from a reference panel; coordinate-based matching is faster when available. |
| {py:meth}`genal.Geno.filter_by_gene` / {py:meth}`genal.Geno.lift` | `CHR`, `POS` | Genomic coordinate operations. |

## Method behaviors (what you get back)

A helpful mental framework:

- Methods that *select/subset variants* create a new working table and return a **new `Geno`** (so you can keep both input and output to chain methods on both Geno objects).
- Many “table-transformer” utilities return a **`pandas.DataFrame`**. If they accept `replace=`, that flag usually controls whether `G.data` is overwritten; the return type stays a DataFrame. If you want to keep chaining Geno methods, wrap the returned DataFrame with `G.copy(df)`.
- Workflow/analysis steps typically **attach results** on the object (and sometimes also return a summary object).

### Cheat sheet (common methods)

| Method | Returns | Behavior / side effects |
|---|---|---|
| `preprocess_data()` | `None` | mutates `G.data` (clean/fill/validate) |
| `clump()` | `Geno` | returns a new `Geno` (instrument set); original unchanged; uses PLINK temp files |
| `update_snpids()` | `pd.DataFrame` | returns updated SNP IDs; `replace=True` overwrites `G.data` |
| `lift()` | `pd.DataFrame` | returns lifted coordinates; `replace=True` overwrites `G.data`; may download chain files / write outputs |
| `update_eaf()` | `pd.DataFrame` | returns updated `EAF`; `replace=True` overwrites `G.data`; uses PLINK |
| `standardize_betas()` | `pd.DataFrame` | returns standardized effects; `replace=True` overwrites `G.data` |
| `set_phenotype()` | `None` | sets `G.phenotype` (phenotype table + metadata) |
| `association_test()` | `None` | runs PLINK `--glm`; mutates `G.data` (`BETA/SE/P`) |
| `query_outcome()` | `None` | sets `G.MR_data` (exposure/outcome tables used by MR) |
| `MR()` | `pd.DataFrame` | sets `G.MR_results` and returns the results table |
| `MR_plot()` | plot object | requires `G.MR_results`; writes `.png` if `filename=...`; supports `use_mrpresso_data=True` for outlier highlighting |
| `MR_funnel()` | plot object | requires `G.MR_results`; writes `.png` if `filename=...`; supports `use_mrpresso_data=True` for outlier highlighting |
| `MR_loo()` | `pd.DataFrame` | sets `G.MR_loo_results` and returns the LOO results table |
| `MR_loo_plot()` | plot object(s) | requires `G.MR_loo_results`; writes `.png` if `filename=...`; may return a list for multi-page output; supports `methods=[...]` overall rows and `use_mrpresso_data=True` for outlier highlighting |
| `MRpresso()` | tuple | sets `G.MRpresso_results` and `G.MRpresso_subset_data` (outlier-removed harmonized table; SNP-indexed) |
| `prs()` | `None` | writes `<name>.csv` and uses PLINK temp files |
| `query_gwas_catalog()` | `pd.DataFrame` | adds an `ASSOC` column (network-bound); `replace=True` overwrites `G.data` |
| `filter_by_gene(replace=False)` | `Geno` | returns a new `Geno` filtered to a locus |
| `filter_by_gene(replace=True)` | `None` | filters `G.data` in place |
| `colocalize()` | dict | returns posterior probabilities; does not modify `G.data` |
| `save()` | `None` | writes `G.name.(h5|csv|txt)` to disk |

## Side effects and files

Be aware of these common side effects:

- `~/.genal/config.json` is created/updated as you configure PLINK, reference folders, or default genotype paths.
- `tmp_GENAL/` is used as a scratch directory for PLINK commands and is **not** automatically deleted.
- Some methods generate output files in your current directory (notably `prs()`, and plot saving in `MR_plot()`, `MR_funnel()`, and `MR_loo_plot()`).

## Resource usage (`ram`, `cpus`)

`Geno` sets defaults for:

- `G.cpus`: derived from `SLURM_CPUS_PER_TASK` when present, otherwise `os.cpu_count() - 1`
- `G.ram`: derived from `SLURM_MEM_PER_CPU` when present, otherwise from detected system RAM

You can override them after initialization:

```python
G.cpus = 8
G.ram = 25_000  # MB
```

Many PLINK commands accept `--memory` and/or `--threads` parameters that are fed by these attributes.
