# Tutorial (detailed)

The `README.md` shows a compact “GWAS → instruments → PRS → MR” example. This page expands that workflow into a **tutorial**, explaining *why* each step exists and which parameters you typically tune. The rules applied to your data at each step are described in {doc}`methods`.

If you haven't configured PLINK/reference data yet, start with {doc}`setup`.

## 0) Mental model: what a `Geno` instance represents

{py:class}`genal.Geno` wraps a SNP-level table in `G.data`. You typically:

1. load and standardize GWAS summary statistics into a `Geno`,
2. derive an instrument set (clumping),
3. run downstream analyses (PRS, MR, MR-PRESSO, colocalization, …).

During a workflow, the object evolves in two main ways:

- the working table `G.data` may be transformed (sometimes into a new `Geno`, e.g. clumping),
- results are attached as attributes (e.g. `G.MR_data` after `query_outcome()`, `G.MR_results` after `MR()`).

Many analysis methods both **store** results on the object and **return** a convenient summary object (for example, `MR()` returns a results table and also populates `G.MR_results`).

See {doc}`concepts` for required columns and a concise “what each method returns / modifies” cheat sheet.

## 1) Create a `Geno` from GWAS summary statistics

You map your input columns to `genal`'s standard names at initialization:

```python
import pandas as pd
import genal

exposure_df = pd.read_csv("exposure_gwas.tsv", sep="\t")
G_exposure = genal.Geno(
    exposure_df,
    CHR="CHR", POS="POS", SNP="SNP",
    EA="EA", NEA="NEA",
    BETA="BETA", SE="SE", P="P", EAF="EAF",
    keep_columns=False,  # keep only standard columns after renaming
)
```

Practical notes:
- If your GWAS uses odds ratios, pass `effect_column="OR"` later in `preprocess_data()` (see below).
- If you have either `SNP` *or* `CHR+POS` (but not both), `preprocess_data()` can often fill the missing identifier using a build-only reference.

## 2) Preprocess and standardize the table

Run {py:meth}`genal.Geno.preprocess_data` early. It validates types/values and can fill missing columns using a build-only reference (`reference_panel="37"` or `"38"`).

```python
G_exposure.preprocess_data(preprocessing="Fill_delete", reference_panel="37")
```

Key arguments you commonly tune:
- `preprocessing`:
  - `"Fill"` keeps rows and sets invalid values to `NaN` (useful for exploratory work).
  - `"Fill_delete"` drops invalid/duplicated/NA rows (recommended before downstream genetics methods).
- `reference_panel="37"` vs `"38"` should match your GWAS build.
- `effect_column="OR"` if the `BETA` column holds odds ratios (they are log-transformed). By default, genal detects the scale from the data and prints its decision.
- `fill_snpids` / `fill_coordinates`: override the default logic if you want to force filling rsIDs from `CHR+POS` or vice-versa.
- `keep_indel` / `keep_dups`: keep indels or duplicated IDs (generally you keep these `False` unless you have a reason).
- `fill_f=True`: force recomputation of the F-statistic (`FSTAT`) column even if it already exists. By default, FSTAT is created if missing or only missing values are filled.

The validation and filling rules are listed in [Methods: preprocessing](methods.md#preprocessing).

## 3) Select independent instruments via LD clumping

{py:meth}`genal.Geno.clump` runs PLINK 2 clumping and returns a **new** `Geno` containing the clumped instruments.

```python
G_instruments = G_exposure.clump(
    p1=5e-8,       # significance threshold for index variants
    r2=0.01,       # LD threshold
    kb=10_000,     # window size in kb
    p2=0.01,       # PLINK's secondary threshold (only affects PLINK's clump listing, not the index SNPs returned)
    reference_panel="EUR_37",  # LD reference panel (not the same as build-only reference)
)
```

Notes:
- `reference_panel` here is an **LD reference panel** (`EUR_37`, `AFR_38`, …) or a path to custom PLINK files.
- If no variants remain after clumping, `clump()` returns `None`.

## 4) Compute a PRS in a target cohort (individual-level genotypes)

{py:meth}`genal.Geno.prs` writes a CSV and stores intermediate PLINK outputs under `tmp_GENAL/`.

```python
G_instruments.prs(
    name="prs_sbp",
    path="my_genotypes_chr$",  # bed/bim/fam or pgen/pvar/psam; use '$' if split by chr
    weighted=True,
)
```

The output `<name>.csv` has one row per individual, with the sum score `SCORE1_SUM` and PLINK's average score `SCORE1_AVG`. `BETA` is used as given: for odds ratios, run `preprocess_data(effect_column="OR")` first.

Tunable options and common pitfalls:
- `weighted=False` gives each variant a weight of 1 on its trait-increasing allele (an allele count).
- **Variant matching**:
  - If `CHR+POS` are available in `G_instruments.data`, genal matches variants to your genotype files by position and alleles (keep `NEA` in your data) and calls `update_snpids()` internally to reduce ID mismatches.
  - If you only have rsIDs, ensure your cohort uses compatible IDs (or add coordinates and rerun preprocessing).
- `proxy=True` enables proxy SNP search when variants are missing in your cohort:
  - `r2`, `kb`, `window_snps`, and the LD `reference_panel` control proxy search scope/strictness.
- `name` can be a path prefix; the output is `<name>.csv`.

Example with proxies:

```python
G_instruments.prs(
    name="prs_sbp_proxy",
    path="my_genotypes_chr$",
    proxy=True,
    reference_panel="EUR_37",
    r2=0.8,
    kb=5000,
    window_snps=1_000_000,
)
```

Matching, scoring and proxy rules: [Methods: polygenic risk scores](methods.md#polygenic-risk-scores).

## 5) Two-sample MR: query an outcome GWAS and run MR

### 5.1) Load and preprocess the outcome GWAS

```python
outcome_df = pd.read_csv("outcome_gwas.tsv", sep="\t")
G_outcome = genal.Geno(
    outcome_df,
    CHR="CHR", POS="POS", SNP="SNP",
    EA="EA", NEA="NEA",
    BETA="BETA", SE="SE", P="P", EAF="EAF",
    keep_columns=False,
)
G_outcome.preprocess_data(preprocessing="Fill_delete", reference_panel="37")
```

### 5.2) Query the outcome at the instrument variants

{py:meth}`genal.Geno.query_outcome` extracts the exposure SNPs from an outcome dataset (matching on rsID) and stores the exposure and outcome tables, with the outcome name, in `G_instruments.MR_data`.

```python
G_instruments.query_outcome(
    G_outcome,
    name="stroke",
    proxy=True,
    reference_panel="EUR_37",
)
```

Tunable options:
- `proxy=True` searches proxies when SNPs are missing from the outcome (see [Methods: proxy SNPs](methods.md#proxy-snps)).
- `r2`, `kb`, `window_snps` control proxy search scope/strictness.
- `name` sets the outcome name used in the MR results.

### 5.3) Run MR

{py:meth}`genal.Geno.MR` harmonizes alleles and runs the requested estimators:

```python
res = G_instruments.MR(
    methods=["IVW", "IVW-FE", "WM", "Simple-mode", "Egger"],
    action=2,
    heterogeneity=True,
    odds=False,
    exposure_name="SBP",
    outcome_name="Stroke",
)
```

Key arguments you commonly tune:
- `methods`: choose estimators/sensitivity analyses; see {doc}`methods` for implementation notes.
- Palindromic SNP handling (`action`):
  - `1`: keep palindromes (assumes alleles already aligned)
  - `2`: frequency-based flipping using `EAF` (conservative default; requires `EAF` in both exposure/outcome)
  - `3`: drop all palindromes
- `heterogeneity=True` adds Cochran’s Q diagnostics for methods that support it.
- `odds=True` adds exponentiated effect summaries (odds ratios + 95% CI strings) when appropriate.
- `seed`: random seed for the bootstrap-based methods (reproducible results).

### 5.4) Plot MR results

After `MR()`, you can generate a scatter plot with method lines:

```python
G_instruments.MR_plot(filename="mr_scatter", figure_size=(10, 6))  # saves mr_scatter.png
```

You can also draw a funnel plot of single-SNP ratio estimates (Wald ratios):

```python
G_instruments.MR_funnel(
    methods=["IVW", "WM", "Egger"],  # vertical reference lines (optional)
    filename="mr_funnel",
    figure_size=(10, 6),
)
```

## 6a) MR-PRESSO (outlier detection and distortion testing)

{py:meth}`genal.Geno.MRpresso` runs a parallel MR-PRESSO implementation (at least 4 instruments are needed).

```python
mod_table, GlobalTest, OutlierTest, BiasTest = G_instruments.MRpresso(
    action=2,
    n_iterations=10_000,
    significance_p=0.05,
    cpus=-1,
)
```

What you typically tune:
- `n_iterations`: more iterations → more stable p-values but slower.
- `significance_p`: threshold for global/outlier tests.
- `outlier_test` / `distortion_test`: disable if you only want the global test.
- `seed`: random seed for the simulated datasets (reproducible results).

Output structure:
- `mod_table`: DataFrame with two rows (`"Raw"` and `"Outlier-corrected"`) for the IVW model before and after outlier removal.
- `GlobalTest`: dictionary with `"RSSobs"` (observed residual sum of squares) and `"global_test_p"` (p-value for the global pleiotropy test).
- `OutlierTest`: DataFrame with per-SNP outlier p-values (`RSSobs`, `Pvalue`); SNP identifiers (rsIDs) are used as row labels. Empty if the global test is not significant.
- `BiasTest`: dictionary containing `"outliers_indices"` (list of SNP IDs), `"distortion_test_coefficient"` (percent), and `"distortion_test_p"`.

To highlight MR-PRESSO outliers on plots, pass `use_mrpresso_data=True` (outliers are colored in red and shown in the legend). The plots show the instruments of the last `MR()` call, so draw them before rerunning MR without the outliers:

```python
G_instruments.MR_plot(filename="mr_scatter_mrpresso", figure_size=(10, 6), use_mrpresso_data=True)
G_instruments.MR_funnel(filename="mr_funnel_mrpresso", figure_size=(10, 6), use_mrpresso_data=True)
```

If outliers are found, you can then rerun MR using the outlier-removed subset:

```python
res_no_outliers = G_instruments.MR(use_mrpresso_data=True, exposure_name="SBP", outcome_name="Stroke")
```

See [Methods: MR-PRESSO](methods.md#mr-presso-summary) for algorithm details, and for when each test runs and what the outputs contain otherwise.

## 6b) Leave-one-out MR (sensitivity analysis)

{py:meth}`genal.Geno.MR_loo` iteratively removes each SNP and re-estimates the causal effect. This helps identify influential variants that may be driving the overall result.

```python
loo_df = G_instruments.MR_loo(
    method="IVW",        # any single MR method key (see MR method map)
    action=2,
    heterogeneity=False, # set True to include Q statistics
    odds=False,          # set True for OR-scale output
)
```

Key arguments:
- `method`: a single MR method key (e.g., `"IVW"`, `"Egger"`, `"WM"`); must not be `"all"`.
- `use_mrpresso_data=True`: use the outlier-removed dataset from MR-PRESSO instead of all instruments.
- `seed`: random seed for bootstrap-based methods.

### Visualizing leave-one-out results

{py:meth}`genal.Geno.MR_loo_plot` creates a forest plot from the stored `MR_loo_results`:

```python
# Default: show top influential instruments
G_instruments.MR_loo_plot(filename="loo_forest", figure_size=(10, 8))
```

```python
# Or paginate all instruments
G_instruments.MR_loo_plot(
    top_influential=False,  # show all, not just influential
    snps_per_page=30,
    page=1,                 # or None for all pages
    filename="loo_forest_all",
    figure_size=(10, 12),
)
```

```python
# Highlight MR-PRESSO outliers (requires MRpresso() first)
G_instruments.MR_loo_plot(filename="loo_forest_mrpresso", figure_size=(10, 8), use_mrpresso_data=True)
```

Key arguments:
- `top_influential=True` (default): select the `snps_per_page` most influential SNPs (largest change in estimate when removed) and render a single compact figure.
- `top_influential=False`: paginate all instruments; use `page=N` to select a specific page or `page=None` to render all pages.
- `snps_per_page`: number of SNPs per page (minimum 5).
- `use_mrpresso_data=True`: color MR-PRESSO outliers in red (requires `MRpresso()` first). When outliers exist and the leave-one-out analysis used all instruments, an extra summary row ("MR-PRESSO corrected") is added using the same MR method.
- `methods=["WM", "Egger"]`: add extra overall estimates for the requested methods (computed on all instruments).


## 7) Additional capabilities (beyond the core pipeline)

### Single-SNP association tests (individual-level data)

Use this when you have **individual-level genotypes** (PLINK files) and a phenotype table for the same cohort, and you want to estimate SNP–phenotype associations in that cohort (linear regression for quantitative traits, logistic regression for binary traits via PLINK 2 `--glm`).

1) Attach phenotype data to the object:

```python
G_instruments.set_phenotype(
    pheno_df,
    PHENO="sbp",
    PHENO_type="quant",  # or "binary"; if omitted, genal tries to infer it
    IID="IID",
    FID="FID",           # optional; defaults to IID
)
```

For binary phenotypes, `set_phenotype()` recodes the two observed values to `0/1` (control/case) and prints the coding it used; use `alternate_control=True` to invert it.

2) Run association tests (`G.data` is replaced by the tested variants, with their cohort estimates):

```python
G_instruments.association_test(
    path="my_genotypes_chr$",  # bed/bim/fam or pgen/pvar/psam; use '$' if split by chr
    covar=["age", "sex"],      # used to specify covariates to adjust for
    standardize=True,          # standardize quantitative phenotypes before regression
)
```

What you typically tune / watch:
- `covar`: covariate names (must be present in `pheno_df`; text or categorical covariates are encoded as indicator columns, and constant covariates are dropped).
- `standardize=True`: for quantitative traits; set `False` if you want raw-scale effects.
- `EA` is required, and `NEA` recommended (it selects the allele pair tested at multi-allelic variants).
- Variant matching: if `CHR+POS` are present in `G.data`, genal maps them to cohort SNP IDs by position and alleles before running PLINK (reduces ID mismatch losses).
- After updating `BETA`, `SE`, and `P`, the F-statistic (`FSTAT`) column is automatically recomputed to remain consistent with the updated estimates.

Allele alignment and filtering rules: [Methods: association tests](methods.md#association-tests).

### Liftover between builds

{py:meth}`genal.Geno.lift` lifts `CHR`/`POS` coordinates across genome builds and returns a **DataFrame** (it does not return a new `Geno`).
The appropriate chain file is downloaded automatically if used for the first time.

```python
lifted_df = G_exposure.lift(start="hg19", end="hg38", replace=False)  # G_exposure.data unchanged
```

Key arguments:
- `start` / `end`: `"hg17"`, `"hg18"`, `"hg19"` or `"hg38"` (`"GRCh37"` and `"GRCh38"` are accepted for `"hg19"` and `"hg38"`).
- `replace=True|False`: whether to overwrite `G.data` (return value is still a DataFrame either way).
- `chain_file=...`: use a local UCSC chain file (useful for offline/reproducible builds).
- `liftover_path=...`: use UCSC `liftOver` executable (faster than Python liftover for large datasets).
- `name=...` / `extraction_file=True`: optionally write the lifted table and/or an extraction file.

Variants that cannot be lifted are removed (their number is printed); see [Methods: liftover](methods.md#liftover) for strand and contig handling.

If you want a lifted `Geno` object (to keep chaining methods), use:

```python
G_lifted = G_exposure.copy(lifted_df)
```

### Update allele frequencies (EAF) from a reference panel

Allele frequencies are useful for MR harmonization with `action=2` (frequency-based palindrome handling) and as a sanity-check for allele alignment.

```python
G_exposure.update_eaf(reference_panel="EUR_37", replace=True, fill=True)
```

Key arguments:
- `reference_panel`: an LD reference panel name (e.g. `"EUR_37"`, `"AFR_38"`) or a custom PLINK fileset (no extension).
- `replace=True|False`: whether to overwrite `G.data` (return value is always a DataFrame).
- `fill=True|False`: keep existing `EAF` values for variants not found in the reference (`True`), or set them to `NaN` (`False`).

Requirements / behavior:
- Requires `EA` and either `SNP` or `CHR+POS` (`NEA` recommended).
- Uses PLINK `--freq` and aligns the returned frequency to the **effect allele** (`EA`). Variants whose alleles do not match the reference panel get `NaN` (or keep their value if `fill=True`); see [Methods: allele frequencies](methods.md#allele-frequencies).

### Colocalization (ABF)

Colocalization is typically run on a locus, not genome-wide. A common pattern is to subset both datasets first:

```python
G_exp_region = G_exposure.filter_by_gene(gene="APOE", window_size=500_000, build="37", replace=False)
G_out_region = G_outcome.filter_by_gene(gene="APOE", window_size=500_000, build="37", replace=False)

coloc = G_exp_region.colocalize(
    G_out_region,
    trait1_type="quant",
    trait2_type="cc",
    n1=700_000,
    n2=40_000,
)
```

Key arguments:
- `trait1_type` / `trait2_type`: `"quant"` or `"cc"`.
- For quantitative traits: provide `sdY1`/`sdY2`, or `n1`/`n2` with an `EAF` column; otherwise $sd_Y = 1$ is assumed (with a warning).
- `p1`, `p2`, `p12`: priors for the ABF model (see {doc}`methods`).
- Matching: by default, genal attempts to merge on `CHR+POS` when available; set `merge_on_snp=True` to force merging by rsID (`SNP`).

Output:
- A dictionary containing `nsnps` and posterior probabilities (`PP.H0.abf` … `PP.H4.abf`).

### GWAS Catalog annotation

{py:meth}`genal.Geno.query_gwas_catalog` annotates variants with traits from the GWAS Catalog REST API (network-bound; can be slow for large SNP sets).

```python
annotated_df = G_exposure.query_gwas_catalog(
    p_threshold=5e-8,
    max_associations=5,
    replace=False,
)
```

Key arguments:
- `p_threshold`: only report associations with $p \le p_\mathrm{threshold}$.
- `max_associations`: cap the number of reported trait-study associations per SNP (the most significant ones are kept).
- `timeout`: per-request timeout; increase for large SNP sets.
- `replace=True` (default) adds the result columns to `G.data`; with `replace=False` they are only added to the returned copy.

The returned table contains:
- `ASSOC`: a list of dictionaries with keys `trait`, `p_value`, and `study_accession`
- `ASSOC_STATUS`: one of `ok`, `no_associations`, `invalid_format`, `not_found`, `timeout`, `rate_limited`, `server_error`, or `client_error`
- `ASSOC_ERROR`: a short error string or `None`

### Gene-window filtering

```python
apoe_region = G_exposure.filter_by_gene(
    gene="APOE",
    id_type="symbol",
    window_size=500_000,
    build="37",
    replace=False,
)
```

Key arguments:
- `gene`: gene identifier (e.g. `"APOE"`; text identifiers are matched case-insensitively).
- `id_type`: identifier type (`"symbol"`, `"HGNC"`, `"Ensembl"`, `"NCBI"`, `"UCSC"`, `"Vega"`, or `"name"`).
- `window_size`: total window size in base pairs (± `window_size/2` around the gene).
- `build`: `"37"` or `"38"` (required; must match your `CHR/POS` coordinates).
- `replace=False|True`: return a new `Geno` (default) or filter this object in place.

Notes:
- Requires `CHR` and `POS`.
- On first use, genal downloads a small gene-mapping file into your reference folder (network required).

### Save intermediate datasets

Use {py:meth}`genal.Geno.save` to persist `G.data` to disk for reuse.

```python
G_instruments.name = "my_instruments"  # controls the filename prefix
G_instruments.save()                   # writes my_instruments.parquet
```

Notes:
- `save()` writes the file using `G.name`; set it explicitly if you want stable filenames.
- Supported formats are `fmt="parquet"` (default), `"csv"`, and `"txt"`.
- You can later load with `pd.read_parquet("my_instruments.parquet")`, or pass the `.parquet` path directly to `query_outcome()`.
