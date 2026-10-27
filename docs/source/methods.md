# Methods (implementation notes)

This page documents what is implemented in the codebase: the rules applied to your data at each step of the workflows, and the statistical models. The {doc}`workflows` tutorial links here for these details.

## Preprocessing

Implementation: {py:meth}`genal.Geno.preprocess_data` (helpers in `genal/geno_tools.py`).

- **Effect scale** (`effect_column`): `"OR"` log-transforms the `BETA` column (odds ratios ≤ 0 give `NaN`). With the default (`None`), a column with values ≤ 0 holds betas; for an all-positive column with `SE` and `P` available, genal checks whether `BETA/SE` or `log(BETA)/SE` reproduces the z-statistics implied by `P` and prints the evidence; without `SE` or `P`, a heuristic is used (median between 0.5 and 1.5 → odds ratios). For odds ratios, `SE` is kept when it is on the log-OR scale (checked against `P` when available) and recomputed from `BETA` and `P` otherwise.
- **Alleles**: `EA` and `NEA` can be filled from the reference only when there is no `BETA` column, since the allele an effect refers to cannot be inferred; with `BETA` and no `EA`, `preprocess_data()` raises an error. With `EA` present, `NEA` can be filled from `CHR`/`POS`.
- **Validation**: placeholder SNP IDs (`.`, `NA`, empty strings, …) are treated as missing IDs; `EAF` values outside [0, 1] are set to missing; `CHR` values `X`, `Y`, `XY`, `MT`/`M` (with or without a `chr` prefix) are converted to 23, 24, 25, 26; `POS = 0` is set to missing. Counts are printed.
- **Identifiers**: when rsIDs are filled from `CHR`/`POS` and `EA`/`NEA` are present, the reference variant must have the same alleles, and variants not found get the ID `CHR:POS:A1:A2` (alleles in alphabetical order). Positions with several reference variants (without alleles to tell them apart) and rsIDs found at several positions are left missing.
- **Duplicates** (removed under `"Fill_delete"`): the first row of each SNP ID is kept; rows with a missing ID are never duplicates.
- **Missing `P` or `SE`**: `P = 2·Φ(−|BETA/SE|)` is computed with the survival function (exact for very small p-values), and `SE = |BETA| / Φ⁻¹(1 − P/2)` is set to `NaN` when it is not identifiable (`BETA = 0`, `P = 0` or `P = 1`).
- **Index**: `Geno()` does not keep the index of the input DataFrame (`G.data` starts with a 0…n−1 index, and rows deleted later keep their labels); rows can be linked to the input through `SNP` or `CHR`/`POS`.

## Per-variant F-statistic (FSTAT)

Implementation: {py:func}`genal.geno_tools.fill_fstatistic`

`genal` computes a per-variant F-statistic (`FSTAT`) during preprocessing and after association testing. This statistic measures instrument strength for each variant.

- **Primary:** `FSTAT = (BETA / SE)²` when `BETA` and `SE` are available and `SE > 0`.
- **Fallback:** `FSTAT = χ²_isf(P, df=1)` (equivalent to `Z²` for a two-sided p-value) when the primary route is not computable (missing `BETA`/`SE`, or `SE ≤ 0`) but `P` is present.

## Polygenic risk scores

Implementation: {py:meth}`genal.Geno.prs` (helpers in `genal/extract_prs.py`).

- **Scores**: PLINK 2 `--score`. The output CSV has one row per individual with the sample IDs (`FID` when present in the genotype files, `IID`), `ALLELE_CT`, `NAMED_ALLELE_DOSAGE_SUM`, `SCORE1_SUM` (Σ BETA × dosage of the scored allele) and `SCORE1_AVG` (`SCORE1_SUM` divided by the number of allele observations, 2 × the number of scored variants).
- **Weights**: `BETA` is used as given. With `weighted=False`, each variant is scored on its trait-increasing allele with weight 1 (`EA` if `BETA ≥ 0`, `NEA` if `BETA < 0`), so the sum score counts trait-increasing alleles; without a `NEA` column, a variant with `BETA < 0` is scored on `EA` with weight −1.
- **Variant matching**: with `CHR`/`POS`, variants are matched to the genotype files by position and alleles (`EA`/`NEA`, also on the opposite strand for non-palindromic SNVs) through `update_snpids()`. A row keeps its own SNP ID when that ID exists in the genotype files, and at multi-allelic sites the record whose alleles match `EA`/`NEA` is used. If fewer than 90% of the rows that have a genotype variant at their position match its alleles, a genome build mismatch warning is printed. Rows with a missing SNP ID are not scored.
- **Genotype files split by chromosome** (`$` placeholder): chromosomes 23–26 can be named with their number or with `X`, `Y`, `XY`, `MT` (e.g. `my_genotypes_chrX`). If the files contain different samples, a warning is printed: a sample absent from a chromosome file gets that chromosome mean-imputed by PLINK. If the extraction of a chromosome fails (after one retry with fewer workers), the analysis stops with an error naming the chromosome and its log file.
- **Proxies** (`proxy=True`): proxies must be present in the genotype files; palindromic proxies (A/T, C/G) and variants already in the score are not used, and each proxy replaces at most one missing SNP (assigned by decreasing $|r|$). A missing SNP whose `EA`/`NEA` do not match the reference-panel alleles is not proxied. A proxy is scored on its reference-panel minor allele, with the weight of the missing SNP (sign aligned through the LD phase).

## Association tests

Implementation: {py:meth}`genal.Geno.set_phenotype` and {py:meth}`genal.Geno.association_test` (helpers in `genal/association.py`).

- **Binary phenotypes**: the two observed values are recoded to 0/1 (control/case). Values coded 0/1 keep 0 as control, values coded 1/2 use 1 as control, and any other two-value coding uses the most frequent value as control; `alternate_control=True` inverts the coding.
- **Model**: PLINK 2 `--glm` (linear regression for quantitative traits, logistic regression for binary traits). Covariates are variance-standardized by PLINK (`--covar-variance-standardize`), which does not change the SNP estimates. Phenotype columns already present in the genotype `.psam` file are ignored; the phenotype set with `set_phenotype()` is used.
- **Covariates**: as in a complete-case regression, the analysed individuals are the samples of the genotype files with a phenotype and all the requested covariates (passed to PLINK with `--keep`, also when no covariate column remains); covariates with a single value among them are dropped. Text or categorical covariates are encoded as $k-1$ indicator (0/1) columns named `<covariate>=<level>`, with the first level in sorted order as the reference (the levels are printed); covariate names are made PLINK-compatible and unique. Numeric covariates are continuous: convert numeric codes, such as assessment centre IDs, to `str` or `category` to treat them as categories. Collinear covariates (e.g. age and uncentered age²) stop PLINK with `VIF_TOO_HIGH`, which rescaling does not fix: center continuous covariates before squaring them or forming interactions.
- **Allele alignment**: each PLINK result is the contrast of the tested allele (A1) against the omitted allele, and it is aligned to `EA`. With `NEA`, the `EA`/`NEA` pair must be that contrast (or its complement for non-palindromic SNVs); at multi-allelic variants, a pair that PLINK did not test (e.g. C vs G when REF is A) is removed. Without `NEA`, only biallelic variants are kept, with `EA` one of their REF/ALT alleles.
- **Result**: `G.data` is replaced by the tested variants, with their cohort `BETA`, `SE`, `P` and a recomputed `FSTAT`. Variants not found in the genotype files, removed during allele alignment, or without a SNP ID are dropped (counts are printed).

## MR harmonization

MR workflows rely on aligning exposure and outcome effects to the same effect allele.

Implementation: {py:func}`genal.MR_tools.harmonize_MR`

High-level steps:

1. Merge exposure and outcome on `SNP` (alleles are compared in upper case). Rows without a valid SNP ID (missing, or a placeholder such as `.`) are excluded before matching, `BETA`/`SE`/`EAF` values stored as text are converted to numbers (values that cannot be converted are treated as missing), and `EAF` values outside [0, 1] are set to missing.
2. Align alleles:
   - If alleles are inverted, swap outcome alleles and flip the outcome effect (and outcome EAF when available).
   - If alleles are neither aligned nor inverted (and not palindromic), flip outcome alleles to their complements (A↔T, C↔G, base by base, so indels are complemented too) and retry.
3. Drop SNPs with remaining allele mismatches.
4. Handle palindromes based on `action`:
   - `action=1`: keep all palindromes (no attempt to resolve strand ambiguity)
   - `action=2`: drop “ambiguous” palindromes (exposure or outcome EAF missing or strictly between `eaf_threshold` and `1 - eaf_threshold`, as in TwoSampleMR) and flip the remaining ones when exposure/outcome frequencies are on opposite sides of 0.5. If either dataset has no `EAF` column, `action=3` is applied (with a warning).
   - `action=3`: drop all palindromes
5. Remove rows with missing or non-finite `BETA`/`SE`, or `SE ≤ 0`, and rows with an exposure effect `BETA_e = 0` (not associated with the exposure: their ratio estimate and their orientation on the exposure effect are undefined; `BETA_o = 0` is kept). Then, if the exposure or the outcome contains several rows for the same SNP ID, keep the first remaining allele-compatible row (a warning lists them), so that `MR`, `MR_loo` and `MRpresso` always use the same harmonized table.

### Proxy SNPs

With {py:meth}`genal.Geno.query_outcome` (`proxy=True`), an exposure SNP missing from the outcome is replaced by its best proxy: the outcome variant in highest LD ($r^2 \geq$ `r2`, assigned by decreasing $|r|$, one target per proxy) that is not itself an instrument found in the outcome, is non-palindromic, is allele-compatible with the reference panel and has valid `BETA`/`SE`. The proxy's effect is expressed for its reference-panel minor allele, its sign is flipped when $r < 0$ (no scaling by $r$), and it is reported for the missing SNP with that SNP's reference-panel minor allele as `EA` and, when the outcome has an `EAF` column, its panel frequency as `EAF` (used by the `action=2` palindrome checks).

Variants identified by coordinates are matched to the reference panel to search for proxies; an identifier given as `CHR:POS` (or `CHR:POS:A1:A2`) that matches several panel variants is not resolved (a message suggests giving `CHR:POS:A1:A2` or the rsID).

## Two-sample MR estimators

Implementation: MR wrappers in {py:func}`genal.MR_tools.MR_func`, estimators in `genal/MR.py`.

For each SNP $i$, define:

- exposure association: $\beta_{e,i}$ with standard error $se_{e,i}$
- outcome association: $\beta_{o,i}$ with standard error $se_{o,i}$

The per-variant ratio estimate is:

```{math}
\hat{\theta}_i = \frac{\beta_{o,i}}{\beta_{e,i}}
```

### Method map (what `methods=[...]` runs)

`genal` accepts the following method names in {py:meth}`genal.Geno.MR`:

| `methods` entry | Estimator | Minimum instruments | Notes |
|---|---:|---:|---|
| `"IVW"` | IVW (WLS) | 2 | Under-dispersion correction in the SE; returns Cochran’s Q |
| `"IVW-RE"` | IVW (WLS) | 2 | No under-dispersion correction; returns Cochran’s Q |
| `"IVW-FE"` | IVW (WLS) | 2 | Fixed-effects-style scaling; returns Cochran’s Q |
| `"UWR"` | Unweighted regression | 2 | OLS through origin; residual-based OLS SE; no Q |
| `"Egger"` | MR-Egger + intercept | 3 | WLS with intercept; returns both slope and intercept (+ Q) |
| `"Egger-boot"` | MR-Egger with bootstrap SE + intercept | 3 | MR-Egger estimates with bootstrap standard errors; returns slope + intercept |
| `"WM"` | Weighted median | 3 | Bootstrap SE |
| `"WM-pen"` | Penalised weighted median | 3 | Downweights outlier ratio estimates using `penk`; bootstrap SE |
| `"Simple-median"` | Simple median | 3 | Unweighted median of ratio estimates; bootstrap SE |
| `"Simple-mode"` | Simple mode | 3 | Mode of ratio estimates via KDE; bandwidth scaled by `phi` |
| `"Weighted-mode"` | Weighted mode | 3 | KDE weighted by inverse-variance-like weights; bandwidth scaled by `phi` |
| `"Sign"` | Sign concordance test | 6 | Binomial test of sign agreement between exposure and outcome effects |

`methods="all"` runs every method of the table. The `method` column of the results uses the labels of `genal.constants.MR_METHODS_NAMES` (e.g. `"IVW-RE"` → “Inverse Variance Weighted (Random Effects)”, `"Egger"` → “MR Egger” and “Egger Intercept”).

When exactly one instrument remains after harmonization, `MR()` returns the **Wald ratio** for that SNP (method “Wald ratio”, whatever `methods` contains): $\hat{\theta} = \beta_o/\beta_e$, $se = se_o/|\beta_e|$ (first-order, as in TwoSampleMR), with a two-sided normal p-value.

Confidence intervals (odds-ratio intervals with `odds=True`, forest and leave-one-out plots) use the same distribution as each method's p-value: a t-distribution with $nSNP-2$ degrees of freedom for the MR-Egger slope and intercept, with $nSNP-1$ degrees of freedom for the mode estimators, and the normal distribution for the other methods.

Bootstrap-based methods (`Egger-boot`, median and mode estimators) draw their bootstrap samples from a random generator seeded with the `seed` argument of {py:meth}`genal.Geno.MR` / {py:meth}`genal.Geno.MR_loo`: for a given seed, results are reproducible. The MR methods run on a single core (the `cpus` argument of `Geno.MR` and `Geno.MR_loo` is ignored). Bootstrap standard errors are standard deviations (`ddof=1`) of the bootstrap estimates, except for the mode estimators, whose standard error is 1.4826 × the median absolute deviation of the bootstrap estimates (as in TwoSampleMR).

### Common quantities used by ratio-based methods

For ratio-based estimators (median/mode), `genal` uses the delta-method standard error of the per-variant ratio:

```{math}
se(\hat{\theta}_i) =
\sqrt{
\frac{se_{o,i}^2}{\beta_{e,i}^2}
\;+\;
\frac{\beta_{o,i}^2 \, se_{e,i}^2}{\beta_{e,i}^4}
}
```

### IVW family (WLS through origin)

`genal` implements IVW as a weighted regression of $\beta_o$ on $\beta_e$ without intercept:

```{math}
\beta_{o,i} = \theta \, \beta_{e,i} + \varepsilon_i
```

with weights $w_i = 1/se_{o,i}^2$.

All IVW variants return heterogeneity diagnostics via Cochran’s Q:

```{math}
Q = \sum_i w_i \left(\beta_{o,i} - \hat{\theta}\beta_{e,i}\right)^2
```

Q is compared with a $\chi^2$ distribution with $J-1$ degrees of freedom ($J-2$ for MR-Egger, from the residuals of its regression with intercept); the `Q`, `Q_df` and `Q_pval` columns are returned with `heterogeneity=True`.

The three IVW variants share the WLS estimate and differ in the standard error. With $se_{FE} = 1/\sqrt{\sum_i w_i \beta_{e,i}^2}$ and $\sigma$ the residual standard error of the regression:
- `"IVW"`: $se = se_{FE} \cdot \max(1, \sigma)$ (multiplicative random effects with an **under-dispersion correction**, as TwoSampleMR `mr_ivw`).
- `"IVW-RE"`: $se = se_{FE} \cdot \sigma$ (raw WLS standard error, as TwoSampleMR `mr_ivw_mre`).
- `"IVW-FE"`: $se = se_{FE}$ (fixed effects).

These closed forms stay valid for perfect fits ($\sigma = 0$). MR-Egger uses the same scaling, $\max(1, \sigma)$, of its fixed-effect standard errors.

Practical note: ratio-based methods divide by $\beta_{e,i}$ and become unstable when it is close to 0. Harmonization removes missing or non-finite `BETA`/`SE`, `SE ≤ 0` and exposure effects equal to 0, but keeps exposure effects close to 0, so instruments should be strongly associated with the exposure.

### MR-Egger (WLS with intercept)

Egger is implemented as a weighted regression with intercept:

```{math}
\beta_{o,i} = \alpha + \theta \, \beta_{e,i} + \varepsilon_i
```

After harmonization, `genal` orients $\beta_e$ to be non-negative and applies the same sign to $\beta_o$. The intercept $\alpha$ is returned as “Egger Intercept”.

`"Egger-boot"` reports the MR-Egger slope and intercept (identical to `"Egger"`) with bootstrap standard errors: exposure/outcome effects are repeatedly sampled from normal distributions using their SEs, the same oriented, weighted MR-Egger regression is refitted in each replicate, and the SE is the standard deviation of the replicates. The p-value is two-sided, from the normal distribution. The bootstrap SE reflects the sampling error of the SNP effects only, not the between-variant heterogeneity included in the regular Egger SE, so it can be smaller than the `"Egger"` SE when Cochran's Q is large.

### Unweighted regression

`"UWR"` regresses $\beta_o$ on $\beta_e$ through the origin without weights. Its standard error is the residual-based OLS standard error. Q is not reported (`NaN`): the unweighted residual sum of squares is not a heterogeneity statistic.

### Median estimators

`genal` implements:

- **Weighted median**: weighted median of $\hat{\theta}_i$ (weights derived from the delta-method variance), with bootstrap SE.
- **Penalised weighted median**: downweights variants whose ratio estimates deviate strongly from the weighted-median estimate, then recomputes the weighted median (controlled by `penk`), with bootstrap SE.
- **Simple median**: unweighted median of $\hat{\theta}_i$, with bootstrap SE.

### Mode estimators

Mode methods estimate the mode of the distribution of $\hat{\theta}_i$ using kernel density estimation (KDE):

- **Simple mode**: unweighted KDE
- **Weighted mode**: KDE weighted by inverse-variance-like weights

The bandwidth uses a modified Silverman rule multiplied by the user-provided factor `phi`. The estimate is the global maximum of the exact (weighted) Gaussian KDE: the density is evaluated on a fine grid (spacing of a tenth of the bandwidth) covering the regions within 4 bandwidths of the ratio estimates, and every local maximum close to the highest one is refined. The grid density is only computed where an upper bound shows that it can come within 0.25% of its highest value, which gives the same estimate as the whole grid. If distinct maxima have the same density (up to rounding errors, e.g. with symmetric ratio estimates), the estimate is ambiguous: the smallest maximum is reported and a warning is printed (TwoSampleMR also reports the first of equal grid maxima). Narrow peaks, peaks located between ratio estimates, and data with a few extreme ratio estimates are handled. In ordinary data this agrees with TwoSampleMR (which evaluates R's `density` on a 512-point grid) to within one grid step. `genal` then uses bootstrap sampling (normal sampling with the per-variant ratio SE) to estimate uncertainty.

### Sign concordance test

The sign method tests whether exposure and outcome effects tend to have the same sign across variants. `genal` performs a binomial test against the null of 50% sign agreement. Zero effects are excluded; `b` is $2\,(x/n - 0.5)$ for $x$ concordant signs among $n$ variants, and `se` is `NaN`.

## Leave-one-out MR

Implementation: {py:func}`genal.MR_tools.MR_loo_func`, wrapped by {py:meth}`genal.Geno.MR_loo`.

Leave-one-out MR iterates over all instruments, sequentially removing each SNP and re-estimating the causal effect using the remaining instruments. This identifies variants that have a disproportionate influence on the overall estimate.

For each SNP $i$ in the instrument set:

1. Remove SNP $i$ from the harmonized data.
2. Re-run the selected MR method on the remaining $J-1$ SNPs.
3. Store the resulting estimate $\hat{\theta}_{-i}$.

The "influence" of SNP $i$ is defined as:

```{math}
\text{influence}_i = \left| \hat{\theta}_{-i} - \hat{\theta}_{\text{all}} \right|
```

where $\hat{\theta}_{\text{all}}$ is the estimate using all instruments. It is used by `MR_loo_plot(top_influential=True)` to select the SNPs shown; it is not a column of the `MR_loo()` table.

Notes:
- Any single MR method can be used (IVW, Egger, weighted median, mode-based, etc.).
- At least 3 instruments are required (so that each LOO subset has ≥2 instruments).

## MR-PRESSO (summary)

Implementation: {py:func}`genal.MRpresso.mr_presso`, wrapped by {py:meth}`genal.Geno.MRpresso`.

At a high level, `genal`'s MR-PRESSO implementation:

- computes an observed leave-one-out residual sum of squares (RSS),
- simulates an empirical RSS distribution by generating random data under the fitted model (`n_iterations` controls this),
- if `outlier_test=True` and the global p-value is significant (`global_test_p < significance_p`), performs:
  - an outlier test (per-variant p-values, Bonferroni-corrected),
  - an optional distortion test (whether the causal estimate changes materially after removing outliers).

(mr-presso-computation)=
### Computation

`genal` computes the same statistics and tests as the MR-PRESSO R package (`getRSS_LOO`, `getRandomData`,
`getRandomBias`), with simulated datasets drawn from the same distributions, and two differences in how they are
computed:

- the leave-one-out IVW estimates use their closed form,
  $\hat{\theta}_{-i} = \sum_{j \neq i} w_j x_j y_j \,/\, \sum_{j \neq i} w_j x_j^2$, instead of one regression per
  left-out variant (the R package refits $J$ regressions on the observed data for every simulated dataset);
- the simulated datasets are processed in vectorized chunks, so that memory use does not grow with `n_iterations`.

These changes do not alter the results, which match those of the MR-PRESSO R package.

As a result, `genal` takes over 99.7% less time than the R package, and runs on a single core (its `cpus` argument is
ignored). Benchmark on one CPU for both, Apple M1 Pro, synthetic data with 3% outliers, outlier and distortion tests
run:

| Instruments | Iterations | MR-PRESSO R package | `genal` |
|---|---|---|---|
| 20 | 10,000 | 112 s | 0.12 s |
| 200 | 1,000 | 120 s | 0.04 s |
| 812 | 1,000 | 521 s | 0.06 s |
| 812 | 10,000 | about 87 min (extrapolated from 1,000 iterations) | 0.34 s |

### Distortion test

The distortion test assesses whether the $k$ detected outliers change the causal estimate more than removing $k$ ordinary (non-outlier) variants would. `genal` follows the MR-PRESSO R implementation (`getRandomBias`):

1. **Observed distortion**: $D_\text{obs} = (\hat{\theta}_\text{all} - \hat{\theta}_\text{no outliers}) / |\hat{\theta}_\text{no outliers}|$, where $\hat{\theta}_\text{no outliers}$ uses the $J-k$ non-outliers. The reported `distortion_test_coefficient` is $100 \cdot D_\text{obs}$ (percent).
2. **Null distribution**: for each iteration, build a set of the same size $J-k$ made of the $k$ outliers and $J-2k$ non-outliers drawn with replacement, fit the IVW model on it to get $\hat{\theta}_\text{exp}$, and compute $D_\text{exp} = (\hat{\theta}_\text{all} - \hat{\theta}_\text{exp}) / |\hat{\theta}_\text{exp}|$. Because the outliers stay in the set and only non-outliers are left out, $D_\text{exp}$ measures the distortion obtained by removing ordinary variants.
3. **P-value**: $p = \#(|D_\text{exp}| > |D_\text{obs}|) / n_\text{iterations}$.

The test is skipped (coefficient and p-value set to `NaN`, with a message) when $k \geq J/2$, since $J-2k$ non-outliers cannot then be drawn, and when the outlier-corrected estimate is exactly 0, since the distortion relative to it is undefined. All other results are still returned.

### Output structure

`Geno.MRpresso()` returns four objects:
- `mod_table`: `Raw` and `Outlier-corrected` rows of the IVW regression through the origin (weights $1/se_{o,i}^2$). As in the MR-PRESSO R package, the SE is the raw WLS standard error $se_{FE} \cdot \sigma$ (as `"IVW-RE"`) and the p-value uses a t-distribution with $nSNP-1$ degrees of freedom,
- `GlobalTest`: RSS and global p-value,
- `OutlierTest`: per-variant outlier p-values (empty if the outlier test did not run); SNP IDs as row labels,
- `BiasTest`: distortion test result dictionary containing `"outliers_indices"` (SNP IDs, or a message when no or all SNPs are outliers), `"distortion_test_coefficient"`, and `"distortion_test_p"` (empty if `distortion_test=False` or the outlier test did not run).

MR-PRESSO requires at least 4 instruments after harmonization, and `n_iterations` must exceed their number; with fewer instruments, no MR-PRESSO results are stored, and those of earlier calls are cleared. When {py:func}`genal.MRpresso.mr_presso` is called directly, it also requires a unique row index and non-zero exposure effects. The random datasets are drawn from generators seeded with the `seed` argument, so results are reproducible for a given seed. Empirical p-values of 0 are reported as 0 with a printed message stating the bound (`< 1/n_iterations`, or `< nSNP/n_iterations` for the Bonferroni-corrected outlier test).

The SNP IDs flagged by the outlier test are stored in `G.MRpresso_outliers` (used to highlight the outliers on plots). Whenever at least one but not all instruments are flagged (whether or not the distortion test runs), `genal` fills the `Outlier-corrected` row of `mod_table` (IVW on the non-outliers) and stores the outlier-removed harmonized table in `G.MRpresso_subset_data`, which allows rerunning MR with `Geno.MR(use_mrpresso_data=True)`. If all instruments are flagged, the `Outlier-corrected` row is `NaN` and `use_mrpresso_data=True` raises an error. These attributes are cleared when `query_outcome()` is called with a new outcome.

## Colocalization (ABF)

Implementation: {py:func}`genal.colocalization.coloc_abf_func`.

For each trait and each variant, `genal` computes an approximate Bayes factor (ABF) from regression estimates. Let:

```{math}
z = \frac{\beta}{\sqrt{\mathrm{var}(\beta)}}
```

and let the prior standard deviation be:

- quantitative trait: $sd_\mathrm{prior} = 0.15 \cdot sd_Y$
- case-control trait: $sd_\mathrm{prior} = 0.2$

Then, with $r = \frac{sd_\mathrm{prior}^2}{sd_\mathrm{prior}^2 + \mathrm{var}(\beta)}$, the log-ABF is:

```{math}
\log(\mathrm{ABF}) = \frac{1}{2}\left(\log(1-r) + r z^2\right)
```

`genal` combines per-trait ABFs to compute posterior probabilities for hypotheses $H_0$–$H_4$ (no association, trait1 only, trait2 only, both traits but different variants, both traits and shared variant) using priors `p1`, `p2`, `p12` (defaults `1e-4`, `1e-4`, `1e-5`, as in `coloc.abf`). The $H_3$ term is computed directly as the sum over pairs of distinct variants, $\sum_{i \neq j} ABF_{1,i} \, ABF_{2,j}$, in a numerically stable way; this equals coloc's formula except where coloc's subtraction loses precision (e.g. one variant dominating both traits with a very small `p12`).

Each dataset is processed separately before merging, as in `coloc.abf`:

- Variants with non-finite `BETA`/`SE`, `SE ≤ 0`, or a missing merge key (missing `CHR/POS`, or a missing/placeholder SNP ID such as `.`) are removed, and duplicated variants (same merge key and, when `EA`/`NEA` are present, same allele pair) are reduced to the row with the largest `|BETA/SE|` (lowest P), with a warning. `coloc.abf` stops with an error on missing values and duplicated SNPs; `genal` removes these rows.
- For quantitative traits without `sdY1` (`sdY2`), $sd_Y$ is estimated from all the variants of the dataset (`sdY.est`, from `SE`, `EAF` strictly between 0 and 1, and the sample size `n1` (`n2`)). If it cannot be estimated (no sample size or usable `EAF`, or a non-finite estimate), $sd_Y = 1$ is assumed (effects on a standardized scale) and a warning is printed; a negative estimate raises an error.
- Priors are adjusted for very large regions: `p1` (`p2`) is set to `1/(m+1)` when `p1·m ≥ 1`, with `m` the number of variants of dataset 1 (2) after these checks; `p12` uses the number of overlapping variants.

When `EA`/`NEA` are available in both datasets, two rows are the same variant only if their allele pairs match (in either orientation, or on the complementary strand for non-palindromic SNVs); other pairs are removed (as are alleles not made only of A, C, G, T). Each variant of either dataset is used in a single pair (direct matches are preferred over strand-complement matches), so multi-allelic sites keep one pair per matched variant. Trait types must be `"quant"` or `"cc"`.

Practical note: colocalization is usually done on a **locus/region**, not genome-wide, and requires careful choice of priors and trait type (`quant` vs `cc`).

## Other utilities

### Liftover

Implementation: {py:meth}`genal.Geno.lift` (`genal/lift.py`).

- By default, positions are lifted in Python through the aligned block of the chain file that contains them, in the chain with the highest score (the mappings of `pyliftover`); `liftover_path` lifts with the UCSC `liftOver` executable instead.
- Positions are 1-based. Variants mapped to the opposite strand (inverted regions) get their `EA`/`NEA` reverse-complemented, and their number is printed.
- Variants that cannot be lifted, or that map to alternative, unplaced or patch contigs (`chrUn_…`, `…_alt`, `…_random`, `…_fix`), are removed. Chromosomes 23, 24 and 26 are lifted as chrX, chrY and chrM; variants on chromosome 25 (XY) and rows with missing `CHR`/`POS` are removed.

### Allele frequencies

Implementation: {py:meth}`genal.Geno.update_eaf`.

- PLINK `--freq` is run on the reference panel and the frequency is aligned to `EA`. A reference variant is used only if its alleles are {`EA`, `NEA`} (or contain `EA` when there is no `NEA` column); unmatched or ambiguous variants get `NaN` (or keep their `EAF` with `fill=True`).
- If fewer than 90% of the variants found in the panel (at their position, or by SNP ID) match its alleles, a genome build mismatch warning is printed. Variants absent from the panel are reported as an informational count.

### GWAS Catalog

Implementation: {py:meth}`genal.Geno.query_gwas_catalog` (`genal/snp_query.py`).

- Only canonical dbSNP rsIDs (`rs` + digits) are queried; all result pages of each SNP are requested.
- `max_associations` keeps the most significant trait-study associations per SNP after filtering and deduplication, ranked on the exact p-value from the catalog's mantissa/exponent (p-values too small for a float are ranked correctly).

### Gene filtering

Implementation: {py:meth}`genal.Geno.filter_by_gene` (`genal/genes.py`).

- The gene table (downloaded on first use into the `v2/` subfolder of the reference folder) has one row per Ensembl gene on chromosomes 1–22, X, Y and MT, with GRCh38 coordinates from Ensembl release 116 and GRCh37 coordinates from the Ensembl GRCh37 annotation (lifted from GRCh38 for genes absent from it; column `coords_37_source`).
- Genes on X, Y and MT match `CHR` values 23, 24 and 26. The window extends `window_size/2` on each side of the gene boundaries, and the `Distance` column gives the distance to the nearest boundary (negative before the gene, 0 inside, positive after).
