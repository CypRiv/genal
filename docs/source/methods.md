# Methods (implementation notes)

This page documents the main statistical models implemented in `genal`. It is intentionally brief and focuses on what is implemented in the codebase.

## Per-variant F-statistic (FSTAT)

Implementation: {py:func}`genal.geno_tools.fill_fstatistic`

`genal` computes a per-variant F-statistic (`FSTAT`) during preprocessing and after association testing. This statistic measures instrument strength for each variant.

- **Primary:** `FSTAT = (BETA / SE)²` when `BETA` and `SE` are available and `SE > 0`.
- **Fallback:** `FSTAT = χ²_isf(P, df=1)` (equivalent to `Z²` for a two-sided p-value) when `BETA/SE` are unavailable but `P` is present.

## MR harmonization

MR workflows rely on aligning exposure and outcome effects to the same effect allele.

Implementation: {py:func}`genal.MR_tools.harmonize_MR`

High-level steps:

1. Merge exposure and outcome on `SNP`.
2. Align alleles:
   - If alleles are inverted, swap outcome alleles and flip the outcome effect (and outcome EAF when available).
   - If alleles are neither aligned nor inverted (and not palindromic), flip outcome alleles to their complements (A↔T, C↔G) and retry.
3. Drop SNPs with remaining allele mismatches.
4. Handle palindromes based on `action`:
   - `action=1`: keep all palindromes (no attempt to resolve strand ambiguity)
   - `action=2`: use allele frequencies to drop “ambiguous” palindromes and flip the remaining ones when exposure/outcome frequencies disagree around 0.5
   - `action=3`: drop all palindromes

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
| `"UWR"` | Unweighted regression | 2 | OLS through origin; SE uses under-dispersion correction |
| `"Egger"` | MR-Egger + intercept | 3 | WLS with intercept; returns both slope and intercept (+ Q) |
| `"Egger-boot"` | MR-Egger bootstrap + intercept | 3 | Bootstrap SE/p-values (normal sampling); returns slope + intercept |
| `"WM"` | Weighted median | 3 | Bootstrap SE |
| `"WM-pen"` | Penalised weighted median | 3 | Downweights outlier ratio estimates using `penk`; bootstrap SE |
| `"Simple-median"` | Simple median | 3 | Unweighted median of ratio estimates; bootstrap SE |
| `"Simple-mode"` | Simple mode | 3 | Mode of ratio estimates via KDE; bandwidth scaled by `phi` |
| `"Weighted-mode"` | Weighted mode | 3 | KDE weighted by inverse-variance-like weights; bandwidth scaled by `phi` |
| `"Sign"` | Sign concordance test | 6 | Binomial test of sign agreement between exposure and outcome effects |

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

Differences between the three IVW variants in `genal`:
- `"IVW"`: WLS estimate, with an **under-dispersion correction** applied to the standard error.
- `"IVW-RE"`: WLS estimate, with the *raw* WLS standard error (no under-dispersion correction).
- `"IVW-FE"`: WLS estimate, with an alternative residual-variance scaling of the standard error.

Practical note: if $\beta_{e,i}$ is near 0 for some SNPs, ratio-based methods become unstable; `genal` filters invalid rows after harmonization.

### MR-Egger (WLS with intercept)

Egger is implemented as a weighted regression with intercept:

```{math}
\beta_{o,i} = \alpha + \theta \, \beta_{e,i} + \varepsilon_i
```

After harmonization, `genal` orients $\beta_e$ to be non-negative and applies the same sign to $\beta_o$. The intercept $\alpha$ is returned as “Egger Intercept”.

`"Egger-boot"` returns bootstrap SE/p-values by repeatedly sampling exposure/outcome effects from normal distributions using their SEs and refitting the regression.

### Median estimators

`genal` implements:

- **Weighted median**: weighted median of $\hat{\theta}_i$ (weights derived from the delta-method variance), with bootstrap SE.
- **Penalised weighted median**: downweights variants whose ratio estimates deviate strongly from the weighted-median estimate, then recomputes the weighted median (controlled by `penk`), with bootstrap SE.
- **Simple median**: unweighted median of $\hat{\theta}_i$, with bootstrap SE.

### Mode estimators

Mode methods estimate the mode of the distribution of $\hat{\theta}_i$ using kernel density estimation (KDE):

- **Simple mode**: unweighted KDE
- **Weighted mode**: KDE weighted by inverse-variance-like weights

The bandwidth uses a modified Silverman rule multiplied by the user-provided factor `phi`. `genal` then uses bootstrap sampling (normal sampling with the per-variant ratio SE) to estimate uncertainty.

### Sign concordance test

The sign method tests whether exposure and outcome effects tend to have the same sign across variants. `genal` performs a binomial test against the null of 50% sign agreement.

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

where $\hat{\theta}_{\text{all}}$ is the estimate using all instruments.

Notes:
- Any single MR method can be used (IVW, Egger, weighted median, mode-based, etc.).
- At least 3 instruments are required (so that each LOO subset has ≥2 instruments).

## MR-PRESSO (summary)

Implementation: {py:func}`genal.MRpresso.mr_presso`, wrapped by {py:meth}`genal.Geno.MRpresso`.

At a high level, `genal`'s MR-PRESSO implementation:

- computes an observed leave-one-out residual sum of squares (RSS),
- simulates an empirical RSS distribution by generating random data under the fitted model (`n_iterations` controls this),
- if the global p-value is significant (`global_test_p < significance_p`), performs:
  - an outlier test (per-variant p-values, Bonferroni-corrected),
  - an optional distortion test (whether the causal estimate changes materially after removing outliers).

### Distortion test

The distortion test assesses whether detected outliers materially bias the causal estimate. `genal` implements the following version:

1. **Observed distortion**: $D_\text{obs} = (\hat{\theta}_\text{all} - \hat{\theta}_\text{no outliers}) / |\hat{\theta}_\text{no outliers}|$
2. **Expected distortion null**: bootstrap resampling is performed *exclusively on the non-outlier subset*. For each iteration:
   - Sample with replacement $J-k$ SNPs from the non-outlier data (where $k$ is the number of detected outliers).
   - Fit the IVW model on the sampled data and record $\hat{\theta}_\text{exp}$.
   - Compute $D_\text{exp} = (\hat{\theta}_\text{all} - \hat{\theta}_\text{exp}) / |\hat{\theta}_\text{exp}|$.
3. **P-value**: $p = \text{mean}(|D_\text{exp}| > |D_\text{obs}|)$.

This differs from the original MR-PRESSO R implementation, which in some cases samples from the full dataset (including outliers) for the expected-bias regressions and was inconsistent with the paper's description.

### Output structure

`Geno.MRpresso()` returns four objects:
- `mod_table`: a small results table (`Raw` and `Outlier-corrected` rows; IVW model),
- `GlobalTest`: RSS and global p-value,
- `OutlierTest`: per-variant outlier p-values (empty if the global test is not significant); SNP IDs as row labels,
- `BiasTest`: distortion test result dictionary containing `"outliers_indices"` (SNP IDs), `"distortion_test_coefficient"`, and `"distortion_test_p"` (empty if distortion test was not run).

If outliers are found, `genal` stores the outlier-removed harmonized table and allows rerunning MR with `Geno.MR(use_mrpresso_data=True)`.

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

`genal` combines per-trait ABFs to compute posterior probabilities for hypotheses $H_0$–$H_4$ (no association, trait1 only, trait2 only, both traits but different variants, both traits and shared variant) using priors `p1`, `p2`, `p12`.

Practical note: colocalization is usually done on a **locus/region**, not genome-wide, and requires careful choice of priors and trait type (`quant` vs `cc`).
