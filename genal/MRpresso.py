import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from typing import NamedTuple
from tqdm import tqdm
from numpy.random import SeedSequence, default_rng

##todo: implement the multivariable option, for the moment we assume only 1 BETA_e column

# Maximum number of (simulated dataset x SNP) values held in one array: the simulations are processed in chunks of
# iterations (arrays of 256 kB, which stay in the CPU caches), so the memory used does not grow with n_iterations
CHUNK_VALUES = 2**15
# Leave-one-out sums of a row are computed without subtracting from the total when one term exceeds DOMINANCE times
# the sum of the absolute values of the others (see loo_sums)
DOMINANCE = 100


# MR-PRESSO main function
def mr_presso(
    data,
    BETA_e_columns=["BETA_e"],
    n_iterations=1000,
    outlier_test=True,
    distortion_test=True,
    significance_p=0.05,
    cpus=None,
    seed=None,
):
    """
    Perform the MR-PRESSO algorithm for detection of horizontal pleiotropy.

    Args:
        data (pd.DataFrame): DataFrame with at least 4 columns: BETA_o (outcome), SE_o, BETA_e (exposure), SE_e.
        BETA_e_columns (list): List of exposure beta columns.
        n_iterations (int): Number of steps performed (random data generation).
        outlier_test (bool): If True, identifies outlier SNPs responsible for horizontal pleiotropy.
        distortion_test (bool): If True, tests significant distortion in the causal estimates.
        significance_p (float): Statistical significance threshold for the detection of horizontal pleiotropy.
        cpus (int, optional): Not used: MR-PRESSO runs on a single core. A message is printed when it is given.
        seed (int, optional): Random seed for the simulations (random data generation and distortion test).
            Default None (non-reproducible).

    Returns:
        mod_table (pd.DataFrame): DataFrame with the original and outlier-corrected inverse variance-weighted MR results.
            The outlier-corrected row is computed whenever the outlier test flagged at least one (but not all) SNPs.
        GlobalTest (dict): Dictionary with p-value of the global MR-PRESSO test.
        OutlierTest (pd.DataFrame): DataFrame with p-value for each SNP for the outlier test.
        BiasTest (dict): Dictionary with results of the distortion test.
        subset_data (pd.DataFrame or None): Data without the outliers whenever the outlier test flagged at least one
            (but not all) SNPs, None otherwise.

    Notes:
        Distortion test (as in the MR-PRESSO R package): the distortion is
        D = (b_all - b_no_outliers) / |b_no_outliers|, with b the IVW estimate (weighted regression through the origin,
        weights 1/SE_o^2); the reported distortion coefficient is 100 * D (percent). Its null distribution is built from
        n_iterations IVW fits on the k outliers plus J - 2k non-outliers drawn at random with replacement (J - k SNPs per
        fit), and p = proportion of |D_null| > |D|.
        The test is skipped if k >= J/2.
        Empirical p-values equal to 0 mean p < 1/n_iterations (global and distortion tests) or
        p < nSNP/n_iterations (Bonferroni-corrected outlier test); a message is printed in these cases.
        Computation (single process): the leave-one-out IVW estimates of the observed and simulated datasets use the
        closed form b_-i = sum_{j != i} w_j x_j y_j / sum_{j != i} w_j x_j^2, and the simulated datasets are processed
        in chunks (bounded memory). Simulated dataset i is drawn from child i of the seed's random stream.
    """
    # The index labels identify the SNPs (outliers are removed by label)
    if not data.index.is_unique:
        raise ValueError(
            "The index of data must be unique (one label per SNP): MR-PRESSO identifies and removes outliers by index label."
        )

    # Transforming the data
    required_cols = ["BETA_o", *BETA_e_columns, "SE_o", "SE_e"]
    keep_cols = required_cols.copy()
    if "SNP" in data.columns:
        keep_cols = ["SNP", *keep_cols]
    data = data[keep_cols].dropna(subset=required_cols)
    if not np.isfinite(data[required_cols].to_numpy(dtype=float)).all():
        raise ValueError(f"The effects and standard errors ({', '.join(required_cols)}) must be finite.")
    if (data[BETA_e_columns[0]] == 0).any():
        raise ValueError(
            "Some exposure effects are 0: the variants are oriented on the sign of the exposure effect, which is undefined "
            "for them (and their ratio estimate). Remove the variants with an exposure effect of 0."
        )
    data[["BETA_o", *BETA_e_columns]] = data[["BETA_o", *BETA_e_columns]].multiply(
        np.sign(data[BETA_e_columns[0]]), axis=0
    )
    data["Weights"] = 1 / (data["SE_o"] ** 2)

    if len(data) <= len(BETA_e_columns) + 2:
        raise ValueError(
            f"Not enough instrumental variables (variants): {len(data)} provided but at least {len(BETA_e_columns) + 3} are required."
        )
    if len(data) >= n_iterations:
        raise ValueError(
            "Not enough elements to compute empirical P-values, increase n_iterations"
        )

    if len(BETA_e_columns) != 1:
        raise ValueError("MR-PRESSO is implemented for a single exposure (one BETA_e column).")
    if (data["SE_e"] < 0).any() or (data["SE_o"] <= 0).any():
        raise ValueError(
            "Standard errors must be non-negative (SE_e) and positive (SE_o: the variants are weighted by 1 / SE_o^2)."
        )
    if cpus is not None:
        print("MR-PRESSO is now single-core only thanks to its optimization: the cpus argument is ignored.")

    # Independent random streams: one per simulated dataset and one for the distortion test
    seed_data, seed_distortion = np.random.SeedSequence(seed).spawn(2)

    beta_e = data[BETA_e_columns[0]].to_numpy(dtype=float)
    beta_o = data["BETA_o"].to_numpy(dtype=float)
    sqrt_w = np.sqrt(data["Weights"].to_numpy(dtype=float))

    print(f"Running the MR-PRESSO algorithm with N = {n_iterations} iterations.")
    # 1- Computing the observed residual sum of squares (RSS) and the observed leave-one-out IVW estimates
    print(f"Computing the observed residual sum of squares...")
    rss_obs, b_obs = loo_rss((beta_e * sqrt_w)[np.newaxis], (beta_o * sqrt_w)[np.newaxis])
    RSSobs, b_loo = rss_obs[0], b_obs[0]
    Dif = beta_o - beta_e * b_loo  # Observed residuals of the outlier test

    # 2- Computing the distribution of expected residual sum of squares (RSS)
    print("Computing the global MR-PRESSO p-value...")
    sim_inputs = SimulationInputs(
        beta_e=beta_e,
        se_e=data["SE_e"].to_numpy(dtype=float),
        se_o=data["SE_o"].to_numpy(dtype=float),
        sqrt_w=sqrt_w,
        b_loo=b_loo,
        abs_dif=np.abs(Dif),
    )
    RSSexp, n_larger = simulate_statistics(sim_inputs, n_iterations, seed_data)

    global_p = np.sum(RSSexp > RSSobs) / n_iterations
    if global_p == 0:
        print(
            f"Global test p-value < 1/n_iterations (= {1 / n_iterations:.2g}): none of the simulated RSS exceeded the observed RSS (reported as 0)."
        )

    GlobalTest = {"RSSobs": RSSobs, "global_test_p": global_p}

    # 3- Computing the single IV outlier test
    if global_p < significance_p and outlier_test:
        print("Global p-value is below the significance threshold. Running the Outlier test.")

        pvals = n_larger / n_iterations

        OutlierTest = pd.DataFrame({"RSSobs": Dif**2, "Pvalue": pvals})

        OutlierTest.index = data.index
        OutlierTest["Pvalue"] = np.minimum(
            OutlierTest["Pvalue"] * len(data), 1
        )  # Bonferroni correction
        n_zero = int((OutlierTest["Pvalue"] == 0).sum())
        if n_zero > 0:
            print(
                f"Outlier test: {n_zero} SNPs have a p-value < nSNP/n_iterations (= {len(data) / n_iterations:.2g}) after Bonferroni correction (reported as 0)."
            )
        if data.shape[0] / n_iterations > significance_p:
            print(
                f"Warning: the Outlier test is unstable. The {significance_p} significance threshold cannot be obtained with {n_iterations} distributions. Increase n_iterations."
            )

    else:
        outlier_test = False
        OutlierTest = pd.DataFrame()

    # 4- Computing the test of the distortion of the causal estimate
    formula = f"BETA_o ~ -1 + {' + '.join(BETA_e_columns)}"
    mod_all = smf.wls(formula, data=data, weights=data["Weights"]).fit()

    BiasTest = {}
    subset_data = None
    mod_no_outliers = None

    if outlier_test:
        ref_outlier = OutlierTest.loc[OutlierTest["Pvalue"] <= significance_p].index

        k = len(ref_outlier)
        J = len(data)

        # Data without the outliers (used by MR methods with use_mrpresso_data=True) and outlier-corrected IVW fit
        if k == 0:
            print("No significant outliers found.")
        elif k == J:
            print("All SNPs considered as outliers.")
        else:
            print(f"{k}/{J} ({k/J*100:.2f}%) outliers found.")
            subset_data = data.drop(index=ref_outlier)
            mod_no_outliers = smf.wls(
                formula, data=subset_data, weights=subset_data["Weights"]
            ).fit()

    if distortion_test and outlier_test:
        if k == 0:
            print("Skipping the Distortion test (no significant outliers).")
            BiasTest = {
                "outliers_indices": "No significant outliers",
                "distortion_test_coefficient": np.nan,
                "distortion_test_p": np.nan,
            }
        elif k == J:
            print("Skipping the Distortion test (all SNPs considered as outliers).")
            BiasTest = {
                "outliers_indices": "All SNPs considered as outliers",
                "distortion_test_coefficient": np.nan,
                "distortion_test_p": np.nan,
            }
        elif J - 2 * k <= 0:
            print(
                f"Skipping the Distortion test: {k} of the {J} SNPs are outliers, but its null distribution requires fewer outliers than non-outliers (J - 2k > 0)."
            )
            BiasTest = {
                "outliers_indices": list(ref_outlier),
                "distortion_test_coefficient": np.nan,
                "distortion_test_p": np.nan,
            }
        else:
            if len(BETA_e_columns) != 1:
                raise ValueError("Distortion test not done for multi MR.")

            print("Running the Distortion test.")
            main_col = BETA_e_columns[0]
            beta_all = float(mod_all.params[main_col])
            beta_no_outliers = float(mod_no_outliers.params[main_col])
            if beta_no_outliers == 0:
                # The distortion coefficient is relative to the outlier-corrected estimate
                print(
                    "Distortion test not computed: the outlier-corrected estimate is 0, so the distortion coefficient (relative to it) is undefined."
                )
                BiasTest = {
                    "outliers_indices": list(ref_outlier),
                    "distortion_test_coefficient": np.nan,
                    "distortion_test_p": np.nan,
                }
            else:
                D_obs = (beta_all - beta_no_outliers) / abs(beta_no_outliers)

                # Null distribution (getRandomBias of MR-PRESSO): IVW fit (through the origin, weights 1/SE_o^2) on the k outliers
                # plus J - 2k non-outliers drawn at random with replacement (subsets of J - k SNPs that include all the outliers)
                is_outlier = data.index.isin(ref_outlier)
                wxy = (data["Weights"] * data[main_col] * data["BETA_o"]).to_numpy()
                wxx = (data["Weights"] * data[main_col] ** 2).to_numpy()
                wxy_outliers, wxx_outliers = wxy[is_outlier].sum(), wxx[is_outlier].sum()
                non_outlier_pos = np.flatnonzero(~is_outlier)
                rng_distortion = default_rng(seed_distortion)
                betas_exp = np.empty(n_iterations)
                rows = max(1, CHUNK_VALUES // (J - 2 * k))
                for start in range(0, n_iterations, rows):
                    stop = min(start + rows, n_iterations)
                    sampled = draw_distortion_samples(rng_distortion, non_outlier_pos, stop - start, J - 2 * k)
                    betas_exp[start:stop] = ivw_with_samples(sampled, wxy, wxx, wxy_outliers, wxx_outliers)

                with np.errstate(divide="ignore", invalid="ignore"):
                    D_exp = (beta_all - betas_exp) / np.abs(betas_exp)
                p_value = float(np.sum(np.abs(D_exp) > abs(D_obs)) / n_iterations)
                if p_value == 0:
                    print(
                        f"Distortion test p-value < 1/n_iterations (= {1 / n_iterations:.2g}) (reported as 0)."
                    )

                BiasTest = {
                    "outliers_indices": list(ref_outlier),
                    "distortion_test_coefficient": float(100 * D_obs),
                    "distortion_test_p": p_value,
                }

    # 5- Format
    main_beta_col = BETA_e_columns[0]
    row_original = {
        "exposure": main_beta_col,
        "method": "Raw",
        "nSNP": len(data),
        "b": mod_all.params[main_beta_col],
        "se": mod_all.bse[main_beta_col],
        "pval": mod_all.pvalues[main_beta_col],
    }
    if mod_no_outliers is not None:
        row_corrected = {
            "exposure": main_beta_col,
            "method": "Outlier-corrected",
            "nSNP": len(data) - len(ref_outlier),
            "b": mod_no_outliers.params[main_beta_col],
            "se": mod_no_outliers.bse[main_beta_col],
            "pval": mod_no_outliers.pvalues[main_beta_col],
        }
    else:
        row_corrected = {
            "exposure": main_beta_col,
            "method": "Outlier-corrected",
            "nSNP": np.nan,
            "b": np.nan,
            "se": np.nan,
            "pval": np.nan,
        }

    mod_table = pd.DataFrame([row_original, row_corrected])

    return mod_table, GlobalTest, OutlierTest, BiasTest, subset_data


## MR-PRESSO helper functions
class SimulationInputs(NamedTuple):
    """Observed per-SNP arrays used to simulate the random datasets of the global and outlier tests."""

    beta_e: np.ndarray  # Exposure effects (oriented so that BETA_e >= 0)
    se_e: np.ndarray
    se_o: np.ndarray
    sqrt_w: np.ndarray  # Square root of the IVW weights 1/SE_o^2
    b_loo: np.ndarray  # Observed leave-one-out IVW estimates
    abs_dif: np.ndarray  # Absolute observed residuals of the outlier test


def loo_sums(t):
    """
    Sums of each row of t leaving out each element in turn (sum over j != i of t_j, for every i): the row total minus t_i.
    In rows where one term exceeds DOMINANCE times the sum of the absolute values of the others, subtracting it from the
    total would lose precision: there, the sums of the terms before and after i are added instead.
    """
    out = t.sum(axis=1, keepdims=True) - t
    abs_t = np.abs(t)
    largest = abs_t.max(axis=1)
    dominated = largest > DOMINANCE * (abs_t.sum(axis=1) - largest)
    if dominated.any():
        rows = t[dominated]
        exact = np.zeros_like(rows)
        np.cumsum(rows[:, :-1], axis=1, out=exact[:, 1:])
        exact[:, :-1] += np.cumsum(rows[:, :0:-1], axis=1)[:, ::-1]
        out[dominated] = exact
    return out


def loo_rss(X, Y):
    """
    Leave-one-out IVW estimates and residual sum of squares (getRSS_LOO of MR-PRESSO) of weighted data
    X = BETA_e * sqrt(w), Y = BETA_o * sqrt(w), with one dataset per row (2-D arrays, columns: SNPs).
    The estimate without SNP i is the regression through the origin of Y on X on the other SNPs,
    b_-i = sum_{j != i} X_j Y_j / sum_{j != i} X_j^2, and RSS = sum_i (Y_i - b_-i X_i)^2 (NaN terms ignored).

    Returns:
        tuple: RSS (one value per dataset) and the leave-one-out estimates (same shape as X).
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        b = loo_sums(X * Y) / loo_sums(X * X)
        rss = np.nansum((Y - b * X) ** 2, axis=1)
    return rss, b


def child_seed(seed_seq, i):
    """Child i of seed_seq.spawn(n) (any n > i) for a SeedSequence that has not spawned children, built directly."""
    return SeedSequence(seed_seq.entropy, spawn_key=seed_seq.spawn_key + (i,), pool_size=seed_seq.pool_size)


def simulate_data(inputs, seed_seq, start, stop):
    """
    Random datasets start, ..., stop - 1 (getRandomData of MR-PRESSO). Dataset i is drawn from its own random generator
    (child i of seed_seq): exposure effects ~ N(BETA_e, SE_e) and outcome effects ~ N(b_-i * BETA_e_i, SE_o_i), with b_-i
    the observed leave-one-out IVW estimate.

    The standard normal draws are those of rng.normal(BETA_e, SE_e) followed by rng.normal(b_-i * BETA_e_i, SE_o_i).

    Returns:
        tuple: exposure and outcome effects, arrays of shape (stop - start, nSNP).
    """
    z = np.empty((stop - start, 2, inputs.beta_e.size))
    for row, i in enumerate(range(start, stop)):
        default_rng(child_seed(seed_seq, i)).standard_normal(out=z[row])
    random_e = inputs.beta_e + inputs.se_e * z[:, 0]
    random_o = inputs.b_loo * inputs.beta_e + inputs.se_o * z[:, 1]
    return random_e, random_o


def random_data_statistics(inputs, random_e, random_o):
    """
    Statistics of random datasets (one per row): the RSS of the global test (with the weights of the observed data)
    and, for each SNP, the number of datasets whose outlier-test residual |BETA_o - BETA_e * b_-i| (b_-i observed) is
    larger than the observed one.
    """
    rss, _ = loo_rss(random_e * inputs.sqrt_w, random_o * inputs.sqrt_w)
    n_larger = np.count_nonzero(np.abs(random_o - random_e * inputs.b_loo) > inputs.abs_dif, axis=0)
    return rss, n_larger


def simulate_statistics(inputs, n_iterations, seed_seq):
    """
    RSS of n_iterations random datasets and, for each SNP, the number of datasets with a larger outlier-test residual
    than observed. The datasets are simulated in chunks of at most CHUNK_VALUES values; dataset i always uses random
    generator i, so the results do not depend on the chunks.
    """
    J = inputs.beta_e.size
    rows = max(1, CHUNK_VALUES // J)
    rss = np.empty(n_iterations)
    n_larger = np.zeros(J, dtype=np.int64)
    with tqdm(total=n_iterations, desc="Generating random data", ncols=100) as progress:
        for start in range(0, n_iterations, rows):
            stop = min(start + rows, n_iterations)
            random_e, random_o = simulate_data(inputs, seed_seq, start, stop)
            rss[start:stop], chunk_larger = random_data_statistics(inputs, random_e, random_o)
            n_larger += chunk_larger
            progress.update(stop - start)
    return rss, n_larger


def draw_distortion_samples(rng, positions, n_rows, size):
    """
    Non-outlier positions drawn with replacement for n_rows null fits of the distortion test (size per fit). One call
    with n_rows rows draws the same positions as n_rows successive calls with one row.
    """
    return rng.choice(positions, size=(n_rows, size), replace=True)


def ivw_with_samples(sampled, wxy, wxx, wxy_fixed, wxx_fixed):
    """
    IVW estimates through the origin, one per row of sampled: on the SNPs whose sums of w*x*y and w*x^2 are wxy_fixed
    and wxx_fixed plus the SNPs at the positions of the row.
    """
    return (wxy_fixed + wxy[sampled].sum(axis=1)) / (wxx_fixed + wxx[sampled].sum(axis=1))
