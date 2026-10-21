import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from concurrent.futures import ProcessPoolExecutor
from sklearn.linear_model import LinearRegression
from tqdm import tqdm
from numpy.random import default_rng
from functools import partial

##todo: implement the multivariable option, for the moment we assume only 1 BETA_e column


# MR-PRESSO main function
def mr_presso(
    data,
    BETA_e_columns=["BETA_e"],
    n_iterations=1000,
    outlier_test=True,
    distortion_test=True,
    significance_p=0.05,
    cpus=5,
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
        cpus (int): Number of CPUs to use for parallel processing.
        seed (int, optional): Random seed for the simulations (random data generation and distortion test).
            Default None (non-reproducible). Each simulated dataset has its own random generator derived from the seed,
            so the results do not depend on the number of CPUs.

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

    # Independent random streams: one per simulated dataset and one for the distortion test
    seed_data, seed_distortion = np.random.SeedSequence(seed).spawn(2)

    print(f"Running the MR-PRESSO algorithm with N = {n_iterations} iterations.")
    # 1- Computing the observed residual sum of squares (RSS)
    print(f"Computing the observed residual sum of squares...")
    RSSobs = getRSS_LOO(data, BETA_e_columns, outlier_test)

    # 2- Computing the distribution of expected residual sum of squares (RSS)
    print("Computing the global MR-PRESSO p-value...")
    partial_parallel_RSS_LOO = partial(
        parallel_RSS_LOO, data=data, BETA_e_columns=BETA_e_columns
    )  # Wrapper function freezing the parallel_RSS_LOO call
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        results = list(
            tqdm(
                executor.map(partial_parallel_RSS_LOO, seed_data.spawn(n_iterations)),
                total=n_iterations,
                desc="Generating random data",
                ncols=100,
            )
        )

    RSSexp = [res[0] for res in results]
    Random_data_e = np.vstack([r[1] for r in results])
    Random_data_o = np.vstack([r[2] for r in results])

    rss_value = RSSobs[0] if outlier_test else RSSobs
    global_p = np.sum([r > rss_value for r in RSSexp]) / n_iterations
    if global_p == 0:
        print(
            f"Global test p-value < 1/n_iterations (= {1 / n_iterations:.2g}): none of the simulated RSS exceeded the observed RSS (reported as 0)."
        )

    if outlier_test:
        GlobalTest = {"RSSobs": RSSobs[0], "global_test_p": global_p}
    else:
        GlobalTest = {"RSSobs": RSSobs, "global_test_p": global_p}

    # 3- Computing the single IV outlier test
    if global_p < significance_p and outlier_test:
        print("Global p-value is below the significance threshold. Running the Outlier test.")

        if len(BETA_e_columns) == 1:
            Dif = data["BETA_o"].values - data["BETA_e"].values * RSSobs[1]
            Exp = Random_data_o - (Random_data_e * RSSobs[1])
        else:
            raise ValueError("Outlier test not done for multi MR.")

        abs_diffs = np.abs(Exp.T) > np.abs(Dif)[:, np.newaxis]
        pvals = np.sum(abs_diffs, axis=1) / Exp.shape[0]

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
                for i in range(n_iterations):
                    sampled = rng_distortion.choice(non_outlier_pos, size=J - 2 * k, replace=True)
                    betas_exp[i] = (wxy_outliers + wxy[sampled].sum()) / (wxx_outliers + wxx[sampled].sum())

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
# Define the matrix power operator
def power_eigen(x, n):
    """Matrix power of a symmetric matrix through its eigendecomposition (orthonormal eigenvectors)."""
    values, vectors = np.linalg.eigh(x)
    return vectors.dot(np.diag(values**n)).dot(vectors.T)


# Function to compute the residual sum of squares in a LOO framework
def getRSS_LOO(data, BETA_e_columns, returnIV):
    dataW = data[["BETA_o"] + BETA_e_columns].multiply(np.sqrt(data["Weights"]), axis=0)
    X = dataW[BETA_e_columns].values
    Y = dataW["BETA_o"].values

    # Matrix operations after LOO
    def loo_calculation(i):
        X_loo = np.delete(X, i, axis=0)
        Y_loo = np.delete(Y, i, axis=0)
        return power_eigen(X_loo.T.dot(X_loo), -1).dot(X_loo.T).dot(Y_loo)

    CausalEstimate_LOO = np.array([loo_calculation(i) for i in range(len(dataW))])

    if len(BETA_e_columns) == 1:
        CausalEstimate_LOO = CausalEstimate_LOO.reshape(-1)
        RSS = np.nansum((Y - CausalEstimate_LOO * X.reshape(-1)) ** 2)
    else:
        raise ValueError("Needs to do the getRSS_LOO for multi exposure.")
        # RSS = np.nansum((Y - np.sum(CausalEstimate_LOO.T * X, axis=1)) ** 2)

    if returnIV:
        return (RSS, CausalEstimate_LOO)
    return RSS


# Generate random data based on normal distributions
def getRandomData(data, BETA_e_columns=["BETA_e"], seed=None):
    rng = default_rng(seed)

    models = []
    for i in range(len(data)):
        lm = LinearRegression(fit_intercept=False)
        mask = np.ones(len(data), dtype=bool)
        mask[i] = False
        data_i = data.iloc[mask]
        lm.fit(
            data_i[BETA_e_columns], data_i["BETA_o"], sample_weight=data_i["Weights"]
        )
        models.append(lm)

    random_data_dict = {}
    for col, sd_col in zip(BETA_e_columns, ["SE_e"]):
        random_data_dict[col] = rng.normal(data[col], data[sd_col])

    random_data_dict["BETA_o"] = [
        rng.normal(
            model.predict(data.iloc[[i]][BETA_e_columns]), data.iloc[i]["SE_o"]
        ).item()
        for i, model in enumerate(models)
    ]
    random_data_dict["Weights"] = data["Weights"].values

    random_data_df = pd.DataFrame(random_data_dict)
    return random_data_df


# Function for the parallel executor in step 2: generate random data (from its own seed) and compute the expected residual sum of squares
def parallel_RSS_LOO(seed, data, BETA_e_columns):
    random_data = getRandomData(data, BETA_e_columns, seed)

    rss_exp = getRSS_LOO(random_data, BETA_e_columns, False)
    return (rss_exp, random_data["BETA_e"].values, random_data["BETA_o"].values)
