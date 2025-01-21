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

    Returns:
        mod_table (pd.DataFrame): DataFrame with the original and outlier-corrected inverse variance-weighted MR results.
        GlobalTest (dict): Dictionary with p-value of the global MR-PRESSO test.
        OutlierTest (pd.DataFrame): DataFrame with p-value for each SNP for the outlier test.
        BiasTest (dict): Dictionary with results of the distortion test.
    """
    # Transforming the data
    data = data[["BETA_o", *BETA_e_columns, "SE_o", "SE_e"]].dropna()
    data[["BETA_o", *BETA_e_columns]] = data[["BETA_o", *BETA_e_columns]].multiply(
        np.sign(data[BETA_e_columns[0]]), axis=0
    )
    data["Weights"] = 1 / (data["SE_o"] ** 2)

    if len(data) <= len(BETA_e_columns) + 2:
        raise Exception("Not enough instrumental variables (variants)")
    if len(data) >= n_iterations:
        raise Exception(
            "Not enough elements to compute empirical P-values, increase n_iterations"
        )

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
                executor.map(partial_parallel_RSS_LOO, range(n_iterations)),
                total=n_iterations,
                desc="Generating random data",
                ncols=100,
            )
        )

    RSSexp = [res[0] for res in results]
    Random_data_e = np.vstack([r[1] for r in results])
    Random_data_o = np.vstack([r[2] for r in results])

    global_p = np.sum([r > RSSobs[0] for r in RSSexp]) / n_iterations
    
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
        if data.shape[0] / n_iterations > significance_p:
            print(
                f"Warning: the Outlier test in unstable. The {significance_p} significance threshold cannot be obtained with {n_iterations} Distributions. Increase n_iterations."
            )

    else:
        outlier_test = False
        OutlierTest = pd.DataFrame()

    # 4- Computing the test of the distortion of the causal estimate
    formula = f"BETA_o ~ -1 + {' + '.join(BETA_e_columns)}"
    mod_all = smf.wls(formula, data=data, weights=data["Weights"]).fit()

    BiasTest = {}
    subset_data = None

    if distortion_test and outlier_test:
        ## Is there an error in the MRPRESSO code? The outlier indices are supposed to be excluded from the expected bias computation (as per the paper).
        def get_random_bias(BETA_e_columns, data, ref_outlier):
            indices = np.concatenate(
                [
                    ref_outlier,
                    np.random.choice(
                        list(set(range(len(data))) - set(ref_outlier)),
                        len(data) - len(ref_outlier),
                    ),
                ]
            )
            subset_data = data.iloc[indices[: -len(ref_outlier)]]
            mod_random = smf.wls(
                f"BETA_o ~ -1 + {' + '.join(BETA_e_columns)}",
                data=subset_data,
                weights=subset_data["Weights"],
            ).fit()
            return mod_random.params[BETA_e_columns]

        ref_outlier = OutlierTest.loc[OutlierTest["Pvalue"] <= significance_p].index

        if len(ref_outlier) > 0:
            if len(ref_outlier) < len(data):
                print(f"{len(ref_outlier)}/{len(data)} ({len(ref_outlier)/len(data)*100:.2f}%) outliers found. Running the Distortion test.")
                BiasExp = [
                    get_random_bias(BETA_e_columns, data, ref_outlier)
                    for _ in range(n_iterations)
                ]
                BiasExp = pd.concat(BiasExp, axis=1).transpose()

                subset_data = data.drop(ref_outlier)
                mod_no_outliers = smf.wls(
                    f"BETA_o ~ -1 + {' + '.join(BETA_e_columns)}",
                    data=subset_data,
                    weights=subset_data["Weights"],
                ).fit()

                BiasObs = (
                    mod_all.params[BETA_e_columns]
                    - mod_no_outliers.params[BETA_e_columns]
                ) / abs(mod_no_outliers.params[BETA_e_columns])
                BiasExp = (mod_all.params[BETA_e_columns] - BiasExp) / abs(BiasExp)

                p_value = np.sum(np.abs(BiasExp) > np.abs(BiasObs)) / n_iterations

                BiasTest = {
                    "outliers_indices": list(ref_outlier),
                    "distortion_test_coefficient": 100 * BiasObs.values[0],
                    "distortion_test_p": p_value.iloc[0],
                }
            else:
                print("All SNPs considered as outliers. Skipping the Distortion test.")
                BiasTest = {
                    "outliers_indices": "All SNPs considered as outliers",
                    "distortion_test_coefficient": np.nan,
                    "distortion_test_p": np.nan,
                }
        else:
            print("No significant outliers found. Skipping the Distortion test.")
            BiasTest = {
                "outliers_indices": "No significant outliers",
                "distortion_test_coefficient": np.nan,
                "distortion_test_p": np.nan,
            }

    # 5- Format
    row_original = {
        "exposure": BETA_e_columns[0],
        "method": "Raw",
        "nSNP": len(data),
        "b": mod_all.params["BETA_e"],
        "se": mod_all.bse["BETA_e"],
        "pval": mod_all.pvalues["BETA_e"],
    }
    if "mod_no_outliers" in locals():
        row_corrected = {
            "exposure": BETA_e_columns[0],
            "method": "Outlier-corrected",
            "nSNP": len(data) - len(ref_outlier),
            "b": mod_no_outliers.params["BETA_e"],
            "se": mod_no_outliers.bse["BETA_e"],
            "pval": mod_no_outliers.pvalues["BETA_e"],
        }
    else:
        row_corrected = {
            "exposure": BETA_e_columns[0],
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
    values, vectors = np.linalg.eig(x)
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
def getRandomData(data, BETA_e_columns=["BETA_e"]):
    rng = default_rng()

    models = []
    for i in range(len(data)):
        lm = LinearRegression(fit_intercept=False)
        data_i = data.drop(i)
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


# Function for the parallel executor in step 2: generate random data and compute the expected residual sum of squares
def parallel_RSS_LOO(i, data, BETA_e_columns):
    random_data = getRandomData(data, BETA_e_columns)

    rss_exp = getRSS_LOO(random_data, BETA_e_columns, False)
    return (rss_exp, random_data["BETA_e"].values, random_data["BETA_o"].values)
