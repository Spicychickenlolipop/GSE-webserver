"""Validation and statistics shared by the three processed-matrix pipelines."""

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def prepare_processed(matrix):
    matrix = matrix.replace([np.inf, -np.inf], np.nan).dropna(axis=1, how="all")
    if matrix.empty:
        raise ValueError(
            "No usable features remain after preprocessing. Check data type, scale, and missing values."
        )
    if matrix.isna().all(axis=1).any():
        raise ValueError("A sample has no usable measurements after preprocessing.")
    return matrix


def welch_test(matrix, labels, id_col="gene", effect_col="log_fc"):
    labels = labels.reindex(matrix.index)
    groups = labels.dropna().unique()
    if len(groups) != 2:
        return None
    a, b = (matrix.loc[labels == group] for group in groups)
    if min(len(a), len(b)) < 2:
        return None
    valid = (a.count() >= 2) & (b.count() >= 2)
    a, b = a.loc[:, valid], b.loc[:, valid]
    if a.shape[1] == 0:
        return None
    stat, p = ttest_ind(a, b, axis=0, equal_var=False, nan_policy="omit")
    result = pd.DataFrame(
        {
            id_col: a.columns,
            effect_col: a.mean() - b.mean(),
            "t_stat": stat,
            "p_value": p,
        }
    ).reset_index(drop=True)
    finite = np.isfinite(result["p_value"])
    result = result.loc[finite].copy()
    if result.empty:
        return None
    result["p_adj"] = multipletests(result["p_value"], method="fdr_bh")[1]
    return result.sort_values("p_value").reset_index(drop=True)
