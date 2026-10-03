"""
DNA methylation pipeline (Illumina 450K / EPIC arrays).

GEO typically provides either raw IDAT files (best, but needs minfi in R to
process properly) or already-processed beta-value matrices (probes x samples,
values in [0,1]) as supplementary files — much more common for a
report-generation webserver to consume directly.

Field-standard QC/normalization is minfi (R/Bioconductor). This module
handles the common case (beta-value matrix already provided) and gives a
hook for full IDAT processing via rpy2 + minfi if raw IDATs are supplied.
"""
import numpy as np
import pandas as pd


def beta_to_mvalue(beta: pd.DataFrame) -> pd.DataFrame:
    """M-values are more statistically valid for downstream stats/ML than beta values
    (better homoscedasticity), per Du et al. 2010."""
    eps = 1e-6
    beta_clipped = beta.clip(eps, 1 - eps)
    return np.log2(beta_clipped / (1 - beta_clipped))


def filter_probes(beta: pd.DataFrame, detection_p: pd.DataFrame = None, p_threshold: float = 0.01) -> pd.DataFrame:
    """Drop probes with poor detection p-values if provided; drop probes with >20% missing."""
    if detection_p is not None:
        bad_probes = (detection_p > p_threshold).mean(axis=1) > 0.05
        beta = beta.loc[~bad_probes]
    missing_frac = beta.isna().mean(axis=1)
    return beta.loc[missing_frac <= 0.2]


def preprocess(beta_matrix: pd.DataFrame, detection_p: pd.DataFrame = None) -> pd.DataFrame:
    """
    beta_matrix: probes x samples, beta values [0,1].
    Returns: samples x probes M-values, ready for ML.
    """
    if ((beta_matrix < 0) | (beta_matrix > 1)).any().any():
        raise ValueError('Methylation input must contain beta values between 0 and 1. For M-values, select already processed input.')
    filtered = filter_probes(beta_matrix, detection_p)
    if filtered.empty:
        raise ValueError('No methylation probes remain after filtering missing values (maximum 20% per probe).')
    filtered = filtered.T.fillna(filtered.mean(axis=1)).T
    mvalues = beta_to_mvalue(filtered)

    # TODO (recommended, if raw IDATs are available instead of a beta matrix):
    # process via minfi in R for proper background correction + normalization
    # (e.g. functional normalization, Noob):
    #   from rpy2.robjects.packages import importr
    #   minfi = importr('minfi')
    #   rgset = minfi.read_metharray_exp(...)
    #   mset = minfi.preprocessFunnorm(rgset)
    #   beta = minfi.getBeta(mset)
    # Also consider cell-type deconvolution (minfi::estimateCellCounts) if the
    # tissue is blood, since methylation is highly cell-composition-dependent.

    return mvalues.T  # samples x probes


def run_differential_methylation(matrix: pd.DataFrame, group_labels: pd.Series):
    """
    Placeholder differential methylation (Welch's t-test per probe on M-values).
    TODO: replace with limma on M-values (field-standard for DMP calling) or
    bumphunter/DMRcate for differentially methylated *regions* rather than
    single probes, via rpy2.
    """
    from analysis_utils import welch_test
    return welch_test(matrix, group_labels, 'probe', 'delta_m')
