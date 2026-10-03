"""
RNA-seq pipeline. Assumes GEO gives a processed counts matrix (genes x samples) —
most GSEs do, as supplementary files. If only raw FASTQ is available, that needs
an alignment/quantification step (STAR/Salmon) upstream of this, which is out of
scope for a report webserver and should be flagged to the user instead of attempted
inline.

Field-standard normalization + DE is DESeq2 or edgeR (R/Bioconductor). This module
gives a Python-only fallback (CPM + log transform, filtering) plus a hook for the
real DESeq2 call via rpy2.
"""
import numpy as np
import pandas as pd


def filter_low_counts(counts: pd.DataFrame, min_count: int = 10, min_samples_frac: float = 0.2) -> pd.DataFrame:
    """counts: genes x samples raw counts."""
    min_samples = max(1, int(counts.shape[1] * min_samples_frac))
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def cpm_log_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    """Fallback normalization: counts-per-million, log2(CPM + 1)."""
    lib_sizes = counts.sum(axis=0)
    cpm = counts.div(lib_sizes.replace(0, 1), axis=1) * 1e6
    return np.log2(cpm + 1)


def preprocess(counts_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    counts_matrix: genes x samples, raw integer counts.
    Returns: samples x genes, normalized, ready for ML.
    """
    if counts_matrix.isna().any().any():
        raise ValueError('Raw RNA-seq counts contain missing values. Supply complete counts; missing counts are not assumed to be zero.')
    if (counts_matrix < 0).any().any():
        raise ValueError('Raw RNA-seq counts must be non-negative. For normalized or log expression, select already processed input.')
    filtered = filter_low_counts(counts_matrix)
    if filtered.empty:
        # Small/low-depth matrices still support exploratory analysis.
        filtered = counts_matrix.loc[counts_matrix.sum(axis=1) > 0]
    if filtered.empty:
        raise ValueError('The RNA-seq matrix contains no positive counts.')
    normalized = cpm_log_normalize(filtered)

    # TODO (recommended): replace CPM/log with DESeq2's median-of-ratios
    # normalization + variance-stabilizing transform via rpy2:
    #   from rpy2.robjects.packages import importr
    #   deseq2 = importr('DESeq2')
    #   dds = deseq2.DESeqDataSetFromMatrix(countData=..., colData=..., design=~condition)
    #   dds = deseq2.DESeq(dds)
    #   vsd = deseq2.vst(dds)

    return normalized.T  # samples x genes


def run_differential_expression(matrix: pd.DataFrame, group_labels: pd.Series):
    """
    Placeholder DE (Welch's t-test on log-CPM per gene).
    TODO: replace with DESeq2's Wald test / negative-binomial GLM via rpy2 —
    this is the field standard and properly models count-data overdispersion,
    which a t-test on log-CPM does not.
    """
    from analysis_utils import welch_test
    return welch_test(matrix, group_labels, 'gene', 'log_fc')
