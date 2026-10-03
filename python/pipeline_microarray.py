"""
Microarray pipeline (Affymetrix / Illumina / Agilent expression arrays).

Real-world normalization and differential expression use the field-standard
Bioconductor `limma` package (via r_bridge.py / rpy2) when R + limma are
installed. If they aren't, this module transparently falls back to a
Python-only implementation (quantile normalization + Welch's t-test) so the
pipeline still runs — just with less statistically rigorous results.

Note: this operates on an already-summarized probe/gene x sample intensity
matrix (e.g. a GEO series matrix), not raw CEL files. Raw CEL processing
(affy::rma()) would need actual CEL file uploads, which isn't supported yet.
"""

import sys

import numpy as np
import pandas as pd

try:
    import r_bridge

    R_AVAILABLE, _r_status_msg = r_bridge.check_available()
except ImportError:
    r_bridge = None
    R_AVAILABLE, _r_status_msg = False, "r_bridge module not found"

if not R_AVAILABLE:
    print(
        f"[microarray] R/limma not available ({_r_status_msg}); using Python fallback.",
        file=sys.stderr,
    )


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Simple quantile normalization fallback (probes x samples)."""
    ranked = df.rank(method="average")
    sorted_df = pd.DataFrame(
        np.sort(df.values, axis=0), index=df.index, columns=df.columns
    )
    mean_sorted = sorted_df.mean(axis=1)
    return ranked.apply(
        lambda col: np.interp(col, np.arange(1, len(mean_sorted) + 1), mean_sorted)
    )


def preprocess(expression_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    expression_matrix: probes x samples, raw intensities.
    Returns: samples x genes, normalized, ready for ML.
    """
    # GEO often supplies already log-transformed expression values.
    observed = expression_matrix.to_numpy().ravel()
    observed = observed[np.isfinite(observed)]
    needs_log = np.quantile(observed, 0.99) > 100 and observed.min() >= 0
    log_data = (
        np.log2(expression_matrix.clip(lower=1))
        if needs_log
        else expression_matrix.copy()
    )
    log_data = log_data.T.fillna(log_data.mean(axis=1)).T

    if R_AVAILABLE:
        try:
            normalized = r_bridge.normalize_between_arrays(log_data)
            print(
                "[microarray] Normalized with limma::normalizeBetweenArrays.",
                file=sys.stderr,
            )
            return normalized.T
        except r_bridge.RBridgeUnavailable as e:
            print(
                f"[microarray] limma normalization failed ({e}); falling back to Python quantile normalization.",
                file=sys.stderr,
            )

    normalized = quantile_normalize(log_data)
    return normalized.T  # samples x probes


def run_differential_expression(matrix: pd.DataFrame, group_labels: pd.Series):
    """
    Differential expression between exactly 2 groups.
    Uses limma (moderated t-statistics, empirical Bayes) when R is available —
    the field standard for microarray DE, and notably better-behaved than a
    plain t-test at small sample sizes. Falls back to Welch's t-test otherwise.
    """
    group_labels = group_labels.reindex(matrix.index)
    groups = group_labels.dropna().unique()
    if len(groups) != 2 or group_labels.value_counts().min() < 2:
        return None  # DE analysis needs exactly 2 groups here; extend for multi-group (ANOVA/limma contrasts)

    if R_AVAILABLE:
        try:
            de_df = r_bridge.limma_differential_expression(
                matrix, group_labels, id_col="gene"
            )
            print(
                "[microarray] Differential expression via limma (lmFit + eBayes).",
                file=sys.stderr,
            )
            return de_df.sort_values("p_value")
        except r_bridge.RBridgeUnavailable as e:
            print(
                f"[microarray] limma DE unavailable ({e}); falling back to Welch's t-test.",
                file=sys.stderr,
            )
        except Exception as e:
            print(
                f"[microarray] limma DE failed ({e}); falling back to Welch's t-test.",
                file=sys.stderr,
            )

    from analysis_utils import welch_test

    return welch_test(matrix, group_labels)
