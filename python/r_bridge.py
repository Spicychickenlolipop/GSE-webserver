"""
Optional R/Bioconductor bridge via rpy2. Currently used by pipeline_microarray.py
for limma (field-standard microarray normalization + differential expression).

Every public function here raises RBridgeUnavailable if R, rpy2, or the
required Bioconductor package isn't installed — callers should catch that
and fall back to the pure-Python implementation rather than crashing the
whole pipeline. This keeps the webserver usable for people who haven't set
up R, while giving correct, field-standard results for people who have.
"""
import pandas as pd


class RBridgeUnavailable(Exception):
    """Raised when R/rpy2/a required R package isn't available."""
    pass


def _get_limma():
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr
    except Exception as e:
        raise RBridgeUnavailable(f"rpy2 is not installed ({e}). Run: pip install rpy2")

    try:
        limma = importr("limma")
    except Exception as e:
        raise RBridgeUnavailable(
            "R and/or the Bioconductor 'limma' package aren't available "
            f"({e}). In R: BiocManager::install('limma')"
        )

    return ro, pandas2ri, localconverter, limma


def check_available() -> tuple:
    """Returns (available: bool, message: str) — used by check_r_setup.py."""
    try:
        _get_limma()
        return True, "R + rpy2 + limma are all available."
    except RBridgeUnavailable as e:
        return False, str(e)


def normalize_between_arrays(log_matrix: pd.DataFrame) -> pd.DataFrame:
    """
    log_matrix: probes x samples, already log2-transformed intensities.
    Runs limma::normalizeBetweenArrays (quantile method) — the field-standard
    between-array normalization for expression microarrays whose data arrives
    as an already-summarized intensity matrix (e.g. a GEO series matrix)
    rather than raw CEL files.
    Returns: probes x samples, normalized.
    """
    ro, pandas2ri, localconverter, limma = _get_limma()

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_matrix = ro.conversion.py2rpy(log_matrix)

    r_normalized = limma.normalizeBetweenArrays(r_matrix, method="quantile")

    with localconverter(ro.default_converter + pandas2ri.converter):
        normalized = ro.conversion.rpy2py(r_normalized)

    return pd.DataFrame(normalized, index=log_matrix.index, columns=log_matrix.columns)


def limma_differential_expression(
    matrix: pd.DataFrame, group_labels: pd.Series, id_col: str = "gene"
) -> pd.DataFrame:
    """
    matrix: samples x probes/genes, normalized log-intensities.
    group_labels: pandas Series indexed by sample, exactly 2 groups.

    Runs limma's moderated t-test (lmFit -> contrasts.fit -> eBayes -> topTable)
    — the field standard for microarray differential expression. The empirical
    Bayes shrinkage is specifically what makes this more reliable than a plain
    t-test at the small sample sizes typical of GEO datasets.

    Returns a dataframe with columns: [id_col, log_fc, t_stat, p_value, p_adj]
    (same schema the existing Python fallback produces, so callers and the
    report template don't need to change).
    """
    ro, pandas2ri, localconverter, limma = _get_limma()
    from rpy2.robjects.packages import importr
    stats_pkg = importr("stats")  # for model.matrix

    groups = list(group_labels.dropna().unique())
    if len(groups) != 2:
        raise ValueError("limma_differential_expression requires exactly 2 groups")

    aligned_labels = group_labels.reindex(matrix.index)
    matrix = matrix.loc[aligned_labels.notna()]
    aligned_labels = aligned_labels.dropna().astype(str)
    expr = matrix.T  # limma expects probes/genes x samples

    with localconverter(ro.default_converter + pandas2ri.converter):
        r_expr = ro.conversion.py2rpy(expr)
        r_group = ro.conversion.py2rpy(aligned_labels.values)

    ro.globalenv["expr"] = r_expr
    ro.globalenv["group"] = ro.r["factor"](r_group, levels=ro.StrVector(groups))

    ro.r("design <- model.matrix(~0 + group)")
    ro.r(f'colnames(design) <- c("g0", "g1")')
    ro.r("contrast <- limma::makeContrasts(g0 - g1, levels = design)")
    ro.r("fit <- limma::lmFit(expr, design)")
    ro.r("fit2 <- limma::contrasts.fit(fit, contrast)")
    ro.r("fit2 <- limma::eBayes(fit2)")
    ro.r("tt <- limma::topTable(fit2, number = Inf, sort.by = 'P')")

    with localconverter(ro.default_converter + pandas2ri.converter):
        tt = ro.conversion.rpy2py(ro.r("tt"))

    tt = tt.reset_index().rename(
        columns={
            "index": id_col,
            "logFC": "log_fc",
            "t": "t_stat",
            "P.Value": "p_value",
            "adj.P.Val": "p_adj",
        }
    )
    return tt[[id_col, "log_fc", "t_stat", "p_value", "p_adj"]]
