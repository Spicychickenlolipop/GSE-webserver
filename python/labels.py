"""Read sample labels and infer unambiguous groups from GEO metadata."""

import pandas as pd

from matrix_io import read_table


def load_labels_from_file(path, sample_ids):
    """
    Load a user-uploaded two-column labels file (sample_id, group).
    Returns a pandas Series indexed by sample ID, aligned to sample_ids.
    """
    df = read_table(path)
    if df.shape[1] != 2:
        raise ValueError("Labels file must have exactly two columns: sample_id, group.")
    if df.iloc[:, 0].isna().any():
        raise ValueError("Labels contain empty sample IDs.")
    ids = df.iloc[:, 0].astype(str).str.strip()
    if ids.duplicated().any():
        raise ValueError("Labels contain duplicate sample IDs.")
    labels = pd.Series(
        df.iloc[:, 1].astype("string").str.strip().replace("", pd.NA).to_numpy(),
        index=ids,
    )
    aligned = labels.reindex(sample_ids)
    if aligned.notna().sum() == 0:
        raise ValueError("No non-empty labels match matrix sample IDs.")
    return aligned


def infer_group_labels(pheno, sample_ids):
    """Use metadata only when exactly one complete grouping is unambiguous."""
    if pheno is None or pheno.empty:
        return None

    if "characteristics" not in pheno:
        return None
    # Only use a single metadata field when its assignment is unambiguous.
    fields = {}
    for sample, characteristics in pheno["characteristics"].items():
        for item in str(characteristics).split(";"):
            key, sep, value = item.partition(":")
            if sep and value.strip():
                fields.setdefault(key.strip().lower(), {})[sample] = value.strip()
    candidates = []
    for values in fields.values():
        candidate = pd.Series(values).reindex(sample_ids)
        if (
            candidate.notna().all()
            and 2 <= candidate.nunique() <= 6
            and candidate.value_counts().min() >= 2
        ):
            candidates.append(candidate)
    return candidates[0] if len(candidates) == 1 else None
