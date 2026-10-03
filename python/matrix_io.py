"""Read processed feature-by-sample tables without relying on file extensions."""

import csv
import gzip
import io
import re
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pandas as pd


def read_table(path):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        lines, size = [], 0
        for line in handle:
            size += len(line)
            if size > 500 * 1024 * 1024:
                raise ValueError("Uncompressed table exceeds the 500 MB size limit.")
            if line.strip() and not line.startswith(("#", "!")):
                lines.append(line)
    if not lines:
        raise ValueError("The file contains no data table.")
    try:
        sep = csv.Sniffer().sniff("".join(lines[:20]), delimiters=",\t;").delimiter
    except csv.Error:
        sep = r"\s+"
    # Check before pandas can silently rename duplicate headers.
    header = (
        re.split(sep, lines[0].strip())
        if sep == r"\s+"
        else next(csv.reader([lines[0]], delimiter=sep))
    )
    names = [name.strip() for name in header[1:]]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate column names; sample IDs must be unique.")
    return pd.read_csv(
        io.StringIO("".join(lines)), sep=sep, engine="python", dtype={header[0]: str}
    )


def validate_matrix(frame):
    if frame.empty:
        raise ValueError("No numeric features or samples were found in the matrix.")
    frame = frame.copy()
    for axis in (frame.index, frame.columns):
        if axis.isna().any() or any(not str(v).strip() for v in axis):
            raise ValueError("Feature and sample IDs must not be empty.")
    frame.index = frame.index.astype(str).str.strip()
    frame.columns = frame.columns.astype(str).str.strip()
    if not frame.index.is_unique or not frame.columns.is_unique:
        raise ValueError("Duplicate feature or sample IDs; provide unique identifiers.")
    numeric = frame.apply(pd.to_numeric, errors="coerce")
    invalid = frame.notna() & numeric.isna()
    if invalid.any().any():
        raise ValueError(
            "The sample matrix contains non-numeric values. Remove annotation columns or use NA for missing values."
        )
    if np.isinf(numeric.to_numpy()).any():
        raise ValueError(
            "The matrix contains infinite values; provide finite measurements or NA."
        )
    if numeric.isna().all(axis=0).any():
        raise ValueError("At least one sample has no measurements.")
    numeric = numeric.dropna(how="all")
    if numeric.empty:
        raise ValueError("No measured features remain in the matrix.")
    return numeric


def load_matrix(path, sample_ids=None, sample_aliases=None):
    table = read_table(path)
    if table.shape[1] < 2:
        raise ValueError("Expected a feature ID column followed by sample columns.")
    frame = table.set_index(table.columns[0])
    frame.columns = frame.columns.astype(str).str.strip()
    frame.index = frame.index.map(lambda v: str(v).strip() if pd.notna(v) else v)
    if sample_aliases:
        frame = frame.rename(index=sample_aliases, columns=sample_aliases)
    if sample_ids is not None:
        samples = set(map(str, sample_ids))
        columns = [c for c in frame.columns if c in samples]
        rows = [r for r in frame.index if r in samples]
        if columns and rows:
            raise ValueError("Ambiguous sample orientation in the matrix.")
        if rows:
            frame = frame.loc[rows].T
        elif columns:
            frame = frame.loc[:, columns]
        else:
            raise ValueError(
                "Matrix sample names do not match GEO GSM IDs or unique sample titles; upload it with an explicit labels file."
            )
    return validate_matrix(frame)


def load_geo_matrix(gse, data_type):
    platforms = {
        tuple(gsm.metadata.get("platform_id", [])) for gsm in gse.gsms.values()
    }
    if len(platforms) > 1:
        raise ValueError(
            "This series contains multiple platforms. Upload a processed matrix for one platform at a time."
        )
    # Different GEO submitters use different measurement column names.
    aliases = {
        "microarray": ["VALUE", "signal", "intensity"],
        "rnaseq": ["VALUE", "count", "counts", "raw_count", "raw_counts"],
        "methylation": ["VALUE", "beta", "beta_value", "beta-value"],
    }[data_type]
    series = []
    for name, gsm in gse.gsms.items():
        table = gsm.table
        columns = {str(c).lower(): c for c in table.columns}
        value = next(
            (columns[a.lower()] for a in aliases if a.lower() in columns), None
        )
        feature = next(
            (
                columns[a]
                for a in ("id_ref", "id", "gene_id", "gene", "probe_id")
                if a in columns
            ),
            None,
        )
        if value is None or feature is None or table.empty:
            break
        item = table.set_index(feature)[value].rename(name)
        if not item.index.is_unique:
            raise ValueError(f"{name} contains duplicate feature IDs.")
        series.append(item)
    if series and len(series) == len(gse.gsms):
        matrix = pd.concat(series, axis=1, join="inner")
        if not matrix.empty:
            return validate_matrix(matrix)

    # Restrict automatic network fetching to GEO's public file hosts.
    urls = []
    for key, values in gse.metadata.items():
        if key.startswith("supplementary_file"):
            urls.extend(values if isinstance(values, list) else [values])
    candidates = []
    for url in dict.fromkeys(urls):
        parsed = urlparse(url)
        if (
            parsed.hostname in {"ftp.ncbi.nlm.nih.gov", "www.ncbi.nlm.nih.gov"}
            and parsed.scheme in {"ftp", "https", "http"}
            and re.search(r"\.(csv|tsv|txt)(\.gz)?$", parsed.path, re.I)
        ):
            candidates.append(url)
    if not candidates:
        raise ValueError(
            "No supported processed matrix found. Upload a CSV/TSV/TXT matrix (optionally gzip compressed). Raw FASTQ, CEL, IDAT and archive files require upstream processing."
        )
    import requests
    import tempfile

    titles = {}
    for sample, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [""])[0].strip()
        if title and title not in gse.gsms:
            titles.setdefault(title, []).append(sample)
    aliases = {
        title: samples[0] for title, samples in titles.items() if len(samples) == 1
    }
    matrices, errors = [], []
    with tempfile.TemporaryDirectory() as directory:
        for i, url in enumerate(candidates):
            safe_url = re.sub(r"^(ftp|http):", "https:", url)
            path = Path(directory) / (
                str(i) + (".txt.gz" if url.lower().endswith(".gz") else ".txt")
            )
            try:
                with requests.get(
                    safe_url, stream=True, timeout=(15, 120), allow_redirects=False
                ) as response:
                    response.raise_for_status()
                    if response.status_code != 200:
                        raise ValueError("Download did not return a matrix.")
                    size = 0
                    with path.open("wb") as output:
                        for chunk in response.iter_content(1024 * 1024):
                            size += len(chunk)
                            if size > 500 * 1024 * 1024:
                                raise ValueError(
                                    "Supplementary matrix exceeds the 500 MB download limit."
                                )
                            output.write(chunk)
                matrix = load_matrix(path, gse.gsms.keys(), aliases)
                if set(matrix.columns) != set(gse.gsms):
                    raise ValueError(
                        "Matrix does not cover every sample; upload the intended sample subset explicitly."
                    )
                matrices.append(matrix)
            except (ValueError, OSError, requests.RequestException) as exc:
                errors.append(f"{Path(urlparse(url).path).name}: {exc}")
    if len(matrices) == 1:
        return matrices[0]
    if len(matrices) > 1:
        raise ValueError(
            "Multiple processed matrices match this series. Upload the matrix you intend to analyze."
        )
    raise ValueError(
        "Could not load supplementary matrices. "
        + "; ".join(errors)
        + " Upload a processed matrix with matching sample IDs."
    )
