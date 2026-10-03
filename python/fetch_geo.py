"""
Fetch a GSE series from GEO and detect its data type.

Detection strategy:
1. Pull series + platform metadata via GEOparse.
2. Inspect platform technology / title / GPL description for keywords.
3. Fall back to inspecting the supplementary file names/extensions
   (e.g. *_series_matrix.txt for arrays, counts.txt/tsv for RNA-seq,
   *_Grn.idat/*_Red.idat or beta-value matrices for methylation).
"""
import re
import gzip
import os
import tempfile
import socket
from contextlib import nullcontext
import zlib
from pathlib import Path

import GEOparse
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DATA_TYPE_MICROARRAY = "microarray"
DATA_TYPE_RNASEQ = "rnaseq"
DATA_TYPE_METHYLATION = "methylation"
DATA_TYPE_UNKNOWN = "unknown"

METHYLATION_KEYWORDS = re.compile(
    r"methylation|infinium|450k|epic|bisulfite|bead ?chip.*methyl", re.I
)
RNASEQ_KEYWORDS = re.compile(
    r"rna-?seq|high.throughput sequencing|illumina.*sequencing|rna sequencing", re.I
)
MICROARRAY_KEYWORDS = re.compile(
    r"expression.*array|microarray|affymetrix|agilent|illumina.*beadchip.*expr",
    re.I,
)


def _validate_soft(path, accession):
    """Check gzip integrity and series identity before publishing a cache file."""
    found = False
    pattern = re.compile(rb'^\^SERIES\s*=\s*' + accession.encode('ascii') + rb'\s*$', re.M)
    tail = b''
    with gzip.open(path, 'rb') as handle:
        while chunk := handle.read(1024 * 1024):
            data = tail + chunk
            found = found or bool(pattern.search(data))
            tail = data[-256:]
    if not found:
        raise ValueError(f'The downloaded file does not contain series {accession}.')


def _is_dns_error(error):
    pending, seen = [error], set()
    while pending:
        item = pending.pop()
        if id(item) in seen:
            continue
        seen.add(id(item))
        if isinstance(item, socket.gaierror):
            return True
        for linked in (getattr(item, '__cause__', None), getattr(item, '__context__', None),
                       getattr(item, 'reason', None), *getattr(item, 'args', ())):
            if isinstance(linked, BaseException):
                pending.append(linked)
    return False


def _download_soft(session, url, cache, accession, compressed=True):
    temporary = None
    try:
        with session.get(url, stream=True, timeout=(15, 120)) as response:
            if response.status_code == 404:
                raise ValueError(f'{accession} has no public GEO SOFT file. Check the accession and whether the study has been released.')
            response.raise_for_status()
            with tempfile.NamedTemporaryFile(dir=cache, prefix=accession + '-', suffix='.part', delete=False) as output:
                temporary = Path(output.name)
                # GEO's query endpoint returns plain SOFT; keep one gzip cache format.
                with nullcontext(output) if compressed else gzip.GzipFile(fileobj=output, mode='wb') as destination:
                    size = 0
                    for chunk in response.iter_content(1024 * 1024):
                        size += len(chunk)
                        if size > 500 * 1024 * 1024:
                            raise ValueError('GEO metadata download exceeds 500 MB. Upload the intended processed matrix instead.')
                        destination.write(chunk)
        _validate_soft(temporary, accession)
        os.replace(temporary, cache / f'{accession}_family.soft.gz')
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def fetch_gse(gse_id: str, destdir=None):
    """Download over HTTPS, atomically cache, and let GEOparse parse locally.

    GEOparse's FTP downloader can swallow a connection failure and then try to
    copy a nonexistent temporary file. Avoid that transport entirely.
    """
    accession = gse_id.strip().upper()
    if not re.fullmatch(r'GSE\d+', accession):
        raise ValueError('Provide a valid GEO Series accession, e.g. GSE12345.')
    cache = Path(destdir) if destdir is not None else Path(__file__).resolve().parent.parent / 'geo_cache'
    cache.mkdir(parents=True, exist_ok=True)
    target = cache / f'{accession}_family.soft.gz'
    valid_cache = False
    if target.is_file():
        try:
            _validate_soft(target, accession)
            valid_cache = True
        except (OSError, EOFError, ValueError, zlib.error):
            pass  # Replace incomplete old cache files only after a valid download.

    if not valid_cache:
        bucket = re.sub(r'\d{1,3}$', 'nnn', accession)
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{accession}/soft/{accession}_family.soft.gz'
        try:
            with requests.Session() as session:
                retries = Retry(total=2, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=['GET'])
                session.mount('https://', HTTPAdapter(max_retries=retries))
                try:
                    _download_soft(session, url, cache, accession)
                except requests.RequestException:
                    # A separate official NCBI hostname can remain reachable when
                    # the file server's DNS or connection is unavailable.
                    fallback = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=all&form=text&view=full'
                    _download_soft(session, fallback, cache, accession, compressed=False)
        except requests.RequestException as exc:
            if _is_dns_error(exc):
                message = (f'Cannot resolve NCBI download servers to download {accession}. '
                           "Check this computer's internet connection and DNS, then retry. "
                           'You can still analyze cached datasets or upload a processed matrix.')
            else:
                message = (f'Could not download {accession} from either NCBI server over HTTPS. '
                           'Check the network connection and retry, or upload a processed matrix.')
            raise RuntimeError(message) from exc
        except (OSError, EOFError, zlib.error) as exc:
            raise RuntimeError(f'Could not save or validate the download for {accession}. Check free disk space and cache permissions, then retry.') from exc

    return GEOparse.get_GEO(filepath=str(target), geotype='GSE', silent=True)


def detect_data_type(gse) -> str:
    """Inspect platform + sample metadata to classify the dataset."""
    text_blobs = []
    for key in ('type', 'title', 'summary', 'supplementary_file'):
        text_blobs.extend(gse.metadata.get(key, []))

    for gpl in gse.gpls.values():
        text_blobs.append(" ".join(str(v) for v in gpl.metadata.get("title", [])))
        text_blobs.append(" ".join(str(v) for v in gpl.metadata.get("technology", [])))

    for gsm in gse.gsms.values():
        text_blobs.append(" ".join(str(v) for v in gsm.metadata.get("library_strategy", [])))
        text_blobs.append(" ".join(str(v) for v in gsm.metadata.get("extract_protocol_ch1", [])))
        for supp in gsm.metadata.get("supplementary_file", []) or []:
            text_blobs.append(supp)

    combined = " ".join(text_blobs)

    if METHYLATION_KEYWORDS.search(combined) or ".idat" in combined.lower():
        return DATA_TYPE_METHYLATION
    # A sequencing instrument alone does not establish an RNA expression assay.
    if re.search(r'ChIP[- ]?seq|ATAC[- ]?seq|Hi-C', combined, re.I) and not re.search(r'RNA[- ]?seq', combined, re.I):
        return DATA_TYPE_UNKNOWN
    if RNASEQ_KEYWORDS.search(combined) or re.search(r'\bcounts?[_. -].*\.(txt|tsv|csv)', combined, re.I):
        return DATA_TYPE_RNASEQ
    if MICROARRAY_KEYWORDS.search(combined):
        return DATA_TYPE_MICROARRAY

    return DATA_TYPE_UNKNOWN


def get_phenotype_table(gse):
    """Build a samples x metadata dataframe (group labels, etc.) for downstream ML."""
    import pandas as pd

    rows = {}
    for gsm_name, gsm in gse.gsms.items():
        rows[gsm_name] = {
            "title": gsm.metadata.get("title", [""])[0],
            "source": gsm.metadata.get("source_name_ch1", [""])[0],
            "characteristics": "; ".join(gsm.metadata.get("characteristics_ch1", [])),
        }
    return pd.DataFrame.from_dict(rows, orient="index")
