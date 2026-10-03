"""
Entry point spawned by the Node.js worker (src/queue/worker.js) for each job.
Prints PROGRESS:<step>[:<data_type>] lines to stdout so Node can report status.

Two input modes:
  1. GEO fetch:      python run_pipeline.py --gse GSE12345 --job-id <uuid> --out report.html
  2. Uploaded file:  python run_pipeline.py --matrix-file data.csv --data-type rnaseq
                        [--labels-file labels.csv] --job-id <uuid> --out report.html

Uploaded matrix file convention: CSV/TSV, first column = feature ID (gene/probe),
first row = sample names as header — i.e. features x samples, same orientation
GEO series matrices use. Values: raw counts (rnaseq), intensities (microarray),
or beta values 0-1 (methylation).

Uploaded labels file convention: CSV/TSV with two columns, "sample_id" and
"group" (header required, any column names accepted as long as there are
exactly two columns) — used for classification + differential analysis.
"""
import argparse
import json
import sys
import traceback

import pandas as pd

import fetch_geo
import ml_models
import report_generator
from matrix_io import load_matrix, load_geo_matrix, read_table, validate_matrix
from analysis_utils import prepare_processed


def progress(step, data_type=None):
    if data_type:
        print(f"PROGRESS:{step}:{data_type}", flush=True)
    else:
        print(f"PROGRESS:{step}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", help="GEO Series accession, e.g. GSE12345")
    parser.add_argument("--matrix-file", help="Path to an uploaded features x samples matrix (csv/tsv)")
    parser.add_argument("--labels-file", help="Path to an uploaded two-column sample_id,group file (csv/tsv)")
    parser.add_argument("--data-type", choices=["microarray", "rnaseq", "methylation"],
                         help="Required with --matrix-file; identifies which pipeline to run")
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--data-scale", choices=["auto", "processed"], default="auto",
                        help="processed skips normalization for already normalized expression or methylation M-values")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if not args.gse and not args.matrix_file:
        print("Either --gse or --matrix-file must be provided.", file=sys.stderr)
        sys.exit(1)
    if args.matrix_file and not args.data_type:
        print("--data-type is required when using --matrix-file.", file=sys.stderr)
        sys.exit(1)

    try:
        gse_label = args.gse if args.gse else args.matrix_file.split("/")[-1].split("\\")[-1]

        if args.gse:
            progress("fetching")
            gse = fetch_geo.fetch_gse(args.gse)

            data_type = args.data_type or fetch_geo.detect_data_type(gse)
            progress("detected", data_type)

            if data_type == fetch_geo.DATA_TYPE_UNKNOWN:
                print(
                    f"Could not confidently detect data type for {args.gse}. "
                    "Select its data type explicitly or upload a processed matrix.",
                    file=sys.stderr,
                )
                sys.exit(1)

            pheno = fetch_geo.get_phenotype_table(gse)

            # Read sample tables or unambiguous series-level processed matrices.
            progress("preprocessing", data_type)
            matrix = load_matrix_for_gse(gse, data_type)
        else:
            data_type = args.data_type
            progress("loading_upload", data_type)
            matrix = load_matrix_from_file(args.matrix_file)
            pheno = None

        if data_type == fetch_geo.DATA_TYPE_MICROARRAY:
            import pipeline_microarray as pipe
        elif data_type == fetch_geo.DATA_TYPE_RNASEQ:
            import pipeline_rnaseq as pipe
        else:
            import pipeline_methylation as pipe

        matrix = validate_matrix(matrix)
        processed = prepare_processed(matrix.T if args.data_scale == "processed" else pipe.preprocess(matrix))

        # --- Try to find a usable group label for supervised steps ---
        if args.labels_file:
            group_labels = load_labels_from_file(args.labels_file, processed.index)
        else:
            group_labels = infer_group_labels(pheno, processed.index)

        # --- ML: PCA, clustering, classification ---
        progress("ml_pca", data_type)
        pca_result = ml_models.run_pca(processed)

        progress("ml_clustering", data_type)
        clustering_result = ml_models.run_clustering(processed)

        classifier_result = None
        de_result = None
        de_id_col, de_x_col = None, None
        if group_labels is not None and group_labels.nunique() >= 2:
            progress("ml_classification", data_type)
            classifier_result = ml_models.train_classifier(processed, group_labels)

            progress("differential_analysis", data_type)
            if data_type == fetch_geo.DATA_TYPE_METHYLATION:
                de_result = pipe.run_differential_methylation(processed, group_labels)
                de_id_col, de_x_col = "probe", "delta_m"
            else:
                de_result = pipe.run_differential_expression(processed, group_labels)
                de_id_col, de_x_col = "gene", "log_fc"

        # --- Report ---
        progress("rendering_report", data_type)
        group_label_map = group_labels.to_dict() if group_labels is not None else None
        report_generator.generate_report(
            out_path=args.out,
            gse_id=gse_label,
            data_type=data_type,
            n_samples=processed.shape[0],
            n_features=processed.shape[1],
            pca_result=pca_result,
            clustering_result=clustering_result,
            de_result_df=de_result,
            de_id_col=de_id_col,
            de_x_col=de_x_col,
            classifier_result=classifier_result,
            group_labels=group_label_map,
            notes=[
                "Normalization skipped: input was marked already processed." if args.data_scale == "processed" else "Type-specific preprocessing was applied.",
                "No unambiguous group labels were available; supervised analysis was skipped." if group_labels is None else "Supervised analyses use only samples with group labels.",
                "Differential analysis requires exactly two groups with at least two measured samples per group and testable features." if de_result is None else "Differential results are exploratory; review the selected groups and preprocessing before interpretation.",
            ],
        )

        progress("done", data_type)

    except Exception as exc:
        print('ERROR:' + json.dumps({'message': str(exc)}), flush=True)
        traceback.print_exc()
        sys.exit(1)


def load_matrix_for_gse(gse, data_type):
    return load_geo_matrix(gse, data_type)


def load_matrix_from_file(path):
    return load_matrix(path)


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
    labels = pd.Series(df.iloc[:, 1].astype('string').str.strip().replace('', pd.NA).to_numpy(), index=ids)
    aligned = labels.reindex(sample_ids)
    if aligned.notna().sum() == 0:
        raise ValueError("No non-empty labels match matrix sample IDs.")
    return aligned


def infer_group_labels(pheno, sample_ids):
    """
    Heuristic: use the 'characteristics' column to derive a group label if it looks
    like a simple 2-3 category variable (e.g. "disease: cancer" vs "disease: normal").
    TODO: replace with a UI step where the user picks which metadata field is the
    label of interest — auto-inference is a convenience default, not reliable enough
    to trust blindly for every dataset.
    """
    import pandas as pd

    if pheno is None or pheno.empty:
        return None

    if "characteristics" not in pheno:
        return None
    # Only use a single metadata field when its assignment is unambiguous.
    fields = {}
    for sample, characteristics in pheno["characteristics"].items():
        for item in str(characteristics).split(';'):
            key, sep, value = item.partition(':')
            if sep and value.strip():
                fields.setdefault(key.strip().lower(), {})[sample] = value.strip()
    candidates = []
    for values in fields.values():
        candidate = pd.Series(values).reindex(sample_ids)
        if candidate.notna().all() and 2 <= candidate.nunique() <= 6 and candidate.value_counts().min() >= 2:
            candidates.append(candidate)
    return candidates[0] if len(candidates) == 1 else None


if __name__ == "__main__":
    main()
