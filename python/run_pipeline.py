"""Run one analysis job and send progress messages to the Node worker.

Input matrices have features in rows and samples in columns. Each assay's
preprocessor returns samples in rows for the shared analysis steps.
"""

import argparse
import importlib
import json
import sys
import traceback

import fetch_geo
import ml_models
import report_generator
from analysis_utils import prepare_processed
from labels import infer_group_labels, load_labels_from_file
from matrix_io import load_geo_matrix, load_matrix, validate_matrix


PIPELINE_MODULES = {
    "microarray": "pipeline_microarray",
    "rnaseq": "pipeline_rnaseq",
    "methylation": "pipeline_methylation",
}


def progress(step, data_type=None):
    message = f"PROGRESS:{step}"
    if data_type:
        message += f":{data_type}"
    print(message, flush=True)


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", help="GEO Series accession, e.g. GSE12345")
    parser.add_argument(
        "--matrix-file", help="Path to an uploaded features x samples matrix (csv/tsv)"
    )
    parser.add_argument(
        "--labels-file",
        help="Path to an uploaded two-column sample_id,group file (csv/tsv)",
    )
    parser.add_argument(
        "--data-type",
        choices=["microarray", "rnaseq", "methylation"],
        help="Required with --matrix-file; identifies which pipeline to run",
    )
    parser.add_argument("--job-id", required=True)
    parser.add_argument(
        "--data-scale",
        choices=["auto", "processed"],
        default="auto",
        help="processed skips normalization for already normalized expression or methylation M-values",
    )
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    if not args.gse and not args.matrix_file:
        print("Either --gse or --matrix-file must be provided.", file=sys.stderr)
        sys.exit(1)
    if args.matrix_file and not args.data_type:
        print("--data-type is required when using --matrix-file.", file=sys.stderr)
        sys.exit(1)

    return args


def load_input(args):
    """Return the dataset name, assay type, feature matrix, and metadata."""
    gse_label = (
        args.gse if args.gse else args.matrix_file.split("/")[-1].split("\\")[-1]
    )

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

    return gse_label, data_type, matrix, pheno


def get_pipeline(data_type):
    # Load only the selected assay; the optional R setup can be expensive.
    return importlib.import_module(PIPELINE_MODULES[data_type])


def build_analysis_notes(data_scale, group_labels, differential_result):
    if data_scale == "processed":
        preprocessing = "Normalization skipped: input was marked already processed."
    else:
        preprocessing = "Type-specific preprocessing was applied."

    if group_labels is None:
        labels = "No unambiguous group labels were available; supervised analysis was skipped."
    else:
        labels = "Supervised analyses use only samples with group labels."

    if differential_result is None:
        differential = (
            "Differential analysis requires exactly two groups with at least two "
            "measured samples per group and testable features."
        )
    else:
        differential = (
            "Differential results are exploratory; review the selected groups "
            "and preprocessing before interpretation."
        )
    return [preprocessing, labels, differential]


def run_analysis(args):
    dataset_name, data_type, matrix, pheno = load_input(args)
    pipeline = get_pipeline(data_type)
    matrix = validate_matrix(matrix)
    processed = prepare_processed(
        matrix.T if args.data_scale == "processed" else pipeline.preprocess(matrix)
    )

    # Labels are optional; PCA and clustering still run without them.
    if args.labels_file:
        group_labels = load_labels_from_file(args.labels_file, processed.index)
    else:
        group_labels = infer_group_labels(pheno, processed.index)

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
            de_result = pipeline.run_differential_methylation(processed, group_labels)
            de_id_col, de_x_col = "probe", "delta_m"
        else:
            de_result = pipeline.run_differential_expression(processed, group_labels)
            de_id_col, de_x_col = "gene", "log_fc"

    progress("rendering_report", data_type)
    group_label_map = group_labels.to_dict() if group_labels is not None else None
    report_generator.generate_report(
        out_path=args.out,
        gse_id=dataset_name,
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
        notes=build_analysis_notes(args.data_scale, group_labels, de_result),
    )

    progress("done", data_type)


def main():
    args = parse_arguments()
    try:
        run_analysis(args)
    except Exception as exc:
        print("ERROR:" + json.dumps({"message": str(exc)}), flush=True)
        traceback.print_exc()
        sys.exit(1)


def load_matrix_for_gse(gse, data_type):
    return load_geo_matrix(gse, data_type)


def load_matrix_from_file(path):
    return load_matrix(path)


if __name__ == "__main__":
    main()
