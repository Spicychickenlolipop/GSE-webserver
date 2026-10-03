"""
Renders the final HTML report from analysis results.
Uses Plotly for interactive charts (PCA scatter, volcano plot) embedded directly
in the HTML — no server-side image rendering needed, works offline once loaded.
"""
import json
from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape

TEMPLATE_DIR = Path(__file__).parent / "templates"


def build_pca_chart_json(pca_result, group_labels=None):
    coords = pca_result["coords"]
    names = pca_result["sample_names"]
    import pandas as pd
    labels = [group_labels.get(n, "unknown") if group_labels else "sample" for n in names]
    labels = ['unknown' if pd.isna(label) else str(label) for label in labels]
    return {
        "x": coords[:, 0].tolist(),
        "y": (coords[:, 1].tolist() if coords.shape[1] > 1 else [0] * len(names)),
        "text": names,
        "labels": labels,
        "explained_variance": pca_result["explained_variance_ratio"],
    }


def build_volcano_chart_json(de_df, x_col="log_fc", id_col=None):
    if de_df is None:
        return None
    import numpy as np
    return {
        "x": de_df[x_col].tolist(),
        "y": (-np.log10(de_df["p_value"].clip(lower=1e-300))).tolist(),
        "text": de_df[id_col].tolist() if id_col else de_df.iloc[:, 0].tolist(),
    }


def generate_report(
    out_path: str,
    gse_id: str,
    data_type: str,
    n_samples: int,
    n_features: int,
    pca_result,
    clustering_result,
    de_result_df,
    de_id_col: str,
    de_x_col: str,
    classifier_result,
    group_labels=None,
    notes=None,
):
    env = Environment(loader=FileSystemLoader(str(TEMPLATE_DIR)), autoescape=select_autoescape(['html']))
    template = env.get_template("report_template.html")

    pca_json = build_pca_chart_json(pca_result, group_labels)
    volcano_json = build_volcano_chart_json(de_result_df, x_col=de_x_col, id_col=de_id_col)

    top_de_table = None
    if de_result_df is not None:
        top_de_table = de_result_df.head(25).to_dict(orient="records")

    html = template.render(
        gse_id=gse_id,
        data_type=data_type,
        n_samples=n_samples,
        n_features=n_features,
        pca_json=json.dumps(pca_json, allow_nan=False).replace('<', '\\u003c'),
        volcano_json=json.dumps(volcano_json, allow_nan=False).replace('<', '\\u003c') if volcano_json else "null",
        clustering=clustering_result,
        top_de_table=top_de_table,
        de_x_col=de_x_col,
        classifier=classifier_result,
        notes=notes or [],
    )

    Path(out_path).write_text(html, encoding="utf-8")
    return out_path
