"""
Shared ML steps applied after type-specific preprocessing/normalization.
Works on a samples x features matrix (rows = samples) plus optional labels.
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import silhouette_score
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from analysis_utils import prepare_processed


def _scale_features(matrix):
    # Filling a very wide DataFrame with a Series iterates over every feature
    # column. GEO matrices often have tens of thousands of columns here.
    values = matrix.to_numpy(dtype=float, copy=True)
    missing = np.isnan(values)
    if missing.any():
        values = np.where(missing, np.nanmean(values, axis=0), values)
    return StandardScaler().fit_transform(values)


def run_pca(matrix: pd.DataFrame, n_components: int = 2):
    """matrix: samples x features. Returns PCA coords + explained variance."""
    matrix = prepare_processed(matrix)
    scaled = _scale_features(matrix)
    if len(matrix) < 2 or not np.any(np.var(scaled, axis=0) > 0):
        return {"coords": np.zeros((len(matrix), 2)), "explained_variance_ratio": [0.0, 0.0],
                "sample_names": matrix.index.tolist()}
    n_components = min(n_components, matrix.shape[0] - 1, matrix.shape[1])
    pca = PCA(n_components=max(n_components, 1))
    coords = pca.fit_transform(scaled)
    return {
        "coords": coords,
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "sample_names": matrix.index.tolist(),
    }


def run_clustering(matrix: pd.DataFrame, k: int = None):
    """Unsupervised clustering with a simple heuristic for k if not given."""
    matrix = prepare_processed(matrix)
    scaled = _scale_features(matrix)
    n_samples = matrix.shape[0]

    if k is None:
        k = min(max(2, n_samples // 5), 6)
    k = min(k, n_samples - 1) if n_samples > 2 else 1
    k = min(k, len(np.unique(scaled, axis=0)))

    if k < 2:
        return {"labels": [0] * n_samples, "k": k, "silhouette": None}

    km = KMeans(n_clusters=k, n_init=10, random_state=42).fit(scaled)
    try:
        sil = silhouette_score(scaled, km.labels_)
    except ValueError:
        sil = None

    return {"labels": km.labels_.tolist(), "k": k, "silhouette": sil}


def train_classifier(matrix: pd.DataFrame, labels: pd.Series):
    """
    Supervised classification when phenotype labels are available (e.g. disease
    vs control). Uses Random Forest with cross-validation — appropriate for the
    classic p >> n bioinformatics setting.
    Returns model performance + top feature importances (candidate biomarkers).
    """
    matrix = prepare_processed(matrix)
    y_raw = labels.reindex(matrix.index)
    mask = y_raw.notna()
    X = matrix.loc[mask]
    y_raw = y_raw.loc[mask].astype(str)

    if y_raw.nunique() < 2 or len(y_raw) < 6:
        return {"trained": False, "reason": "insufficient labeled samples or classes"}

    le = LabelEncoder()
    y = le.fit_transform(y_raw)

    min_class_count = int(np.bincount(y).min())
    if min_class_count < 2:
        smallest = le.inverse_transform([np.bincount(y).argmin()])[0]
        return {
            "trained": False,
            "reason": f"class '{smallest}' has only {min_class_count} sample(s); "
                      f"need at least 2 per class for cross-validation",
        }

    # Fit missing-value imputation within each training fold.
    clf = make_pipeline(SimpleImputer(keep_empty_features=True),
                        RandomForestClassifier(n_estimators=500, random_state=42, n_jobs=-1))
    cv = StratifiedKFold(n_splits=min(5, min_class_count), shuffle=True, random_state=42)

    try:
        scores = cross_val_score(clf, X, y, cv=cv, scoring="accuracy")
    except ValueError as e:
        return {"trained": False, "reason": str(e)}

    clf.fit(X, y)
    importances = pd.Series(clf[-1].feature_importances_, index=X.columns).sort_values(ascending=False)

    return {
        "trained": True,
        "classes": le.classes_.tolist(),
        "cv_accuracy_mean": float(scores.mean()),
        "cv_accuracy_std": float(scores.std()),
        "top_features": importances.head(25).to_dict(),
    }
