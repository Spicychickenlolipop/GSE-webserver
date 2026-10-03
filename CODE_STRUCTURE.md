# Code structure

The browser submits an analysis request to Express. Express stores the job in
Redis, and a BullMQ worker runs the Python pipeline. Python writes an HTML report
and sends progress messages back to the worker.

## Frontend

- `frontend/index.html`: page markup and input controls.
- `frontend/styles.css`: layout, colours, and responsive styles.
- `frontend/app.js`: event listeners, uploads, requests, and status polling.

## API and worker

- `src/index.js`: starts Express and, by default, the worker.
- `src/config.js`: shared paths, Python executable, and job options.
- `src/routes/analyze.js`: request validation, job submission, status, and reports.
- `src/middleware/uploads.js`: accepted file formats and upload storage.
- `src/queue/queue.js`: Redis connection and BullMQ queue.
- `src/queue/worker.js`: worker registration and lifecycle logging.
- `src/queue/process_job.js`: Python arguments, execution, and progress parsing.

## Analysis

- `python/run_pipeline.py`: command-line arguments and orchestration of one job.
- `python/fetch_geo.py`: GEO downloads, cache validation, and assay detection.
- `python/matrix_io.py`: matrix parsing, validation, and supplementary-file loading.
- `python/labels.py`: uploaded labels and metadata-based group inference.
- `python/pipeline_microarray.py`, `pipeline_rnaseq.py`, `pipeline_methylation.py`:
  assay-specific preprocessing and differential analysis.
- `python/analysis_utils.py`: shared validation and statistical helpers.
- `python/ml_models.py`: PCA, clustering, and Random Forest classification.
- `python/r_bridge.py`: optional R/limma integration.
- `python/report_generator.py` and `python/templates/`: HTML report generation.

Matrices enter as features by samples and are transposed by preprocessing.
Shared analysis functions expect samples by features. Keep this distinction
explicit when adding a new method.

## Running and checking changes

Start Redis, then run `npm start` to launch the API and worker. Open
`http://localhost:4000`. Set `RUN_WORKER=false` on the API process if using a
separate `npm run worker` process.

Run JavaScript tests with `npm test` and Python tests with `npm run test:python`.
Tests cover requests, worker messages, browser event wiring, downloads, and
complete report generation using synthetic matrices. Download tests mock NCBI
responses and do not require internet access.

Python files use Black formatting. If Black is installed, run
`python -m black python` after Python edits. `.editorconfig` supplies indentation
settings for both languages.

## Current analysis limits

This cleanup does not add scientific methods. RNA-seq still uses the CPM/Welch
fallback; DESeq2, raw IDAT/minfi processing, and single-cell analysis are not
implemented. Microarrays use limma when available and a Python fallback otherwise.
The classifier uses retained features independently of differential testing.
Significant-feature selection and multiple-model comparison remain future work.
