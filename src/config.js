const path = require("path");

const PROJECT_ROOT = path.resolve(__dirname, "..");

module.exports = {
  UPLOADS_DIR: path.join(PROJECT_ROOT, "uploads"),
  REPORTS_DIR: path.join(PROJECT_ROOT, "reports"),
  PIPELINE_ENTRY: path.join(PROJECT_ROOT, "python", "run_pipeline.py"),
  PYTHON_BIN: process.env.PYTHON_BIN || (process.platform === "win32" ? "python" : "python3"),
  JOB_OPTIONS: { removeOnComplete: false, removeOnFail: false, attempts: 1 },
};
