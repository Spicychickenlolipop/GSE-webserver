const express = require("express");
const fs = require("fs");
const path = require("path");
const upload = require("../middleware/uploads");
const { REPORTS_DIR, JOB_OPTIONS } = require("../config");
const { analysisQueue } = require("../queue/queue");

const router = express.Router();

const GSE_REGEX = /^GSE\d+$/i;
const VALID_DATA_TYPES = ["microarray", "rnaseq", "methylation"];

// POST /api/analyze  { gseId: "GSE12345" }
router.post("/analyze", async (req, res) => {
  const { gseId, dataType, dataScale = "auto" } = req.body || {};

  if (typeof gseId !== "string" || !GSE_REGEX.test(gseId.trim())) {
    return res.status(400).json({ error: "Provide a valid GSE accession, e.g. GSE12345" });
  }
  if ((dataType && !VALID_DATA_TYPES.includes(dataType)) || !["auto", "processed"].includes(dataScale)) {
    return res.status(400).json({ error: "Invalid data type or data scale." });
  }

  const job = await analysisQueue.add(
    "analyze-gse",
    { inputType: "gse", gseId: gseId.trim().toUpperCase(), dataType, dataScale },
    JOB_OPTIONS
  );

  res.json({ jobId: job.id, gseId: gseId.trim().toUpperCase() });
});

// POST /api/upload  (multipart/form-data)
// fields: matrixFile (required), labelsFile (optional), dataType (required)
router.post(
  "/upload",
  upload.fields([
    { name: "matrixFile", maxCount: 1 },
    { name: "labelsFile", maxCount: 1 },
  ]),
  async (req, res) => {
    try {
      const { dataType, dataScale = "auto" } = req.body || {};
      const matrixFile = req.files?.matrixFile?.[0];
      const labelsFile = req.files?.labelsFile?.[0];

      if (!matrixFile) {
        return res.status(400).json({ error: "A matrix file is required." });
      }
      if (!dataType || !VALID_DATA_TYPES.includes(dataType)) {
        return res.status(400).json({
          error: `dataType must be one of: ${VALID_DATA_TYPES.join(", ")}`,
        });
      }
      if (!["auto", "processed"].includes(dataScale)) {
        return res.status(400).json({ error: "Invalid data scale." });
      }

      const job = await analysisQueue.add(
        "analyze-upload",
        {
          inputType: "upload",
          dataType,
          dataScale,
          matrixPath: matrixFile.path,
          labelsPath: labelsFile ? labelsFile.path : null,
          originalMatrixName: matrixFile.originalname,
        },
        JOB_OPTIONS
      );

      res.json({ jobId: job.id, dataType, fileName: matrixFile.originalname });
    } catch (err) {
      res.status(400).json({ error: err.message });
    }
  }
);

// GET /api/status/:jobId
router.get("/status/:jobId", async (req, res) => {
  const job = await analysisQueue.getJob(req.params.jobId);
  if (!job) return res.status(404).json({ error: "Job not found" });

  const state = await job.getState(); // waiting|active|completed|failed|delayed
  const progress = job.progress || {}; // custom step info the worker reports
  let message = null;
  let workerAvailable = null;
  if (["waiting", "paused", "prioritized"].includes(state)) {
    workerAvailable = (await analysisQueue.getWorkers()).length > 0;
    if (await analysisQueue.isPaused()) {
      message = "The analysis queue is paused.";
    } else if (!workerAvailable) {
      message = "Waiting for an analysis worker. Start the app with npm start, or start a separate worker with npm run worker.";
    } else {
      message = "Queued — waiting for a worker to finish its current analysis.";
    }
  }

  res.json({
    jobId: job.id,
    state,
    step: progress.step || null,
    dataType: progress.dataType || null,
    error: state === "failed" ? job.failedReason : null,
    workerAvailable,
    message,
  });
});

// GET /api/report/:jobId
router.get("/report/:jobId", async (req, res) => {
  const reportPath = path.join(REPORTS_DIR, `${req.params.jobId}.html`);
  if (!fs.existsSync(reportPath)) {
    return res.status(404).json({ error: "Report not ready yet" });
  }
  res.sendFile(reportPath);
});

// Handles multer errors (file too large, bad field name) thrown before the route body runs
router.use((err, req, res, next) => {
  res.status(400).json({ error: err.message });
});

module.exports = router;
