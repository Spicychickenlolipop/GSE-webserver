const express = require("express");
const fs = require("fs");
const path = require("path");
const multer = require("multer");
const { v4: uuidv4 } = require("uuid");
const { analysisQueue } = require("../queue/queue");

const router = express.Router();

const GSE_REGEX = /^GSE\d+$/i;
const VALID_DATA_TYPES = ["microarray", "rnaseq", "methylation"];

// --- File upload setup ---
const UPLOADS_DIR = path.join(__dirname, "..", "..", "uploads");
if (!fs.existsSync(UPLOADS_DIR)) fs.mkdirSync(UPLOADS_DIR, { recursive: true });

const ALLOWED_EXTENSIONS = new Set([".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz"]);
const inputExtension = (name) => name.toLowerCase().match(/\.(?:csv|tsv|txt)(?:\.gz)?$/)?.[0];

const storage = multer.diskStorage({
  destination: (req, file, cb) => cb(null, UPLOADS_DIR),
  filename: (req, file, cb) => {
    const ext = inputExtension(file.originalname);
    cb(null, `${uuidv4()}${ext}`);
  },
});

const upload = multer({
  storage,
  limits: { fileSize: 500 * 1024 * 1024 }, // 500MB — expression/methylation matrices can be large
  fileFilter: (req, file, cb) => {
    const ext = inputExtension(file.originalname);
    if (!ALLOWED_EXTENSIONS.has(ext)) {
      return cb(new Error('Use .csv, .tsv, or .txt, optionally compressed as .gz.'));
    }
    cb(null, true);
  },
});

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
    {
      removeOnComplete: false,
      removeOnFail: false,
      attempts: 1,
    }
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
        {
          removeOnComplete: false,
          removeOnFail: false,
          attempts: 1,
        }
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
  const reportPath = path.join(__dirname, "..", "..", "reports", `${req.params.jobId}.html`);
  if (!fs.existsSync(reportPath)) {
    return res.status(404).json({ error: "Report not ready yet" });
  }
  res.sendFile(reportPath);
});

// Handles multer errors (file too large, bad field name) thrown before the route body runs
router.use((err, req, res, next) => {
  if (err instanceof multer.MulterError || err) {
    return res.status(400).json({ error: err.message });
  }
  next();
});

module.exports = router;
