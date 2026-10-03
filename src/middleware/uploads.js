const fs = require("fs");
const multer = require("multer");
const { v4: uuidv4 } = require("uuid");
const { UPLOADS_DIR } = require("../config");

fs.mkdirSync(UPLOADS_DIR, { recursive: true });

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

module.exports = upload;
