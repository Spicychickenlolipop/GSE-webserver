const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const readline = require("readline");
const { REPORTS_DIR, PYTHON_BIN, PIPELINE_ENTRY } = require("../config");

const STDERR_LIMIT = 64 * 1024;

fs.mkdirSync(REPORTS_DIR, { recursive: true });

function buildArguments(job) {
  const { inputType, gseId, dataType, dataScale = "auto", matrixPath, labelsPath } = job.data;
  const jobId = job.id;

  const args = [
    PIPELINE_ENTRY,
    "--job-id", jobId,
    "--out", path.join(REPORTS_DIR, `${jobId}.html`),
    "--data-scale", dataScale,
  ];

  if (inputType === "upload") {
    args.push("--matrix-file", matrixPath, "--data-type", dataType);
    if (labelsPath) args.push("--labels-file", labelsPath);
  } else {
    args.push("--gse", gseId);
    if (dataType) args.push("--data-type", dataType);
  }

  return args;
}

function runPython(job, args) {
  const jobId = job.id;
  return new Promise((resolve, reject) => {
    const proc = spawn(PYTHON_BIN, args);
    proc.on("error", (err) => reject(new Error(`Could not start Python: ${err.message}. Check PYTHON_BIN.`)));

    let stderr = "";
    let errorMessage = null;
    let progressUpdates = Promise.resolve();

    readline.createInterface({ input: proc.stdout }).on("line", (line) => {
      // A stream chunk can contain half a line or several complete messages.
      console.log(`[job ${jobId}]`, line);
      if (line.startsWith("ERROR:")) {
        try {
          errorMessage = JSON.parse(line.slice(6)).message;
        } catch {
          // If the message is incomplete, use the stderr summary instead.
        }
      }
      const match = line.match(/^PROGRESS:(\w+)(?::(\w+))?/);
      if (match) {
        progressUpdates = progressUpdates
          .then(() => job.updateProgress({ step: match[1], dataType: match[2] || null }))
          .catch((err) => console.error(`[job ${jobId}] Could not update progress:`, err.message));
      }
    });

    proc.stderr.on("data", (chunk) => {
      stderr = (stderr + chunk.toString()).slice(-STDERR_LIMIT);
      console.error(`[job ${jobId}][stderr]`, chunk.toString());
    });

    proc.on("close", async (code) => {
      await progressUpdates;
      if (code === 0) {
        resolve({ reportPath: `${jobId}.html` });
      } else {
        const lastError = stderr.trim().split(/\r?\n/).filter(Boolean).pop();
        reject(new Error(errorMessage || lastError || `Pipeline exited with code ${code}`));
      }
    });
  });
}

async function processJob(job) {
  const isUpload = job.data.inputType === "upload";
  await job.updateProgress({
    step: isUpload ? "loading_upload" : "fetching",
    dataType: isUpload ? job.data.dataType : null,
  });
  return runPython(job, buildArguments(job));
}

module.exports = processJob;
