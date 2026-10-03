const { Worker } = require("bullmq");
const { spawn } = require("child_process");
const path = require("path");
const fs = require("fs");
const readline = require("readline");
const { connection } = require("./queue");

const REPORTS_DIR = path.join(__dirname, "..", "..", "reports");
if (!fs.existsSync(REPORTS_DIR)) fs.mkdirSync(REPORTS_DIR, { recursive: true });

const PYTHON_BIN = process.env.PYTHON_BIN || (process.platform === "win32" ? "python" : "python3");
const PIPELINE_ENTRY = path.join(__dirname, "..", "..", "python", "run_pipeline.py");

const worker = new Worker(
  "gse-analysis",
  async (job) => {
    const { inputType, gseId, dataType, dataScale = "auto", matrixPath, labelsPath } = job.data;
    const jobId = job.id;

    const args = [
      PIPELINE_ENTRY,
      "--job-id", jobId,
      "--out", path.join(REPORTS_DIR, `${jobId}.html`),
      "--data-scale", dataScale,
    ];

    if (inputType === "upload") {
      await job.updateProgress({ step: "loading_upload", dataType });
      args.push("--matrix-file", matrixPath, "--data-type", dataType);
      if (labelsPath) args.push("--labels-file", labelsPath);
    } else {
      await job.updateProgress({ step: "fetching", dataType: null });
      args.push("--gse", gseId);
      if (dataType) args.push("--data-type", dataType);
    }

    return new Promise((resolve, reject) => {
      // run_pipeline.py handles: [fetch GEO | load upload] -> detect/confirm type -> preprocess -> ML -> report
      const proc = spawn(PYTHON_BIN, args);
      proc.on("error", (err) => reject(new Error(`Could not start Python: ${err.message}. Check PYTHON_BIN.`)));

      let stderr = "";
      let errorMessage = null;
      let progressUpdates = Promise.resolve();

      readline.createInterface({ input: proc.stdout }).on("line", (line) => {
        // A stream chunk can contain half a line or several complete messages.
        console.log(`[job ${jobId}]`, line);
        if (line.startsWith("ERROR:")) {
          try { errorMessage = JSON.parse(line.slice(6)).message; } catch (_) { /* Keep stderr fallback. */ }
        }
        const match = line.match(/^PROGRESS:(\w+)(?::(\w+))?/);
        if (match) {
          progressUpdates = progressUpdates
            .then(() => job.updateProgress({ step: match[1], dataType: match[2] || null }))
            .catch((err) => console.error(`[job ${jobId}] Could not update progress:`, err.message));
        }
      });

      proc.stderr.on("data", (chunk) => {
        stderr = (stderr + chunk.toString()).slice(-65536);
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
  },
  { connection, concurrency: 2 }
);

worker.on("failed", (job, err) => {
  console.error(`Job ${job?.id} failed:`, err.message);
});

worker.on("error", (err) => {
  console.error("Analysis worker error:", err.message);
});

worker.on("ready", () => {
  console.log("GSE analysis worker connected to Redis, ready for jobs.");
});
