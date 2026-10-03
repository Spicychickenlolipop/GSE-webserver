const { Worker } = require("bullmq");
const { connection } = require("./queue");
const processJob = require("./process_job");

const worker = new Worker("gse-analysis", processJob, {
  connection,
  concurrency: 2,
});

worker.on("failed", (job, err) => {
  console.error(`Job ${job?.id} failed:`, err.message);
});

worker.on("error", (err) => {
  console.error("Analysis worker error:", err.message);
});

worker.on("ready", () => {
  console.log("GSE analysis worker connected to Redis, ready for jobs.");
});
