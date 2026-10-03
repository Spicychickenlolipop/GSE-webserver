const { Queue, QueueEvents } = require("bullmq");
const IORedis = require("ioredis");

const connection = new IORedis(process.env.REDIS_URL || "redis://localhost:6379", {
  maxRetriesPerRequest: null,
});

const analysisQueue = new Queue("gse-analysis", { connection });
const analysisQueueEvents = new QueueEvents("gse-analysis", { connection });

module.exports = { analysisQueue, analysisQueueEvents, connection };
