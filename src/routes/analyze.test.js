const { test, after } = require('node:test');
const assert = require('node:assert/strict');
const express = require('express');
const fs = require('node:fs');
const path = require('node:path');

const jobs = [];
let jobState = 'waiting';
let connectedWorkers = [];
let queuePaused = false;
const queuePath = require.resolve('../queue/queue');
require.cache[queuePath] = { id: queuePath, filename: queuePath, loaded: true,
  exports: { analysisQueue: {
    add: async (name, data) => { jobs.push(data); return { id: 'test' }; },
    getJob: async () => ({ id: 'test', progress: {}, getState: async () => jobState }),
    getWorkers: async () => connectedWorkers,
    isPaused: async () => queuePaused,
  } } };
const app = express();
app.use(express.json());
app.use('/api', require('./analyze'));
const server = app.listen(0, '127.0.0.1');
const ready = new Promise(resolve => server.on('listening', resolve));
after(async () => {
  for (const job of jobs) {
    for (const file of [job.matrixPath, job.labelsPath].filter(Boolean)) {
      assert.equal(path.dirname(file), path.resolve(__dirname, '../../uploads'));
      fs.unlinkSync(file);
    }
  }
  await new Promise(resolve => server.close(resolve));
});

async function request(endpoint, body) {
  await ready;
  return fetch(`http://127.0.0.1:${server.address().port}/api/${endpoint}`, {
    method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
  });
}

test('invalid GEO request returns validation error', async () => {
  assert.equal((await request('analyze', { gseId: 123 })).status, 400);
  assert.equal((await request('analyze', { gseId: 'GSE123', dataScale: 'invalid' })).status, 400);
});

test('GEO type and scale override reach the job queue', async () => {
  const response = await request('analyze', { gseId: ' gse123 ', dataType: 'rnaseq', dataScale: 'processed' });
  assert.equal(response.status, 200);
  assert.equal(jobs.at(-1).gseId, 'GSE123');
  assert.equal(jobs.at(-1).dataScale, 'processed');
  assert.equal(jobs.at(-1).dataType, 'rnaseq');
});

test('gzip uploads preserve compression suffix and scale', async () => {
  await ready;
  const form = new FormData();
  form.append('matrixFile', new Blob([require('node:zlib').gzipSync('gene,s1\na,1\n')]), 'counts.csv.gz');
  form.append('dataType', 'rnaseq');
  form.append('dataScale', 'processed');
  const response = await fetch(`http://127.0.0.1:${server.address().port}/api/upload`, { method: 'POST', body: form });
  assert.equal(response.status, 200);
  assert.ok(jobs.at(-1).matrixPath.endsWith('.csv.gz'));
  assert.equal(jobs.at(-1).dataScale, 'processed');
});

async function status() {
  await ready;
  const response = await fetch(`http://127.0.0.1:${server.address().port}/api/status/test`);
  assert.equal(response.status, 200);
  return response.json();
}

test('waiting status explains when no worker is connected', async () => {
  const data = await status();
  assert.equal(data.state, 'waiting');
  assert.equal(data.workerAvailable, false);
  assert.match(data.message, /npm start/);
});

test('waiting status distinguishes a busy worker and a paused queue', async () => {
  connectedWorkers = [{ id: 'worker' }];
  assert.match((await status()).message, /finish its current analysis/);
  queuePaused = true;
  assert.match((await status()).message, /queue is paused/);
  queuePaused = false;
});

test('completed status does not report a worker warning', async () => {
  jobState = 'completed';
  connectedWorkers = [];
  const data = await status();
  assert.equal(data.state, 'completed');
  assert.equal(data.message, null);
});
