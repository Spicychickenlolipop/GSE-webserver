const { test, mock } = require('node:test');
const assert = require('node:assert/strict');
const { EventEmitter } = require('node:events');
const { PassThrough } = require('node:stream');
const childProcess = require('node:child_process');

let processJob;
let child;
const bullPath = require.resolve('bullmq');
require.cache[bullPath] = { id: bullPath, filename: bullPath, loaded: true, exports: { Worker: class extends EventEmitter {
  constructor(name, handler) { super(); processJob = handler; }
} } };
mock.method(childProcess, 'spawn', () => {
  child = new EventEmitter();
  child.stdout = new PassThrough();
  child.stderr = new PassThrough();
  return child;
});
const queuePath = require.resolve('./queue');
require.cache[queuePath] = { id: queuePath, filename: queuePath, loaded: true, exports: { connection: {} } };
require('./worker');

function job() {
  const updates = [];
  return { id: 'worker-test', data: { inputType: 'gse', gseId: 'GSE12345' }, updates,
    updateProgress: async value => { await new Promise(resolve => setImmediate(resolve)); updates.push(value); } };
}

async function started(current) {
  const result = processJob(current);
  await new Promise(resolve => setImmediate(resolve));
  return { result };
}

test('worker parses split and combined progress lines in order', async () => {
  const current = job();
  const { result } = await started(current);
  child.stdout.write('PRO');
  child.stdout.write('GRESS:preprocessing:microarray\nPROGRESS:done:microarray\n');
  child.stdout.end();
  child.emit('close', 0);
  await result;
  assert.deepEqual(current.updates.map(value => value.step), ['fetching', 'preprocessing', 'done']);
});

test('worker reports the structured error instead of a Python traceback', async () => {
  const { result } = await started(job());
  const rejected = assert.rejects(result, { message: 'GEO download failed. Please retry.' });
  child.stderr.write('Traceback (most recent call last):\nlarge internal error\n');
  child.stdout.end('ERROR:{"message":"GEO download failed. Please retry."}\n');
  child.emit('close', 1);
  await rejected;
});

test('worker retains a readable fallback for failures without the protocol', async () => {
  const { result } = await started(job());
  const rejected = assert.rejects(result, { message: 'ValueError: invalid matrix' });
  child.stderr.write('Traceback (most recent call last):\nstack frame\nValueError: invalid matrix\n');
  child.stdout.end();
  child.emit('close', 1);
  await rejected;
});
