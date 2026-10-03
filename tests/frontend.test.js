const { test } = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const frontend = path.resolve(__dirname, "../frontend");
const html = fs.readFileSync(path.join(frontend, "index.html"), "utf8");
const script = fs.readFileSync(path.join(frontend, "app.js"), "utf8");

function loadPage(responses) {
  const elements = new Map();
  const requests = [];
  const timers = [];
  for (const [, id] of html.matchAll(/id="([^"]+)"/g)) {
    const classes = new Set();
    elements.set(id, {
      value: "", files: [], disabled: false, innerHTML: "", textContent: "",
      listeners: {},
      addEventListener(event, handler) { this.listeners[event] = handler; },
      click() { return this.listeners.click?.(); },
      classList: {
        add(name) { classes.add(name); },
        remove(name) { classes.delete(name); },
        toggle(name, enabled) { enabled ? classes.add(name) : classes.delete(name); },
        contains(name) { return classes.has(name); },
      },
    });
  }
  vm.runInNewContext(script, {
    document: { getElementById: id => elements.get(id) },
    FormData,
    fetch: async (url, options) => {
      requests.push({ url, options });
      const response = responses.shift();
      assert.ok(response, "Unexpected request");
      return { ok: response.ok !== false, json: async () => response.body };
    },
    setTimeout: callback => timers.push(callback),
  });
  return { elements, requests, timers };
}

const flush = () => new Promise(resolve => setImmediate(resolve));

test("page loads separate CSS and deferred JavaScript assets", () => {
  assert.match(html, /href="styles.css"/);
  assert.match(html, /src="app.js" defer/);
  assert.ok(fs.statSync(path.join(frontend, "styles.css")).size > 0);
  assert.doesNotMatch(html, /onclick=|<style>/);
});

test("GEO button submits settings and shows the completed report", async () => {
  const page = loadPage([
    { body: { jobId: "42" } },
    { body: { state: "completed", dataType: "microarray" } },
  ]);
  page.elements.get("gse-input").value = " GSE12345 ";
  page.elements.get("gse-scale").value = "processed";
  await page.elements.get("submit-btn-gse").click();
  await flush();
  assert.equal(page.requests[0].url, "/api/analyze");
  assert.deepEqual(JSON.parse(page.requests[0].options.body), {
    gseId: "GSE12345", dataScale: "processed",
  });
  assert.equal(page.requests[1].url, "/api/status/42");
  assert.match(page.elements.get("report-link").innerHTML, /\/api\/report\/42/);
  assert.equal(page.elements.get("submit-btn-gse").disabled, false);
  assert.equal(page.timers.length, 0);
});

test("upload button submits matrix and labels and keeps polling queued jobs", async () => {
  const page = loadPage([
    { body: { jobId: "43" } },
    { body: { state: "waiting", message: "Queued" } },
  ]);
  page.elements.get("matrix-file-input").files = [new Blob(["gene,sample\na,1\n"])];
  page.elements.get("labels-file-input").files = [new Blob(["sample,group\nsample,control\n"])];
  page.elements.get("datatype-select").value = "rnaseq";
  page.elements.get("upload-scale").value = "auto";
  await page.elements.get("submit-btn-upload").click();
  await flush();
  assert.equal(page.requests[0].url, "/api/upload");
  const form = page.requests[0].options.body;
  assert.equal(form.get("dataType"), "rnaseq");
  assert.ok(form.get("matrixFile"));
  assert.ok(form.get("labelsFile"));
  assert.equal(page.elements.get("status-text").textContent, "Queued");
  assert.equal(page.timers.length, 1);
});

test("submission errors restore the Analyze button", async () => {
  const page = loadPage([{ ok: false, body: { error: "Invalid accession" } }]);
  page.elements.get("gse-input").value = "invalid";
  await page.elements.get("submit-btn-gse").click();
  assert.match(page.elements.get("status-text").textContent, /Invalid accession/);
  assert.equal(page.elements.get("submit-btn-gse").disabled, false);
  assert.equal(page.timers.length, 0);
});
