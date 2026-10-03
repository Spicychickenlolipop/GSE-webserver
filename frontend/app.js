"use strict";

const API_BASE = "/api";
const POLL_INTERVAL_MS = 2000;

function switchTab(tab) {
  document.getElementById("tab-gse").classList.toggle("active", tab === "gse");
  document.getElementById("tab-upload").classList.toggle("active", tab === "upload");
  document.getElementById("panel-gse").classList.toggle("active", tab === "gse");
  document.getElementById("panel-upload").classList.toggle("active", tab === "upload");
}

function wireDropzone(zoneId, inputId, primaryTextId, defaultText) {
  const zone = document.getElementById(zoneId);
  const input = document.getElementById(inputId);
  const primary = document.getElementById(primaryTextId);

  function showFile(file) {
    if (file) {
      primary.innerHTML = `<span class="filename">${file.name}</span>`;
      zone.classList.add("has-file");
    } else {
      primary.textContent = defaultText;
      zone.classList.remove("has-file");
    }
  }

  input.addEventListener("change", () => showFile(input.files[0]));

  zone.addEventListener("dragover", (event) => {
    event.preventDefault();
    zone.classList.add("drag-over");
  });
  zone.addEventListener("dragleave", () => zone.classList.remove("drag-over"));
  zone.addEventListener("drop", (e) => {
    e.preventDefault();
    zone.classList.remove("drag-over");
    if (e.dataTransfer.files.length) {
      input.files = e.dataTransfer.files;
      showFile(input.files[0]);
    }
  });
}

function setStatus(text, cls) {
  const card = document.getElementById("status-card");
  const el = document.getElementById("status-text");
  card.classList.add("visible");
  el.textContent = text;
  el.className = cls || "";
}

async function requestJson(endpoint, options, fallbackMessage) {
  const response = await fetch(`${API_BASE}${endpoint}`, options);
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || fallbackMessage);
  return data;
}

async function submitJob(endpoint, options, button, message) {
  button.disabled = true;
  setStatus(message);

  try {
    const data = await requestJson(endpoint, options, "Failed to submit");
    pollStatus(data.jobId, button);
  } catch (error) {
    setStatus(`Error: ${error.message}`, "error");
    button.disabled = false;
  }
}

async function submitGseJob() {
  const gseId = document.getElementById("gse-input").value.trim();
  const button = document.getElementById("submit-btn-gse");
  document.getElementById("report-link").innerHTML = "";
  if (!gseId) return;

  const payload = {
    gseId,
    dataType: document.getElementById("gse-datatype").value || undefined,
    dataScale: document.getElementById("gse-scale").value,
  };
  return submitJob("/analyze", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  }, button, "Submitting job...");
}

async function submitUploadJob() {
  const button = document.getElementById("submit-btn-upload");
  const matrixInput = document.getElementById("matrix-file-input");
  const labelsInput = document.getElementById("labels-file-input");
  document.getElementById("report-link").innerHTML = "";

  if (!matrixInput.files.length) {
    setStatus("Choose a matrix file first.", "error");
    return;
  }

  const formData = new FormData();
  formData.append("matrixFile", matrixInput.files[0]);
  formData.append("dataType", document.getElementById("datatype-select").value);
  formData.append("dataScale", document.getElementById("upload-scale").value);
  if (labelsInput.files.length) {
    formData.append("labelsFile", labelsInput.files[0]);
  }

  return submitJob("/upload", {
    method: "POST",
    body: formData,
  }, button, "Uploading and submitting job...");
}

async function pollStatus(jobId, button) {
  const linkEl = document.getElementById("report-link");

  const poll = async () => {
    try {
      const data = await requestJson(`/status/${jobId}`, {}, "Could not read job status");

      if (data.state === "completed") {
        setStatus(`Done — ${data.dataType || "analysis complete"}`, "done");
        linkEl.innerHTML = `<a class="report-link" href="${API_BASE}/report/${jobId}" target="_blank">Open report →</a>`;
        button.disabled = false;
        return;
      }
      if (data.state === "failed") {
        const error = data.error || "unknown error";
        const summary = error.includes("Traceback (most recent call last)")
          ? error.trim().split(/\r?\n/).filter(Boolean).pop() : error;
        setStatus(`Failed: ${summary}`, "error");
        button.disabled = false;
        return;
      }

      const step = data.step ? ` — ${data.step}` : "";
      const dataType = data.dataType ? ` (${data.dataType})` : "";
      setStatus(data.message || `${data.state}${step}${dataType}`);
    } catch (err) {
      setStatus(`Status unavailable: ${err.message}. Retrying...`, "error");
    }
    setTimeout(poll, POLL_INTERVAL_MS);
  };

  poll();
}

function initializePage() {
  document.getElementById("tab-gse").addEventListener("click", () => switchTab("gse"));
  document.getElementById("tab-upload").addEventListener("click", () => switchTab("upload"));
  document.getElementById("submit-btn-gse").addEventListener("click", submitGseJob);
  document.getElementById("submit-btn-upload").addEventListener("click", submitUploadJob);

  for (const name of ["matrix", "labels"]) {
    const input = document.getElementById(`${name}-file-input`);
    document.getElementById(`${name}-dropzone`).addEventListener("click", () => input.click());
    wireDropzone(
      `${name}-dropzone`,
      `${name}-file-input`,
      `${name}-primary-text`,
      `Drag & drop your ${name} file here`,
    );
  }
}

initializePage();
