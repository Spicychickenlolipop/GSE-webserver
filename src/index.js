const express = require("express");
const cors = require("cors");
const path = require("path");
const analyzeRouter = require("./routes/analyze");

// A default launch must consume the jobs that the API accepts.
// Set RUN_WORKER=false only when running a separate worker service.
if (process.env.RUN_WORKER !== "false") {
  require("./queue/worker");
}

const app = express();
app.use(cors());
app.use(express.json());

app.use("/api", analyzeRouter);

// Serve the frontend itself — visit http://localhost:4000 in the browser
app.use(express.static(path.join(__dirname, "..", "frontend")));

// Serve generated reports as static files
app.use("/reports", express.static(path.join(__dirname, "..", "reports")));

const PORT = process.env.PORT || 4000;
app.listen(PORT, () => {
  console.log(`GSE webserver API listening on http://localhost:${PORT}`);
});
