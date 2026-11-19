// Serves the probe page with plain static headers (no COOP/COEP), launches headless
// Chromium on it and prints what the page reports.
import { createServer } from 'node:http';
import { readFileSync } from 'node:fs';
import { spawn } from 'node:child_process';
const root = new URL('../', import.meta.url);
const files = {
  '/': ['probe/bench.html', 'text/html'], '/bench-worker.mjs': ['probe/bench-worker.mjs', 'text/javascript'],
  '/xtb.mjs': ['probe/xtb.mjs', 'text/javascript'],
  '/app/src/engine/xtb-binding.js': ['../app/src/engine/xtb-binding.js', 'text/javascript'],
  '/xtb_wasm.wasm': ['xtb-wasm/target/wasm32-wasip1/release/xtb_wasm.wasm', 'application/wasm'],
};
for (const c of ['h2o', 'bh3_nh3', 'hexane']) files[`/cases/${c}.xyz`] = [`cases/${c}.xyz`, 'text/plain'];
let browser;
const server = createServer((req, res) => {
  if (req.method === 'POST') {
    let body = '';
    req.on('data', (d) => (body += d)).on('end', () => {
      res.end('ok');
      console.log(body);
      browser.kill();
      server.close();
    });
    return;
  }
  const f = files[req.url];
  if (!f) return res.writeHead(404).end();
  res.writeHead(200, { 'content-type': f[1] }).end(readFileSync(new URL(f[0], root)));
}).listen(0, '127.0.0.1', () => {
  const dir = process.argv[2];
  browser = spawn('chromium', ['--headless=new', '--no-first-run', '--disable-gpu', `--user-data-dir=${dir}`,
    `http://127.0.0.1:${server.address().port}/`], { stdio: 'ignore' });
  setTimeout(() => { console.log('{"error":"timeout"}'); browser.kill(); server.close(); }, 120000).unref();
});
