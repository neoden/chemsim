// Whole-page screenshot of the app:  node tests/screenshot.mjs "<query>" <out.png> [WxH] [--angle=vulkan]
// Drives headless Chromium over the DevTools protocol and waits until the app reports that
// it is ready (Chromium's own --screenshot fires at the load event, before the first result).
import { spawn } from 'node:child_process';
import { mkdtempSync, rmSync, readFileSync, writeFileSync, existsSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join, resolve } from 'node:path';
import { startServer } from '../tools/serve.mjs';

const args = process.argv.slice(2).filter((a) => !a.startsWith('--'));
const [query = '', out = 'tests/out/app.png', size = '1280x800'] = args;
const angle = process.argv.find((a) => a.startsWith('--angle='))?.slice(8) ?? 'gl';
const profile = mkdtempSync(join(tmpdir(), 'chemsim-chrome-'));
const server = await startServer(0);
const browser = spawn('chromium', ['--headless=new', '--no-first-run', '--enable-gpu', '--ignore-gpu-blocklist', `--use-angle=${angle}`,
  ...(angle === 'vulkan' ? ['--enable-features=Vulkan'] : []), `--user-data-dir=${profile}`, `--window-size=${size.replace('x', ',')}`,
  '--remote-debugging-port=0', 'about:blank'], { stdio: 'ignore' });
try {
  const portFile = join(profile, 'DevToolsActivePort');
  for (let i = 0; i < 200 && !existsSync(portFile); i++) await new Promise((r) => setTimeout(r, 50));
  const port = readFileSync(portFile, 'utf8').split('\n')[0];
  const targets = await (await fetch(`http://127.0.0.1:${port}/json`)).json();
  const ws = new WebSocket(targets.find((t) => t.type === 'page').webSocketDebuggerUrl);
  await new Promise((r, j) => { ws.onopen = r; ws.onerror = j; });
  let id = 0;
  const pending = new Map();
  ws.onmessage = ({ data }) => { const m = JSON.parse(data); pending.get(m.id)?.(m); pending.delete(m.id); };
  const send = (method, params = {}) => new Promise((r) => { pending.set(++id, r); ws.send(JSON.stringify({ id, method, params })); });
  await send('Page.enable');
  await send('Page.navigate', { url: `http://127.0.0.1:${server.address().port}/app/index.html?${query}` });
  const ready = await send('Runtime.evaluate', {
    awaitPromise: true, returnByValue: true,
    expression: `(async () => { for (let i = 0; i < 400 && !window.chemsim; i++) await new Promise((r) => setTimeout(r, 25));
      await Promise.race([window.chemsim.ready, new Promise((r) => setTimeout(r, 20000))]);
      await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
      return document.getElementById('status').textContent; })()`,
  });
  if (ready.result?.exceptionDetails) console.error('page error:', ready.result.exceptionDetails.exception?.description);
  const shot = await send('Page.captureScreenshot', { format: 'png' });
  writeFileSync(resolve(out), Buffer.from(shot.result.data, 'base64'));
  ws.close();
} finally {
  browser.kill();
  server.close(); server.closeAllConnections();
  await new Promise((r) => setTimeout(r, 200));
  try { rmSync(profile, { recursive: true, force: true }); } catch {}
}
