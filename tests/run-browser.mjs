// Opens a check page of tests/ in headless Chromium and waits until it reports "done.json".
//   node tests/run-browser.mjs render-check.html [--angle=vulkan] [--scale=1.875] [--software] [--log]
// Everything the page reports is saved to tests/out/.
import { spawn } from 'node:child_process';
import { mkdtempSync, rmSync, openSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { startServer } from '../tools/serve.mjs';

const page = process.argv[2], software = process.argv.includes('--software');
// --default-gpu: let Chromium pick its GPU backend, as a normal browser window does; --scale=N sets devicePixelRatio
const defaultGpu = process.argv.includes('--default-gpu');
const scaleArg = process.argv.find((a) => a.startsWith('--scale='));
// --angle=gl|vulkan|...: ANGLE backend (a normal Chromium window on Linux uses gl)
const angle = process.argv.find((a) => a.startsWith('--angle='))?.slice(8) ?? 'gl';
const profile = mkdtempSync(join(tmpdir(), 'chemsim-chrome-'));
let browser;
const finish = (code) => { browser?.kill(); server.close(); server.closeAllConnections(); try { rmSync(profile, { recursive: true, force: true }); } catch {} process.exitCode = code; };
const server = await startServer(0, (name, body) => {
  if (name === 'log.txt') console.log(body.toString());
  if (name === 'done.json') { console.log(body.toString()); finish(JSON.parse(body).ok === false ? 1 : 0); }
});
const flags = ['--headless=new', '--enable-logging=stderr', '--no-first-run', '--window-size=2100,1400', `--user-data-dir=${profile}`,
  ...(scaleArg ? [`--force-device-scale-factor=${scaleArg.slice(8)}`] : []),
  ...(defaultGpu ? [] : software ? ['--disable-gpu', '--enable-unsafe-swiftshader'] : ['--enable-gpu', '--ignore-gpu-blocklist', `--use-angle=${angle}`, ...(angle === 'vulkan' ? ['--enable-features=Vulkan'] : [])])];
browser = spawn('chromium', [...flags, `http://127.0.0.1:${server.address().port}/tests/${page}`], 
  // --log: keep the browser's own log (GPU process messages included) in tests/out/browser.log
  process.argv.includes('--log') ? { stdio: ['ignore', 'ignore', openSync(new URL('./out/browser.log', import.meta.url), 'w')] } : { stdio: 'ignore' });
// the software renderer is very slow: generous timeout
setTimeout(() => { console.log('{"ok":false,"error":"timeout"}'); finish(1); }, 15 * 60e3).unref();
