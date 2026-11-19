// Static server for development and for the browser checks.
//   node tools/serve.mjs [port]      serves the project root; the app is at /app/
// POST /__report/<name> saves the body to tests/out/<name> (png data URLs are decoded).
import { createServer } from 'node:http';
import { readFile, writeFile, mkdir } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { join, normalize, extname } from 'node:path';

const root = fileURLToPath(new URL('../', import.meta.url));
const TYPES = { '.html': 'text/html; charset=utf-8', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.json': 'application/json', '.wasm': 'application/wasm', '.png': 'image/png', '.xyz': 'text/plain' };

export function startServer(port = 0, onReport = () => {}) {
  const server = createServer(async (req, res) => {
    const url = new URL(req.url, 'http://x');
    if (req.method === 'POST' && url.pathname.startsWith('/__report/')) {
      const chunks = [];
      for await (const c of req) chunks.push(c);
      let body = Buffer.concat(chunks);
      const name = url.pathname.slice('/__report/'.length).replace(/[^\w.-]/g, '_');
      const text = body.toString('latin1', 0, 30);
      if (text.startsWith('data:image/png;base64,')) body = Buffer.from(body.toString('latin1').slice(22), 'base64');
      await mkdir(join(root, 'tests/out'), { recursive: true });
      await writeFile(join(root, 'tests/out', name), body);
      res.end('ok');
      onReport(name, body);
      return;
    }
    let path = normalize(decodeURIComponent(url.pathname));
    if (path.endsWith('/')) path += 'index.html';
    try {
      const data = await readFile(join(root, path));
      res.writeHead(200, { 'content-type': TYPES[extname(path)] || 'application/octet-stream', 'cache-control': 'no-store' }).end(data);
    } catch {
      res.writeHead(404).end('not found');
    }
  });
  return new Promise((resolve) => server.listen(port, '127.0.0.1', () => resolve(server)));
}

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  const server = await startServer(+(process.argv[2] ?? 8080));
  console.log(`http://127.0.0.1:${server.address().port}/app/`);
}
