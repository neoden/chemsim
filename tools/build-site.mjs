// Builds the site for showing:  node tools/build-site.mjs
//   1. app/licenses/  licence texts that must travel with the WebAssembly module: GPL and
//      LGPL with axtb's NOTICE, and the notices of every Rust crate linked into it;
//   2. app/source/    the source of the wrapper and how to rebuild the module, so that anyone
//      can relink it with a changed axtb (LGPL-3.0, section 4);
//   3. dist/          a clean copy of app/: plain static files, relative links only, no
//      server-side requirements. Upload its contents anywhere, into any subfolder.
import { execFileSync } from 'node:child_process';
import { cpSync, rmSync, mkdirSync, readdirSync, writeFileSync, readFileSync, statSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { join } from 'node:path';

const root = fileURLToPath(new URL('../', import.meta.url));
const app = join(root, 'app'), axtb = join(root, 'third_party/axtb'), wrapper = join(root, 'engine/xtb-wasm');
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;');
const page = (title, body) => `<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">\n<title>${esc(title)}</title>\n` +
  `<style>body{font:15px/1.5 system-ui,sans-serif;max-width:820px;margin:32px auto;padding:0 16px;background:#0d0f15;color:#e7e9ee}a{color:#6aa9ff}td,th{padding:3px 14px 3px 0;text-align:left;vertical-align:top}th{color:#9299a8;font-weight:500}code{background:#1b1f2b;padding:1px 5px;border-radius:4px}h1{font-size:22px}h2{font-size:17px;margin-top:28px}</style>\n</head><body>\n${body}\n</body></html>\n`;

// ---- which crates are in the module
const meta = JSON.parse(execFileSync('cargo', ['metadata', '--format-version', '1', '--filter-platform', 'wasm32-wasip1'], { cwd: wrapper, maxBuffer: 1 << 28 }));
const linked = new Set(meta.resolve.nodes.map((n) => n.id));
const crates = meta.packages.filter((p) => linked.has(p.id) && p.name !== 'xtb-wasm').sort((a, b) => a.name.localeCompare(b.name));
const axtbCommit = execFileSync('git', ['rev-parse', 'HEAD'], { cwd: axtb, encoding: 'utf8' }).trim();
const axtbVersion = crates.find((c) => c.name === 'axtb').version;

// ---- licences
const lic = join(app, 'licenses');
rmSync(lic, { recursive: true, force: true });
mkdirSync(join(lic, 'third-party'), { recursive: true });
for (const [from, to] of [['COPYING', 'GPL-3.0.txt'], ['COPYING.LESSER', 'LGPL-3.0.txt'], ['NOTICE', 'axtb-NOTICE.txt'], ['LICENSE', 'axtb-LICENSE.txt']])
  cpSync(join(axtb, from), join(lic, to));
cpSync(join(root, 'LICENSE'), join(lic, 'chemsim-LICENSE.txt'));
// Some MIT crates do not ship the licence text in the package. For those the standard MIT
// text is written out with the authors named in the crate's own manifest.
const MIT = (who) => `MIT License

Copyright (c) ${who}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
`;
const rows = [], missing = [], written = [];
for (const c of crates) {
  const dir = c.manifest_path.replace(/\/Cargo\.toml$/, ''), own = c.name.startsWith('axtb');
  const files = own ? [] : readdirSync(dir).filter((f) => /^(LICEN[SC]E|COPYING|NOTICE|UNLICENSE)/i.test(f));
  const links = [];
  for (const f of files) {
    const name = `${c.name}-${c.version}-${f.replace(/\.(md|txt)$/i, '')}.txt`;
    cpSync(join(dir, f), join(lic, 'third-party', name));
    links.push(`<a href="third-party/${esc(name)}">${esc(f)}</a>`);
  }
  if (!own && !files.length) {
    if (c.license === 'MIT' && (c.authors?.length || c.repository)) {
      const name = `${c.name}-${c.version}-MIT.txt`;
      writeFileSync(join(lic, 'third-party', name), `${c.name} ${c.version}${c.repository ? `, ${c.repository}` : ''}\nThe package does not include a licence file; its manifest declares the MIT licence and these authors.\n\n` + MIT(c.authors?.length ? c.authors.join(', ') : `the ${c.name} authors (${c.repository})`));
      links.push(`<a href="third-party/${esc(name)}">MIT</a> (no file in the package; authors from the manifest)`);
      written.push(c.name);
    } else missing.push(`${c.name} ${c.version} (${c.license})`);
  }
  rows.push(`<tr><td>${esc(c.name)}</td><td>${esc(c.version)}</td><td>${esc(c.license)}</td><td>${own ? '<a href="LGPL-3.0.txt">LGPL</a>, <a href="GPL-3.0.txt">GPL</a>, <a href="axtb-NOTICE.txt">NOTICE</a>' : links.join(', ') || 'the package ships no license text'}</td></tr>`);
}
writeFileSync(join(lic, 'index.html'), page('chemsim: licenses', `<h1>Licenses</h1>
<p><a href="../index.html">← back to the app</a></p>
<h2>The calculation module</h2>
<p>The file <code>engine/xtb_wasm.wasm</code> contains the <a href="https://github.com/atomsforhumanity/axtb">axtb</a> library, version ${esc(axtbVersion)} (the GFN2-xTB method), which is distributed under the
<a href="LGPL-3.0.txt">GNU Lesser General Public License, version 3</a> or later; it supplements the <a href="GPL-3.0.txt">GNU General Public License, version 3</a>.
Where the parameters and algorithms of axtb come from is stated in its <a href="axtb-NOTICE.txt">NOTICE</a>; how the licenses apply to the parts of axtb is in its <a href="axtb-LICENSE.txt">LICENSE</a>.</p>
<p>The library is used unmodified. Its source code: <a href="https://github.com/atomsforhumanity/axtb/tree/${axtbCommit}">github.com/atomsforhumanity/axtb</a>, commit <code>${axtbCommit.slice(0, 12)}</code> (tag v${esc(axtbVersion)}).
The source code of the wrapper it is built into the module with, and instructions for rebuilding with a modified version of the library, are <a href="../source/index.html">next to the site</a>.</p>
<h2>chemsim</h2>
<p>The code of chemsim itself (the application, the shaders and the wrapper of the module) is under the <a href="chemsim-LICENSE.txt">MIT license</a>. The axtb library inside the module stays under the LGPL.</p>
<h2>What the module contains</h2>
<table><tr><th>Package</th><th>Version</th><th>License</th><th>Texts</th></tr>
${rows.join('\n')}
</table>`));

// ---- source of the wrapper, for relinking
const src = join(app, 'source');
rmSync(src, { recursive: true, force: true });
mkdirSync(join(src, 'xtb-wasm'), { recursive: true });
for (const f of ['Cargo.toml', 'Cargo.lock', 'src']) cpSync(join(wrapper, f), join(src, 'xtb-wasm', f), { recursive: true });
cpSync(join(root, 'LICENSE'), join(src, 'xtb-wasm', 'LICENSE'));
cpSync(join(app, 'src/engine/xtb-binding.js'), join(src, 'xtb-binding.js'));
const list = (dir, base = '') => readdirSync(join(dir, base)).flatMap((f) => (statSync(join(dir, base, f)).isDirectory() ? list(dir, join(base, f)) : [join(base, f)]));
writeFileSync(join(src, 'index.html'), page('chemsim: source code of the calculation module', `<h1>Source code of the calculation module</h1>
<p><a href="../index.html">← back to the app</a> · <a href="../licenses/index.html">licenses</a></p>
<p>The module <code>engine/xtb_wasm.wasm</code> is the axtb library, version ${esc(axtbVersion)}, and a small wrapper with a flat C interface. Everything needed to build it again is here, including with a modified version of axtb.</p>
<p>The wrapper is under the <a href="xtb-wasm/LICENSE">MIT license</a>. You may change it and rebuild the module, also with a modified version of axtb.</p>
<h2>Files</h2>
<ul>${list(src).filter((f) => f !== 'index.html').map((f) => `<li><a href="${esc(f)}">${esc(f)}</a></li>`).join('')}</ul>
<h2>Building</h2>
<ol>
<li>Install Rust 1.89 or later and the <code>wasm32-wasip1</code> target: <code>rustup target add wasm32-wasip1</code>.</li>
<li>Get axtb: <code>git clone https://github.com/atomsforhumanity/axtb</code>, commit <code>${axtbCommit}</code> (or your own modified version).</li>
<li>In <code>xtb-wasm/Cargo.toml</code>, set the path to the <code>crates/axtb</code> directory of that clone in the line <code>axtb = { path = "…" }</code>.</li>
<li>In the <code>xtb-wasm</code> directory, run <code>cargo build --release --lib --target wasm32-wasip1</code>.</li>
<li>Replace <code>engine/xtb_wasm.wasm</code> on the site with the file <code>target/wasm32-wasip1/release/xtb_wasm.wasm</code>.</li>
</ol>
<p>The module calls six WASI functions; their implementation and a description of its interface are in <a href="xtb-binding.js">xtb-binding.js</a>.</p>`));

// ---- dist
const dist = join(root, 'dist');
rmSync(dist, { recursive: true, force: true });
cpSync(app, dist, { recursive: true });
if (!existsSync(join(dist, 'engine/xtb_wasm.wasm'))) throw new Error('app/engine/xtb_wasm.wasm is missing: run tools/build-engine.sh first');
const files = list(dist), bytes = files.reduce((s, f) => s + statSync(join(dist, f)).size, 0);
const absolute = files.filter((f) => /\.(html|js|css)$/.test(f) && !f.startsWith('source/')).filter((f) => /(src|href)="\/|fetch\('\/|from '\//.test(readFileSync(join(dist, f), 'utf8')));
console.log(`dist/: ${files.length} files, ${(bytes / 1e6).toFixed(2)} MB (module ${(statSync(join(dist, 'engine/xtb_wasm.wasm')).size / 1e6).toFixed(2)} MB)`);
console.log(`crates in the module: ${crates.length}; licence text written from the manifest for ${written.length} MIT crates (${written.join(', ')}); still without a text: ${missing.length ? missing.join(', ') : 'none'}`);
console.log(`files with absolute links: ${absolute.length ? absolute.join(', ') : 'none'}`);
