// Compares the wasm build (run here in node) with the native build of the same crate.
import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { loadXtb, parseXyz, bench, unpack } from './xtb.mjs';
const root = new URL('../', import.meta.url);
const wasm = readFileSync(new URL('xtb-wasm/target/wasm32-wasip1/release/xtb_wasm.wasm', root));
const xtb = await loadXtb(wasm);
const maxdiff = (a, b) => a.reduce((m, v, i) => Math.max(m, Math.abs(v - b[i])), 0);
for (const name of ['h2o', 'bh3_nh3', 'hexane']) {
  const file = new URL(`cases/${name}.xyz`, root);
  const w = bench(xtb, parseXyz(readFileSync(file, 'utf8')));
  const n = JSON.parse(execFileSync(new URL('xtb-wasm/target/release/native', root).pathname, [file.pathname]));
  for (const key of ['out', 'last']) {
    const a = unpack(Float64Array.from(w[key])), b = unpack(Float64Array.from(n[key]));
    // orbital sign is arbitrary: compare coefficients up to a sign per orbital
    let dc = 0;
    for (let k = 0; k < a.nmo; k++) {
      const ca = a.coeff.subarray(k * a.nao, (k + 1) * a.nao), cb = b.coeff.subarray(k * a.nao, (k + 1) * a.nao);
      dc = Math.max(dc, Math.min(maxdiff(ca, cb), maxdiff(ca, cb.map((v) => -v))));
    }
    console.log(`${name.padEnd(8)} ${key.padEnd(4)} E=${a.energy.toFixed(10)} Eh  |dE|=${Math.abs(a.energy - b.energy).toExponential(1)}` +
      `  max|dgrad|=${maxdiff(a.gradient, b.gradient).toExponential(1)}  max|deps|=${maxdiff(a.eps, b.eps).toExponential(1)}` +
      `  max|dC|=${dc.toExponential(1)}  nao=${a.nao} iter=${a.iterations}/${b.iterations}`);
  }
  console.log(`${name.padEnd(8)} time, ms: wasm(node) cold ${w.cold_ms.toFixed(2)} warm ${w.warm_median_ms.toFixed(2)} | native cold ${n.cold_ms.toFixed(2)} warm ${n.warm_median_ms.toFixed(2)}`);
}
