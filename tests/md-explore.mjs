// Exploration for stage 3: energy drift against SCF accuracy and time step.
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { Dynamics } from '../app/src/md.js';
import { BOHR_ANGSTROM, NUMBERS } from '../app/src/units.js';
const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
export function setup(key, { accuracy = 1, etemp = 300, seed = 1, wallRadius = Infinity } = {}) {
  const m = library.find((x) => x.key === key), numbers = m.elements.map((el) => NUMBERS[el]);
  let iterations = 0;
  const force = (x, orbitals) => { const r = xtb.compute({ numbers, positions: x, charge: 0, etemp, forces: true, orbitals, accuracy }); if (r.ok) iterations += r.iterations; return r; };
  const md = new Dynamics({ numbers, positions: Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM), force, seed, wallRadius });
  return { md, m, iterations: () => iterations };
}
// slope of a least-squares line through (t, y)
export function slope(t, y) {
  const n = t.length, mt = t.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let a = 0, b = 0;
  for (let i = 0; i < n; i++) { a += (t[i] - mt) * (y[i] - my); b += (t[i] - mt) ** 2; }
  return a / b;
}
if (process.argv[1].endsWith('md-explore.mjs')) {
  for (const key of ['h2o', 'bh3nh3']) for (const accuracy of [1, 0.1, 0.01]) for (const dt of [0.5, 0.25]) {
    const { md, iterations } = setup(key, { accuracy });
    md.params.thermostat = false; md.params.dt = dt;
    md.init(false); md.thermalize(300);
    const t = [], e = [], ps = 2, t0 = performance.now();
    for (let i = 0; i < ps * 1000 / dt; i++) { if (!md.step()) throw new Error(md.error.error); t.push(md.time / 1000); e.push(md.etot); }
    const ms = performance.now() - t0, drift = slope(t, e) / md.n;
    const mean = e.reduce((a, b) => a + b) / e.length, rms = Math.sqrt(e.reduce((a, b) => a + (b - mean) ** 2, 0) / e.length);
    console.log(`${key.padEnd(7)} acc ${String(accuracy).padEnd(5)} dt ${dt}  drift ${drift.toExponential(2)} Eh/atom/ps  rms fluct ${rms.toExponential(1)} Eh  ` +
      `${(ms / md.steps).toFixed(3)} ms/step  ${(iterations() / md.steps).toFixed(1)} SCF it/step`);
  }
}
