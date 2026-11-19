// Engine checks that need no browser:  node tests/engines.test.mjs
//  1. the engine describes its basis correctly through the common interface
//     (the density integrates to the number of valence electrons);
//  2. the xTB frontier orbitals have the expected character on the library geometries.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { fixSigns, frontier, describe } from '../app/src/orbitals.js';
import { BOHR_ANGSTROM, HARTREE_EV, NUMBERS } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
const mol = (key) => {
  const m = library.find((x) => x.key === key);
  return { ...m, numbers: m.elements.map((el) => NUMBERS[el]), positions: Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM) };
};
const run = (engine, m) => {
  const r = engine.compute({ numbers: m.numbers, positions: m.positions, charge: 0, etemp: 300, forces: true });
  assert.ok(r.ok, r.error);
  fixSigns(r);
  return r;
};
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };

console.log('basis description through the interface: electrons on a grid');
function electrons(m, r, h = 0.14, L = 8) {
  const aos = [];
  for (const sh of r.shells) for (const k of sh.l ? [0, 1, 2] : [-1]) aos.push({ sh, k, c: m.positions.subarray(3 * sh.atom, 3 * sh.atom + 3) });
  const n = Math.round(2 * L / h), ao = new Float64Array(aos.length);
  let total = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) {
    const p = [-L + (i + 0.5) * h, -L + (j + 0.5) * h, -L + (k + 0.5) * h];
    aos.forEach(({ sh, k, c }, q) => {
      const d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]], r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      let v = 0;
      for (let t = 0; t < sh.alpha.length; t++) v += sh.coeff[t] * Math.exp(-sh.alpha[t] * r2);
      ao[q] = k < 0 ? v : v * d[k];
    });
    for (let o = 0; o < r.nmo; o++) {
      if (r.occ[o] < 1e-12) continue;
      let psi = 0;
      for (let q = 0; q < r.nao; q++) psi += r.coeff[o * r.nao + q] * ao[q];
      total += r.occ[o] * psi * psi;
    }
  }
  return total * h ** 3;
}
for (const [name, engine] of [['xtb', xtb]]) {
  const m = mol('h2o'), ne = electrons(m, run(engine, m));
  ok(Math.abs(ne - 8) < 2e-3, `${name}, water: ${ne.toFixed(4)} electrons (expected 8)`);
}

console.log('character of the xTB frontier orbitals');
const eV = (r, m) => (r.eps[m] * HARTREE_EV).toFixed(2);
const show = (r, m, k) => describe(r, m.elements, k, 0.2).map(([c, l]) => `${l}:${c.toFixed(2)}`).join(' ');
// unit normal of the plane through atoms i, j, k
function normal(p, i, j, k) {
  const a = [0, 1, 2].map((d) => p[3 * j + d] - p[3 * i + d]), b = [0, 1, 2].map((d) => p[3 * k + d] - p[3 * i + d]);
  const n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], l = Math.hypot(...n);
  return n.map((v) => v / l);
}
// for orbital k: the p coefficient vector on each atom projected on `axis`, and the squared
// weight of everything else
function piAnalysis(r, k, axis) {
  const c = r.coeff.subarray(k * r.nao, (k + 1) * r.nao), along = {};
  let rest = 0, ao = 0;
  for (const sh of r.shells) {
    if (sh.l === 0) { rest += c[ao] ** 2; ao += 1; continue; }
    const v = [c[ao], c[ao + 1], c[ao + 2]], a = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
    along[sh.atom] = a;
    rest += v[0] ** 2 + v[1] ** 2 + v[2] ** 2 - a * a;
    ao += 3;
  }
  return { along, rest };
}
{
  const m = mol('h2o'), r = run(xtb, m), { homo } = frontier(r), a = piAnalysis(r, homo, normal(m.positions, 0, 1, 2));
  ok(Math.abs(a.along[0]) > 0.99 && a.rest < 1e-6, `water HOMO (${eV(r, homo)} eV) is the oxygen p orbital perpendicular to the molecular plane: ${show(r, m, homo)}`);
}
{
  const m = mol('bh3'), r = run(xtb, m), { lumo } = frontier(r), a = piAnalysis(r, lumo, normal(m.positions, 1, 2, 3));
  ok(Math.abs(a.along[0]) > 0.99 && a.rest < 1e-6 && r.occ[lumo] < 1e-6, `borane LUMO (${eV(r, lumo)} eV) is the empty boron p orbital perpendicular to the plane: ${show(r, m, lumo)}`);
}
{
  const m = mol('c2h4'), r = run(xtb, m), { homo, lumo } = frontier(r), axis = normal(m.positions, 0, 2, 3);
  const h = piAnalysis(r, homo, axis), l = piAnalysis(r, lumo, axis);
  ok(h.rest < 1e-6 && h.along[0] * h.along[1] > 0.2, `ethylene HOMO (${eV(r, homo)} eV) is pi, both carbon p in phase: ${show(r, m, homo)}`);
  ok(l.rest < 1e-6 && l.along[0] * l.along[1] < -0.2, `ethylene LUMO (${eV(r, lumo)} eV) is pi*, carbon p in antiphase: ${show(r, m, lumo)}`);
}
{
  const m = mol('bh3nh3'), r = run(xtb, m);
  const f = Math.max(...r.forces.map(Math.abs));
  ok(f < 5e-5, `library geometry of BH3NH3 is a stationary point of the engine: max force ${f.toExponential(1)} Eh/bohr`);
}
console.log(`${checks} checks passed`);

// ---- appended for stage 3: overlap matrix and orbital continuity
{
  const { overlapMatrix, followOrbitals } = await import('../app/src/orbitals.js');
  console.log('overlap matrix and orbital continuity');
  for (const [name, engine] of [['xtb', xtb]]) {
    const m = mol('bh3nh3'), r = run(engine, m), S = overlapMatrix(r.shells, m.positions), n = r.nao;
    let err = 0;
    for (let i = 0; i < r.nmo; i++) for (let j = 0; j < r.nmo; j++) {
      let t = 0;
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) t += r.coeff[i * n + a] * S[a * n + b] * r.coeff[j * n + b];
      err = Math.max(err, Math.abs(t - (i === j ? 1 : 0)));
    }
    ok(err < 1e-8, `${name}: orbitals are orthonormal in the analytic overlap matrix (max deviation ${err.toExponential(1)})`);
  }
  // water with one O-H stretched a little: orbitals of the new geometry, with scrambled
  // signs, must come back aligned with the old ones
  const m = mol('h2o'), a = run(xtb, m), moved = { ...m, positions: Float64Array.from(m.positions) };
  moved.positions[3] += 0.05; moved.positions[7] -= 0.03;
  const b = xtb.compute({ numbers: m.numbers, positions: moved.positions, charge: 0, etemp: 300 });
  for (let j = 0; j < b.nmo; j += 2) for (let i = 0; i < b.nao; i++) b.coeff[j * b.nao + i] *= -1;
  const S = overlapMatrix(b.shells, moved.positions), sel = followOrbitals(a, b, S, 3);
  let worst = 1;
  for (let j = 0; j < b.nmo; j++) {
    let t = 0;
    for (let p = 0; p < b.nao; p++) for (let q = 0; q < b.nao; q++) t += a.coeff[j * b.nao + p] * S[p * b.nao + q] * b.coeff[j * b.nao + q];
    worst = Math.min(worst, t);
  }
  ok(worst > 0.9 && sel === 3, `after a small move every orbital keeps its sign and index (smallest overlap with its predecessor ${worst.toFixed(3)})`);
  console.log(`${checks} checks passed in total`);
}
