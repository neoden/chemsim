// Engine checks that need no browser:  node tests/engines.test.mjs
//  1. the engine describes its basis correctly through the common interface
//     (the density integrates to the number of valence electrons);
//  2. the xTB frontier orbitals have the expected character on the library geometries.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { fixSigns, frontier, describe, shellValues, overlapMatrix } from '../app/src/orbitals.js';
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
  const n = Math.round(2 * L / h), ao = new Float64Array(r.nao);
  let total = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) {
    const p = [-L + (i + 0.5) * h, -L + (j + 0.5) * h, -L + (k + 0.5) * h];
    let at = 0;
    for (const sh of r.shells) at += shellValues(sh, p[0] - m.positions[3 * sh.atom], p[1] - m.positions[3 * sh.atom + 1], p[2] - m.positions[3 * sh.atom + 2], ao, at);
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
  // sulfur brings a d shell
  const s = mol('h2s'), nd = electrons(s, run(engine, s));
  ok(Math.abs(nd - 8) < 2e-3, `${name}, hydrogen sulfide (s, p and d functions): ${nd.toFixed(4)} electrons (expected 8)`);
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
    if (sh.l !== 1) { for (let i = 0; i < 2 * sh.l + 1; i++) rest += c[ao + i] ** 2; ao += 2 * sh.l + 1; continue; }
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
// squared weight of the d functions of atom `atom` in orbital k
function dWeight(r, k, atom) {
  let w = 0, ao = 0;
  for (const sh of r.shells) { const n = 2 * sh.l + 1; if (sh.l === 2 && sh.atom === atom) for (let i = 0; i < n; i++) w += r.coeff[k * r.nao + ao + i] ** 2; ao += n; }
  return w;
}
{
  // a d metal: the five lowest empty orbitals of TiCl4 are the titanium d orbitals, split
  // by the tetrahedron into a pair (e) below a triple (t2)
  const m = mol('ticl4'), r = run(xtb, m), { lumo } = frontier(r);
  const w = [0, 1, 2, 3, 4].map((i) => dWeight(r, lumo + i, 0)), e = (i) => r.eps[lumo + i] * HARTREE_EV;
  ok(Math.min(...w) > 0.5, `the five lowest empty orbitals of TiCl4 are mostly titanium d: weights ${w.map((v) => v.toFixed(2)).join(' ')}`);
  ok(e(1) - e(0) < 0.02 && e(4) - e(2) < 0.02 && e(2) - e(1) > 0.3, `they split into e below t2: ${[0, 1, 2, 3, 4].map((i) => e(i).toFixed(2)).join(' ')} eV`);
}
for (const key of ['sf6', 'ticl4', 'feco5']) {
  const m = mol(key), r = run(xtb, m), S = overlapMatrix(r.shells, m.positions), n = r.nao;
  let err = 0;
  for (let i = 0; i < r.nmo; i++) for (let j = 0; j <= i; j++) {
    let t = 0;
    for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) t += r.coeff[i * n + a] * S[a * n + b] * r.coeff[j * n + b];
    err = Math.max(err, Math.abs(t - (i === j ? 1 : 0)));
  }
  ok(err < 1e-8, `${key}: orbitals with d functions are orthonormal in the analytic overlap matrix (max deviation ${err.toExponential(1)})`);
}
console.log(`${checks} checks passed`);

// ---- appended for stage 3: overlap matrix and orbital continuity
{
  const { followOrbitals } = await import('../app/src/orbitals.js');
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
}
{
  // every element of the table, H to Rn: the engine takes the atom with the number of valence
  // electrons the table gives it and refuses the other parity; its orbitals are orthonormal
  // under the overlap matrix the page computes from the shells
  const { SYMBOLS, VALENCE, MASS_AMU, COVALENT } = await import('../app/src/units.js');
  const wrong = [];
  let count = 0, worstNorm = 0;
  for (let z = 1; z in SYMBOLS; z++, count++) {
    const v = VALENCE[z], atom = (unpaired, orbitals) => xtb.compute({ numbers: [z], positions: [0, 0, 0], charge: 0, etemp: 300, forces: false, orbitals, unpaired });
    const r = atom(v % 2, true);
    if (!r.ok || Math.abs(Array.from(r.occ).reduce((a, b) => a + b, 0) - v) > 1e-6 || atom((v + 1) % 2, false).ok || !(MASS_AMU[z] > 0) || !(COVALENT[z] > 0)) { wrong.push(SYMBOLS[z]); continue; }
    const S = overlapMatrix(r.shells, new Float64Array(3)), n = r.nao;
    for (let i = 0; i < r.nmo; i++) {
      let s = 0;
      for (let a = 0; a < n; a++) for (let b = 0; b < n; b++) s += r.coeff[i * n + a] * S[a * n + b] * r.coeff[i * n + b];
      worstNorm = Math.max(worstNorm, Math.abs(s - 1));
    }
  }
  ok(count === 86 && wrong.length === 0, `all ${count} elements, H to ${SYMBOLS[count]}: the engine counts the valence electrons the table says${wrong.length ? ` (not: ${wrong.join(' ')})` : ''}`);
  ok(worstNorm < 1e-8, `and their orbitals are normalized under the page's own overlap matrix (within ${worstNorm.toExponential(1)})`);
  console.log(`${checks} checks passed in total`);
}
