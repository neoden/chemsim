// The electrostatic potential (app/src/esp.js), without a browser:  node tests/esp.test.mjs
//  1. charge sites with none between the atoms are the engine's own atomic multipoles
//     (two independent sets of integrals), and the sites between the atoms change neither
//     the total charge nor the dipole moment;
//  2. on the surface the app colours, the potential of the sites is the potential of the
//     valence density itself (integrated on a grid) plus the cores;
//  3. it has the signs chemistry expects, also where the atomic multipoles alone get them wrong.
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { shellValues } from '../app/src/orbitals.js';
import { chargeSites, potential, potentialScale, ESP_DENSITY } from '../app/src/esp.js';
import { BOHR_ANGSTROM, NUMBERS, VALENCE } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
function solve(key) {
  const m = library.find((x) => x.key === key);
  const numbers = m.elements.map((el) => NUMBERS[el]), positions = Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM);
  const result = xtb.compute({ numbers, positions, charge: m.charge, etemp: 300, forces: false, unpaired: m.unpaired });
  assert.ok(result.ok, `${key}: ${result.error}`);
  return { m, numbers, positions, result, sites: chargeSites(numbers, positions, result) };
}
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };

console.log('charge sites against the multipoles of the engine');
{
  let worst = 0, worstTotal = 0, bonds = 0;
  const keys = ['h2o', 'nh3', 'c6h6', 'h2co', 'bh3nh3', 'o2', 'ch3', 'nacl', 'sf6', 'ticl4', 'ferrocene', 'nh4+', 'oh-', 'mno4-'];
  for (const key of keys) {
    const { m, numbers, positions, result, sites } = solve(key), nat = numbers.length;
    const atoms = chargeSites(numbers, positions, result, { bondShare: Infinity });
    assert.equal(atoms.n, nat);
    for (let a = 0; a < nat; a++) {
      worst = Math.max(worst, Math.abs(atoms.multipoles[10 * a] - result.charges[a]));
      for (let c = 0; c < 3; c++) worst = Math.max(worst, Math.abs(atoms.multipoles[10 * a + 1 + c] - result.dipoles[3 * a + c]));
      for (let c = 0; c < 6; c++) worst = Math.max(worst, Math.abs(atoms.multipoles[10 * a + 4 + c] - result.quadrupoles[6 * a + c]));
    }
    // charge and dipole moment of the whole: sum of q, sum of q r + d
    const total = (s) => {
      const t = [0, 0, 0, 0];
      for (let k = 0; k < s.n; k++) { t[0] += s.multipoles[10 * k]; for (let c = 0; c < 3; c++) t[1 + c] += s.multipoles[10 * k] * s.positions[3 * k + c] + s.multipoles[10 * k + 1 + c]; }
      return t;
    };
    const ta = total(atoms), ts = total(sites);
    worstTotal = Math.max(worstTotal, Math.abs(ts[0] - (m.charge || 0)), ...ts.map((v, i) => Math.abs(v - ta[i])));
    bonds += sites.n - nat;
  }
  ok(worst < 1e-9, `${keys.length} molecules, d metals, radicals and ions among them: sites on the atoms alone are the engine's charges, dipoles and quadrupoles (largest difference ${worst.toExponential(1)})`);
  ok(worstTotal < 1e-9 && bonds > 0, `with ${bonds} sites between the atoms the total charge and the dipole moment stay the same (largest difference ${worstTotal.toExponential(1)})`);
}

// The potential of the valence density, integrated on a grid of spacing h, plus the cores, at
// the given points (all outside the dense part of the cloud). Returns the potentials and the
// number of electrons the grid found.
function reference({ numbers, positions, result }, points, h) {
  const nat = numbers.length, { nao, nmo, shells, coeff, occ } = result, ao = new Float64Array(nao);
  const occupied = [];
  for (let o = 0; o < nmo; o++) if (occ[o] > 1e-9) occupied.push(o);
  let reach = 0;
  for (let a = 0; a < nat; a++) reach = Math.max(reach, Math.hypot(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]));
  const L = reach + 7, n = Math.round(2 * L / h), w = h ** 3, v = new Float64Array(points.length);
  let electrons = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) {
    const x = -L + (i + 0.5) * h, y = -L + (j + 0.5) * h, z = -L + (k + 0.5) * h;
    const d = density(x, y, z);
    if (d < 1e-10) continue;
    electrons += d * w;
    for (let q = 0; q < points.length; q++) { const p = points[q]; v[q] -= d * w / Math.hypot(p[0] - x, p[1] - y, p[2] - z); }
  }
  let expected = 0;
  for (const z of numbers) expected += VALENCE[z];
  for (let q = 0; q < points.length; q++) {
    v[q] *= expected / electrons; // what the grid lost or gained of the electrons would otherwise show as a net charge
    for (let a = 0; a < nat; a++) v[q] += VALENCE[numbers[a]] / Math.hypot(points[q][0] - positions[3 * a], points[q][1] - positions[3 * a + 1], points[q][2] - positions[3 * a + 2]);
  }
  return { v, electrons, expected };
  function density(x, y, z) {
    let at = 0;
    for (const sh of shells) at += shellValues(sh, x - positions[3 * sh.atom], y - positions[3 * sh.atom + 1], z - positions[3 * sh.atom + 2], ao, at);
    let t = 0;
    for (const o of occupied) { let psi = 0; for (let q = 0, b = o * nao; q < nao; q++) psi += coeff[b + q] * ao[q]; t += occ[o] * psi * psi; }
    return t;
  }
}
// points where the valence density is ESP_DENSITY: outwards from every nucleus in a few
// directions, kept if that nucleus is the nearest one (the others are inside the molecule)
function surfacePoints({ numbers, positions, result }, perAtom = 10) {
  const nat = numbers.length, ao = new Float64Array(result.nao), points = [];
  const density = (p) => {
    let at = 0;
    for (const sh of result.shells) at += shellValues(sh, p[0] - positions[3 * sh.atom], p[1] - positions[3 * sh.atom + 1], p[2] - positions[3 * sh.atom + 2], ao, at);
    let t = 0;
    for (let o = 0; o < result.nmo; o++) { if (result.occ[o] < 1e-9) continue; let psi = 0; for (let q = 0; q < result.nao; q++) psi += result.coeff[o * result.nao + q] * ao[q]; t += result.occ[o] * psi * psi; }
    return t;
  };
  for (let a = 0; a < nat; a++) for (let t = 0; t < perAtom; t++) {
    const z = 1 - (2 * t + 1) / perAtom, s = Math.sqrt(1 - z * z), ph = t * 2.399963 + a, dir = [s * Math.cos(ph), s * Math.sin(ph), z];
    const at = (d) => [0, 1, 2].map((c) => positions[3 * a + c] + d * dir[c]);
    let d = 9;
    while (d > 0.5 && density(at(d)) < ESP_DENSITY) d -= 0.05;
    const p = at(d + 0.025);
    let nearest = 0, best = Infinity;
    for (let b = 0; b < nat; b++) { const r = Math.hypot(p[0] - positions[3 * b], p[1] - positions[3 * b + 1], p[2] - positions[3 * b + 2]); if (r < best) { best = r; nearest = b; } }
    if (nearest === a) points.push(p);
  }
  return points;
}

console.log(`the potential on the surface of density ${ESP_DENSITY} against the potential of the density itself`);
const benzene = solve('c6h6');
// benzene: over the middle of the ring; beyond a hydrogen; in the plane of the ring between two hydrogens
function dist(p, a, b) { return Math.hypot(p[3 * a] - p[3 * b], p[3 * a + 1] - p[3 * b + 1], p[3 * a + 2] - p[3 * b + 2]); }
const special = (() => {
  const p = benzene.positions, at = (k) => [p[3 * k], p[3 * k + 1], p[3 * k + 2]];
  const a = at(0).map((v, i) => at(1)[i] - v), b = at(0).map((v, i) => at(2)[i] - v);
  const n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], l = Math.hypot(...n);
  const h1 = at(6), r = Math.hypot(...h1);
  // the hydrogen next to the first one round the ring
  let h2 = 7;
  for (let h = 8; h < 12; h++) if (dist(p, h, 6) < dist(p, h2, 6)) h2 = h;
  const mid = at(h2).map((v, i) => 0.5 * (v + h1[i])), ml = Math.hypot(...mid);
  return { face: n.map((v) => 3.4 * v / l), end: h1.map((v) => v * (r + 2.6) / r), between: mid.map((v) => v * (r + 1.2) / ml) };
})();
for (const [key, h, limit] of [['h2o', 0.08, 0.004], ['h2co', 0.09, 0.004], ['c6h6', 0.11, 0.004]]) {
  const mol = key === 'c6h6' ? benzene : solve(key), points = surfacePoints(mol);
  if (key === 'c6h6') points.push(special.face, special.end, special.between);
  const ref = reference(mol, points, h), atoms = chargeSites(mol.numbers, mol.positions, mol.result, { bondShare: Infinity });
  let sum = 0, sumAtoms = 0, lo = Infinity, hi = -Infinity;
  points.forEach((p, q) => {
    sum += (potential(mol.sites, p) - ref.v[q]) ** 2;
    // the engine's multipoles as bare points, without the size of the clouds
    sumAtoms += (potential({ ...atoms, pen: new Float32Array(atoms.pen.length) }, p) - ref.v[q]) ** 2;
    lo = Math.min(lo, ref.v[q]); hi = Math.max(hi, ref.v[q]);
  });
  const rms = Math.sqrt(sum / points.length), rmsAtoms = Math.sqrt(sumAtoms / points.length);
  ok(Math.abs(ref.electrons - ref.expected) < 0.02, `${key}: the grid finds ${ref.electrons.toFixed(3)} of ${ref.expected} electrons`);
  ok(rms < limit && rms < 0.6 * rmsAtoms, `${key}: ${points.length} points, potential from ${lo.toFixed(3)} to ${hi.toFixed(3)} Eh; the sites are off by ${rms.toFixed(4)} Eh (rms), the atomic multipoles alone by ${rmsAtoms.toFixed(4)}`);
  if (key === 'c6h6') {
    const [face, end, between] = ref.v.subarray(points.length - 3), model = [special.face, special.end, special.between].map((p) => potential(mol.sites, p));
    const bare = potential({ ...atoms, pen: new Float32Array(atoms.pen.length) }, special.between);
    ok(face < -0.005 && model[0] < -0.005, `benzene, over the ring: negative (density ${face.toFixed(4)}, sites ${model[0].toFixed(4)} Eh)`);
    ok(end > 0.005 && model[1] > 0.005, `benzene, beyond a hydrogen: positive (density ${end.toFixed(4)}, sites ${model[1].toFixed(4)} Eh)`);
    ok(between > 0 && model[2] > 0 && bare < 0, `benzene, between two hydrogens: positive (density ${between.toFixed(4)}, sites ${model[2].toFixed(4)} Eh); the atomic multipoles alone make it negative (${bare.toFixed(4)})`);
  }
}

console.log('signs');
{
  const along = (mol, from, to, d) => { // the point d bohr beyond atom `to` on the line from atom `from`
    const p = mol.positions, v = [0, 1, 2].map((c) => p[3 * to + c] - p[3 * from + c]), l = Math.hypot(...v);
    return v.map((x, c) => p[3 * to + c] + x * d / l);
  };
  const water = solve('h2o');
  const lone = potential(water.sites, [0, water.positions[1] + 3.2, 0]), acid = potential(water.sites, along(water, 0, 1, 2.6));
  ok(lone < -0.05 && acid > 0.03, `water: ${lone.toFixed(3)} Eh behind the oxygen, ${acid.toFixed(3)} beyond a hydrogen`);
  const bh3 = solve('bh3'), p = bh3.positions;
  const a = [0, 1, 2].map((c) => p[3 + c] - p[c]), b = [0, 1, 2].map((c) => p[6 + c] - p[c]);
  const n = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]], l = Math.hypot(...n);
  const hole = potential(bh3.sites, [0, 1, 2].map((c) => p[c] + 3 * n[c] / l)), hydride = potential(bh3.sites, along(bh3, 0, 1, 2.8));
  ok(hole > 0.02 && hole > hydride + 0.02, `borane: ${hole.toFixed(3)} Eh over the boron, where its empty orbital is; ${hydride.toFixed(3)} beyond a hydrogen`);
  const nh4 = solve('nh4+'), oh = solve('oh-');
  const cation = potential(nh4.sites, along(nh4, 0, 1, 2.6)), anion = potential(oh.sites, along(oh, 1, 0, 3.2));
  ok(cation > 0.2 && anion < -0.2, `ions: ${cation.toFixed(3)} Eh at the ammonium ion, ${anion.toFixed(3)} at the hydroxide ion`);
  const methane = solve('ch4');
  ok(potentialScale(methane.numbers, methane.positions, methane.sites) === 0.02 && potentialScale(water.numbers, water.positions, water.sites) > 0.06,
    `scale of the colours: ${potentialScale(water.numbers, water.positions, water.sites).toFixed(3)} Eh for water, the smallest allowed (0.02) for methane`);
}
console.log(`${checks} checks passed`);
