// The molecule editor's model and the relaxation of what it builds:  node tests/editor.test.mjs
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { minimize, modelGuess, modelHessian } from '../app/src/relax.js';
import * as E from '../app/src/editor.js';
import { composition } from '../app/src/scene.js';
import { BOHR_ANGSTROM, NUMBERS } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };

const at = (m, k) => m.positions.slice(3 * k, 3 * k + 3);
const dist = (m, i, j) => Math.hypot(...at(m, i).map((v, c) => v - m.positions[3 * j + c])) * BOHR_ANGSTROM;
const angle = (m, i, k, j) => {
  const u = at(m, i).map((v, c) => v - m.positions[3 * k + c]), w = at(m, j).map((v, c) => v - m.positions[3 * k + c]);
  return Math.acos(u.reduce((s, v, c) => s + v * w[c], 0) / Math.hypot(...u) / Math.hypot(...w)) * 180 / Math.PI;
};
const formula = (m) => composition({ elements: m.elements, positions: Float64Array.from(m.positions) });
const heavy = (m) => m.elements.map((el, k) => (el === 'H' ? -1 : k)).filter((k) => k >= 0);
const firstH = (m, k) => E.hydrogens(m, k)[0];
const closest = (m) => { let d = 9; for (let i = 0; i < m.elements.length; i++) for (let j = i + 1; j < m.elements.length; j++) d = Math.min(d, dist(m, i, j)); return d; };
// every angle between the bonds of atom k
const angles = (m, k) => { const nb = E.neighbours(m, k).map((n) => n.atom), out = []; for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++) out.push(angle(m, nb[a], k, nb[b])); return out; };
const within = (list, want, tol) => list.every((v) => Math.abs(v - want) < tol);
const bond = (m, i, j) => m.bonds.findIndex(([a, b]) => (a === i && b === j) || (a === j && b === i));

function relax(m, unpaired) {
  const numbers = m.elements.map((el) => NUMBERS[el]);
  let calls = 0;
  const walk = minimize((x) => {
    calls++;
    const r = xtb.compute({ numbers, positions: x, charge: m.charge, etemp: 300, forces: true, orbitals: false, unpaired });
    return r.ok ? { e: r.energy, g: r.forces.map((f) => -f) } : { e: Infinity };
  }, m.positions, { gtol: 5e-5, maxIter: 600, guess: (x) => modelGuess(numbers, x) });
  for (;;) { const s = walk.next(); if (s.done) return { ...s.value, calls, mol: { ...m, positions: Array.from(s.value.x) } }; }
}
const libraryEnergy = (key) => library.find((x) => x.key === key).energy;

console.log('building by hydrogens');
let m = E.placeAtom(E.emptyMolecule(), 'C', [0, 0, 0]);
ok(formula(m) === 'CH4' && within(angles(m, 0), 109.47, 0.5) && Math.abs(dist(m, 0, 1) - 1.07) < 0.01, `a carbon arrives as CH4: tetrahedral, C-H ${dist(m, 0, 1).toFixed(2)} A`);
const methane = m;
m = E.replaceHydrogen(m, firstH(m, 0), 'C');
let [c1, c2] = heavy(m);
ok(formula(m) === 'C2H6' && Math.abs(dist(m, c1, c2) - 1.52) < 0.01, `clicking a hydrogen with carbon gives ethane, C-C ${dist(m, c1, c2).toFixed(2)} A`);
ok(within(angles(m, c1), 109.47, 3) && within(angles(m, c2), 109.47, 3), 'both carbons stay tetrahedral');
{
  // staggered: every H-C-C-H dihedral is 60 or 180 degrees
  const axis = at(m, c2).map((v, c) => v - m.positions[3 * c1 + c]), l = Math.hypot(...axis), a = axis.map((v) => v / l);
  const perp = (k, c0) => { const v = at(m, k).map((x, c) => x - m.positions[3 * c0 + c]), s = v.reduce((t, x, c) => t + x * a[c], 0); const p = v.map((x, c) => x - s * a[c]), n = Math.hypot(...p); return p.map((x) => x / n); };
  const dihedrals = E.hydrogens(m, c1).flatMap((h1) => E.hydrogens(m, c2).map((h2) => Math.acos(Math.max(-1, Math.min(1, perp(h1, c1).reduce((t, x, c) => t + x * perp(h2, c2)[c], 0)))) * 180 / Math.PI));
  ok(dihedrals.every((v) => Math.abs(v - 60) < 12 || Math.abs(v - 180) < 12), `the new hydrogens are staggered against the old ones (dihedrals ${dihedrals.map((v) => v.toFixed(0)).join(' ')})`);
}
const ethane = m;
m = E.cycleBond(m, bond(m, c1, c2));
[c1, c2] = heavy(m);
ok(formula(m) === 'C2H4' && m.bonds[bond(m, c1, c2)][2] === 2, 'clicking the bond makes it double and takes a hydrogen from each end: ethene');
{
  const n = E.hydrogens(m, c1).map((h) => at(m, h).map((v, c) => v - m.positions[3 * c1 + c]));
  const normal = [n[0][1] * n[1][2] - n[0][2] * n[1][1], n[0][2] * n[1][0] - n[0][0] * n[1][2], n[0][0] * n[1][1] - n[0][1] * n[1][0]];
  const out = [...E.hydrogens(m, c2), c2].map((k) => Math.abs(at(m, k).reduce((s, v, c) => s + (v - m.positions[3 * c1 + c]) * normal[c], 0)) / Math.hypot(...normal) * BOHR_ANGSTROM);
  ok(Math.max(...out) < 0.08 && within(angles(m, c1), 120, 4), `ethene comes out flat (atoms within ${Math.max(...out).toFixed(3)} A of one plane), angles near 120`);
}
const ethene = m;
m = E.cycleBond(m, bond(m, c1, c2));
[c1, c2] = heavy(m);
ok(formula(m) === 'C2H2' && within(angles(m, c1), 180, 3) && within(angles(m, c2), 180, 3), 'once more: ethyne, linear');
m = E.cycleBond(m, bond(m, c1, c2));
ok(formula(m) === 'C2H6', 'and once more: back to a single bond, the hydrogens return');
assert.throws(() => E.cycleBond(E.replaceHydrogen(E.placeAtom(E.emptyMolecule(), 'O', [0, 0, 0]), 1, 'F'), 0), /hydrogen/);
ok(true, 'a bond whose end has no hydrogen to give (HO-F) refuses to become double, with a message');

{
  // hydrogens erased by hand leave free valences, and those make bonds too
  let e = ethane;
  for (const c of [0, 1]) e = E.erase(e, firstH(e, heavy(e)[c]));
  e = E.cycleBond(e, bond(e, ...heavy(e)));
  ok(formula(e) === 'C2H4' && e.bonds[bond(e, ...heavy(e))][2] === 2, 'ethane with a hydrogen erased from each carbon: the bond becomes double without taking more hydrogens');
  let bare = ethane;
  while (bare.elements.includes('H')) bare = E.erase(bare, bare.elements.indexOf('H'));
  bare = E.cycleBond(E.cycleBond(bare, 0), 0);
  ok(bare.elements.join('') === 'CC' && bare.bonds[0][2] === 3, 'two bare carbons can be joined by a double and a triple bond');
  const nh4 = E.addHydrogen(E.setElement(methane, 0, 'N'), 0);
  ok(formula(nh4) === 'H4N' && within(angles(nh4, 0), 109.47, 1), 'one more hydrogen on ammonia: four, tetrahedral (NH4+ once the charge is set)');
}
console.log('changing and erasing');
m = E.setElement(methane, 0, 'O');
ok(formula(m) === 'H2O' && Math.abs(angles(m, 0)[0] - 109.5) < 3, 'carbon of methane changed to oxygen: water, bent');
m = E.setElement(methane, 0, 'N');
ok(formula(m) === 'H3N' && within(angles(m, 0), 109.5, 3), 'changed to nitrogen: ammonia, a pyramid');
m = E.setElement(ethane, heavy(ethane)[1], 'O');
ok(formula(m) === 'CH4O', 'a carbon of ethane changed to oxygen: methanol');
m = E.erase(m, heavy(m)[1]);
ok(formula(m) === 'CH4', 'erasing the oxygen leaves methane: the carbon gets its hydrogen back');
m = E.erase(methane, 1);
ok(formula(m) === 'CH3' && within(angles(m, 0), 109.47, 0.5), 'erasing a hydrogen leaves the methyl radical, the rest untouched');
m = E.replaceHydrogen(E.setElement(ethene, heavy(ethene)[1], 'O'), 2, 'C');
m = E.erase(m, heavy(m).find((k) => m.elements[k] === 'O'));
ok(formula(m) === 'C2H6', 'erasing the oxygen of acetaldehyde returns two hydrogens for its double bond: ethane');
assert.throws(() => E.erase(E.placeAtom(E.emptyMolecule(), 'Fe', [0, 0, 0]), 0), /last atom/);
ok(E.placeAtom(E.emptyMolecule(), 'Fe', [0, 0, 0]).elements.length === 1, 'an element without a usual number of bonds arrives bare; the last atom cannot be erased');
m = E.removeBond(ethane, bond(ethane, ...heavy(ethane)));
ok(m.elements.length === 10 && m.bonds.length === 8, 'removing the C-C bond of ethane gives two methanes');
m = E.moveAtom(ethane, heavy(ethane)[0], [1, 0, 0]);
ok(E.hydrogens(ethane, heavy(ethane)[0]).every((h) => Math.abs(m.positions[3 * h] - ethane.positions[3 * h] - 1) < 1e-12) && m.positions[3 * heavy(ethane)[1]] === ethane.positions[3 * heavy(ethane)[1]], 'a moved atom takes its hydrogens along and nothing else');

console.log('closing a ring');
let chain = methane;
for (let i = 1; i < 6; i++) { const last = heavy(chain).at(-1); chain = E.replaceHydrogen(chain, E.hydrogens(chain, last).at(-1), 'C'); }
ok(formula(chain) === 'C6H14' && closest(chain) > 0.95, `six clicks give hexane; no two atoms closer than ${closest(chain).toFixed(2)} A`);
{
  const cs = heavy(chain), far = dist(chain, cs[0], cs[5]);
  const ring = E.connect(chain, cs[0], cs[5]), rc = heavy(ring);
  const cc = ring.bonds.filter(([i, j]) => ring.elements[i] === 'C' && ring.elements[j] === 'C').map(([i, j]) => dist(ring, i, j));
  ok(formula(ring) === 'C6H12' && Math.max(...cc) < 1.75 && Math.min(...cc) > 1.3, `joining the ends (${far.toFixed(1)} A apart) gives cyclohexane with C-C between ${Math.min(...cc).toFixed(2)} and ${Math.max(...cc).toFixed(2)} A`);
  ok(closest(ring) > 0.9 && rc.every((k) => E.hydrogens(ring, k).length === 2), `every carbon has two hydrogens, nothing collides (closest pair ${closest(ring).toFixed(2)} A)`);
  const r = relax(ring);
  const cc2 = ring.bonds.filter(([i, j]) => ring.elements[i] === 'C' && ring.elements[j] === 'C').map(([i, j]) => dist(r.mol, i, j));
  ok(r.converged && Math.max(...cc2) - Math.min(...cc2) < 0.02 && formula(r.mol) === 'C6H12', `the engine relaxes it in ${r.it} steps to a ring with all C-C ${cc2[0].toFixed(3)} A`);
}

console.log('relaxing what the editor builds against the library');
for (const [name, key, mol] of [['methane', 'ch4', methane], ['ethane', 'c2h6', ethane], ['ethene', 'c2h4', ethene],
  ['water', 'h2o', E.setElement(methane, 0, 'O')], ['ammonia', 'nh3', E.setElement(methane, 0, 'N')],
  ['methanol', 'ch3oh', E.setElement(ethane, heavy(ethane)[1], 'O')]]) {
  const r = relax(mol), de = r.e - libraryEnergy(key);
  ok(r.converged && Math.abs(de) < 2e-5, `${name}: ${r.it} steps, ${r.calls} calculations, energy ${de >= 0 ? '+' : ''}${de.toExponential(1)} Eh from the library's`);
}

console.log('geometries from elsewhere');
for (const key of ['c2h4', 'c2h2', 'co2', 'h2co', 'hcn', 'c6h6', 'ch3cooh', 'n2']) {
  const lib = library.find((x) => x.key === key), mol = E.fromXYZ(E.toXYZ(E.withBonds(lib.elements, lib.xyz.map((v) => v / BOHR_ANGSTROM))));
  const orders = mol.bonds.filter(([i, j]) => mol.elements[i] !== 'H' && mol.elements[j] !== 'H').map((b) => b[2]).sort().join('');
  const want = { c2h4: '2', c2h2: '3', co2: '22', h2co: '2', hcn: '3', c6h6: '111111', ch3cooh: '112', n2: '3' }[key];
  ok(orders === want && mol.elements.length === lib.elements.length, `${key}: through XYZ and back, bond orders read from the distances: ${orders}`);
}
{
  const lib = library.find((x) => x.key === 'h2co'), mol = E.withBonds(lib.elements, lib.xyz.map((v) => v / BOHR_ANGSTROM));
  const m2 = E.erase(mol, mol.elements.indexOf('O'));
  ok(formula(m2) === 'CH4', 'library formaldehyde, oxygen erased: methane');
  assert.throws(() => E.fromXYZ('Xx 0 0 0'), /not among/);
  assert.throws(() => E.fromXYZ('nothing here'), /No atoms/);
  ok(E.fromXYZ('3\nwater\nO 0 0 0\nH 0.96 0 0\nh -0.24 0.93 0').bonds.length === 2, 'XYZ with a header, any letter case; text without atoms and unknown elements are refused');
}

console.log('the model Hessian the relaxation starts from');
{
  // hexane as the editor built it; acetylene for the straight angles the model leaves out
  for (const [name, m] of [['hexane', chain], ['acetylene', E.fromXYZ('C 0 0 0.6\nC 0 0 -0.6\nH 0 0 1.66\nH 0 0 -1.66')]]) {
    const numbers = m.elements.map((el) => NUMBERS[el]), x = m.positions, nat = numbers.length, n = 3 * nat, H = modelHessian(numbers, x);
    const times = (v) => Float64Array.from({ length: n }, (_, i) => { let t = 0; for (let j = 0; j < n; j++) t += H[i * n + j] * v[j]; return t; });
    const largest = (v) => v.reduce((t, w) => Math.max(t, Math.abs(w)), 0);
    // the molecule moved or turned as a whole: three shifts and three turns
    let rigid = 0;
    for (let p = 0; p < 3; p++) {
      const q = (p + 1) % 3, r = (p + 2) % 3, shift = new Float64Array(n), turn = new Float64Array(n);
      for (let a = 0; a < nat; a++) { shift[3 * a + p] = 1; turn[3 * a + q] = -x[3 * a + r]; turn[3 * a + r] = x[3 * a + q]; }
      rigid = Math.max(rigid, largest(times(shift)), largest(times(turn)));
    }
    // one atom moved along its bond: the stiffness of a bond
    const [i, j] = m.bonds[0], l = dist(m, i, j) / BOHR_ANGSTROM, pull = new Float64Array(n);
    for (let p = 0; p < 3; p++) pull[3 * i + p] = (x[3 * i + p] - x[3 * j + p]) / l;
    const stiff = times(pull).reduce((t, w, k) => t + w * pull[k], 0);
    const G = modelGuess(numbers, x), back = times(G[0]).map((w, k) => w + 0.01 * G[0][k] - (k === 0 ? 1 : 0));
    ok(rigid < 1e-12 && stiff > 0.3 && stiff < 2 && largest(back) < 1e-10,
      `${name}: moving or turning it as a whole costs nothing (${rigid.toExponential(0)}), a bond is ${stiff.toFixed(2)} Eh/bohr² stiff, the guess is the inverse`);
  }
}
console.log(`\n${checks} checks passed`);
