// Stage 4 acceptance, in node: BH3 and NH3 brought together by a mouse-like pull form the
// adduct.  node tests/adduct.test.mjs
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { Dynamics } from '../app/src/md.js';
import { buildScene } from '../app/src/scene.js';
import { fixSigns, frontier, overlapMatrix, followOrbitals } from '../app/src/orbitals.js';
import { BOHR_ANGSTROM, HARTREE_EV } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules
  .map((m) => ({ ...m, bohr: Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM) }));
const lib = (key) => library.find((m) => m.key === key);
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };
const dist = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]) * BOHR_ANGSTROM;

const scene = buildScene(lib('bh3'), lib('nh3'));
const B = 0, N = 4;
ok(scene.nFragments === 2 && scene.fragments.join('') === '00001111' && scene.elements.join('') === 'BHHHNHHH', 'the scene holds BH3 and NH3 as two fragments');
ok(Math.abs(dist(scene.positions, B, N) - 8 * BOHR_ANGSTROM) < 0.2, `they start ${dist(scene.positions, B, N).toFixed(2)} A apart (B to N)`);

const compute = (x, orbitals) => xtb.compute({ numbers: scene.numbers, positions: x, charge: 0, etemp: 300, forces: true, orbitals, accuracy: 0.1 });
const md = new Dynamics({ numbers: scene.numbers, positions: scene.positions, force: compute, seed: 5, wallRadius: 9 });
md.params.temperature = 300;
md.init(true);
md.thermalize(300);
let res = md.result;
fixSigns(res);
let selected = frontier(res).homo;
const weightOn = (r, m, atom) => {
  let w = 0, ao = 0;
  for (const sh of r.shells) { const n = sh.l ? 3 : 1; if (sh.atom === atom) for (let i = 0; i < n; i++) w += r.coeff[m * r.nao + ao + i] ** 2; ao += n; }
  return w;
};
const wN0 = weightOn(res, selected, N), wB0 = weightOn(res, selected, B);
ok(wN0 > 0.5 && wB0 < 0.005, `apart, the HOMO (${(res.eps[selected] * HARTREE_EV).toFixed(2)} eV) is the nitrogen lone pair: weight ${wN0.toFixed(2)} on N, ${wB0.toFixed(3)} on B`);
const e0 = md.epot;
// how far boron sits out of the plane of its three hydrogens
const pyramid = (x) => {
  const p = (k) => [x[3 * k], x[3 * k + 1], x[3 * k + 2]], sub = (a, b) => a.map((v, i) => v - b[i]);
  const u = sub(p(2), p(1)), v = sub(p(3), p(1)), n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const l = Math.hypot(...n), d = sub(p(B), p(1));
  return Math.abs(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]) / l * BOHR_ANGSTROM;
};
ok(pyramid(md.x) < 0.02, 'apart, borane is flat');

// one frame of the app: a step with orbitals, followed like the app follows them
const frameStep = () => {
  assert.ok(md.step(true), 'SCF failed');
  const next = md.result;
  selected = followOrbitals(res, next, overlapMatrix(next.shells, md.x), selected);
  res = next;
};
// the "mouse" holds nitrogen and carries it to 1.9 A from boron, then lets go
md.grab = { atom: N, target: [0, 0, 0] };
let carried = 0, onTheWay = null;
while (dist(md.x, B, N) > 2.0 && md.time < 3000) {
  const b = [md.x[0], md.x[1], md.x[2]], n = [md.x[3 * N], md.x[3 * N + 1], md.x[3 * N + 2]], l = Math.hypot(n[0] - b[0], n[1] - b[1], n[2] - b[2]);
  md.grab.target = b.map((v, i) => v + (n[i] - v) / l * 3.4); // a point 1.8 A from boron on the B-N line
  frameStep();
  carried = md.time;
  // the picture half-way in: is the followed orbital still the lone pair, and how much of it is on boron?
  if (process.env.TRACE && md.steps % 6 === 0) console.log(`    ${md.time.toFixed(0)} fs  B-N ${dist(md.x, B, N).toFixed(2)}  sel ${selected} homo ${frontier(res).homo}  eps ${[3, 4, 5, 6].map((m) => (res.eps[m] * HARTREE_EV).toFixed(2)).join(' ')}  sel wN ${weightOn(res, selected, N).toFixed(2)} wB ${weightOn(res, selected, B).toFixed(2)}  | per orb wN ${[4, 5, 6].map((m) => weightOn(res, m, N).toFixed(2)).join(' ')}`);
  if (!onTheWay && dist(md.x, B, N) < 2.6) onTheWay = { d: dist(md.x, B, N), wN: weightOn(res, selected, N), wB: weightOn(res, selected, B), homo: selected === frontier(res).homo };
}
ok(onTheWay && onTheWay.homo && onTheWay.wN > 0.5 && onTheWay.wB > 0.03,
  `on the way in (B-N ${onTheWay.d.toFixed(2)} A) the HOMO is still the lone pair and already reaches onto boron: weight ${onTheWay.wN.toFixed(2)} on N, ${onTheWay.wB.toFixed(2)} on B (was ${wB0.toFixed(3)})`);
ok(dist(md.x, B, N) <= 2.0, `the pull brings nitrogen within 2 A of boron in ${carried.toFixed(0)} fs without tearing anything (N-H max ${Math.max(...[5, 6, 7].map((h) => dist(md.x, N, h))).toFixed(2)} A)`);
md.grab = null;
let far = 0, esum = 0, psum = 0, n = 0;
const t0 = md.time;
while (md.time < t0 + 2000) {
  frameStep();
  if (md.time > t0 + 1000) { far = Math.max(far, dist(md.x, B, N)); esum += md.epot; psum += pyramid(md.x); n++; }
}
ok(far < 2.1, `released, the adduct holds: over the second picosecond B-N never exceeds ${far.toFixed(2)} A`);
ok(psum / n > 0.15, `borane bends out of its plane: boron sits ${(psum / n).toFixed(2)} A off the plane of its hydrogens on average`);
ok(esum / n - e0 < -0.02, `the potential energy fell by ${((e0 - esum / n) * HARTREE_EV).toFixed(2)} eV (${(e0 - esum / n).toFixed(4)} Eh) on average`);
// in the finished adduct the donated pair is the B-N bonding orbital; it has sunk below the
// two B-H orbitals, so it is no longer the HOMO
const { homo } = frontier(res);
let dative = -1;
for (let m = 0; m <= homo; m++) if (weightOn(res, m, N) > 0.3 && weightOn(res, m, B) > 0.05 && (dative < 0 || weightOn(res, m, B) > weightOn(res, dative, B))) dative = m;
ok(dative >= 0 && dative < homo, `in the adduct the pair shared by N and B is orbital ${dative} at ${(res.eps[dative] * HARTREE_EV).toFixed(2)} eV (weight ${weightOn(res, dative, N).toFixed(2)} on N, ${weightOn(res, dative, B).toFixed(2)} on B), below the HOMO at ${(res.eps[homo] * HARTREE_EV).toFixed(2)} eV`);
console.log(`  note: the orbital followed from the start ended as orbital ${selected} at ${(res.eps[selected] * HARTREE_EV).toFixed(2)} eV (weight ${weightOn(res, selected, N).toFixed(2)} on N, ${weightOn(res, selected, B).toFixed(2)} on B)`);
console.log(`${checks} checks passed`);

// ---- pieces: what hangs together now
{
  const { pieces } = await import('../app/src/scene.js');
  const count = (key, second) => pieces(buildScene(lib(key), second && lib(second))).count;
  ok(count('bh3nh3') === 1 && count('c6h6') === 1 && count('b2h6') === 1 && count('h2o') === 1, 'intact library molecules are one piece each (adduct, benzene, diborane, water)');
  ok(count('bh3', 'nh3') === 2 && count('pair') === 2, 'BH3 and NH3 side by side are two pieces, also when they come as one library entry');
  const parts = pieces({ elements: scene.elements, positions: md.x });
  ok(parts.count === 1, 'after the adduct has formed in the run above, the two former molecules are one piece');
}
