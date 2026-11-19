// Screens candidate "experiments": pairs of library molecules brought together the way the
// mouse does it (a spring on one nucleus), or simply left alone or heated, and reports what
// the scene consists of afterwards.   node tests/experiments-screen.mjs [filter]
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { buildScene, composition } from '../app/src/scene.js';
import { Dynamics } from '../app/src/md.js';
import { BOHR_ANGSTROM as A, HARTREE_EV } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const lib = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules.map((m) => ({ ...m, bohr: Float64Array.from(m.xyz, (v) => v / A) }));
const L = (k) => { const m = lib.find((x) => x.key === k); if (!m) throw new Error('no molecule ' + k); return m; };

// [first, second, pull: [element in first, element in second, stop distance A] | null, options]
const CASES = [
  ['bh3', 'nh3', ['B', 'N', 2.0]], ['bf3', 'nh3', ['B', 'N', 2.0]], ['bh3', 'bh3', ['B', 'B', 2.0]], ['bh3', 'co', ['B', 'C', 1.9]],
  ['bh3', 'h2o', ['B', 'O', 1.9]], ['c2h4', 'bh3', ['C', 'B', 1.9]], ['bh3', 'f-', ['B', 'F', 1.8]],
  ['h3o+', 'nh3', null], ['h3o+', 'oh-', null], ['h3o+', 'oh-', null, { etemp: 1500 }], ['nh4+', 'oh-', null, { etemp: 1500 }], ['h3o+', 'f-', null, { etemp: 1500 }],
  ['hf', 'oh-', ['H', 'O', 1.3]], ['hf', 'nh3', ['H', 'N', 1.3]], ['h3o+', 'cn-', null, { etemp: 1500 }], ['h3o+', 'c2h4', ['O', 'C', 2.4]], ['h3o+', 'h2o', ['O', 'O', 2.6]],
  ['ch3cooh', 'oh-', ['O', 'O', 2.6]], ['hcooh', 'oh-', ['O', 'O', 2.6]], ['ch3cooh', 'nh3', ['O', 'N', 2.7]], ['hno3', 'nh3', ['O', 'N', 2.7]], ['hno3', 'oh-', ['N', 'O', 2.8]],
  ['lih', 'h2o', ['H', 'H', 1.0]], ['nah', 'h2o', ['H', 'H', 1.0]], ['lih', 'hf', ['H', 'H', 1.0]], ['nah', 'hf', ['H', 'H', 1.0]], ['bh4-', 'h3o+', null, { etemp: 1500 }], ['lih', 'co2', ['H', 'C', 1.4]],
  ['naoh', 'hf', ['O', 'H', 1.3]], ['lioh', 'hf', ['O', 'H', 1.3]], ['naoh', 'co2', ['O', 'C', 1.7]], ['oh-', 'co2', ['O', 'C', 1.7]], ['h2o', 'co2', ['O', 'C', 1.7]], ['nh3', 'co2', ['N', 'C', 1.8]],
  ['naoh', 'hcooh', ['O', 'O', 2.6]], ['naoh', 'h3o+', null, { etemp: 1500 }],
  ['lif', 'lif', ['Li', 'F', 1.9]], ['naf', 'naf', ['Na', 'F', 2.1]], ['lih', 'lih', ['Li', 'H', 1.8]], ['h2o', 'h2o', ['O', 'O', 2.9]], ['hf', 'hf', ['F', 'F', 2.7]], ['naf', 'h2o', ['Na', 'O', 2.3]],
  ['ch3f', 'oh-', ['C', 'O', 1.9]], ['ch3f', 'f-', ['C', 'F', 1.9]], ['ch3f', 'cn-', ['C', 'C', 1.9]],
  ['li2', 'f2', null, { etemp: 1500 }], ['h2', 'f2', ['H', 'F', 1.1], { etemp: 1500 }], ['beh2', 'beh2', ['Be', 'H', 1.6]], ['hcn', 'h3o+', ['N', 'O', 2.6]],
  // oxygen
  ['h2', 'o2', null, { T: 3000, etemp: 3000, ps: 5 }], ['h2', 'o2', ['H', 'O', 1.0], { etemp: 1500 }], ['li2', 'o2', null, { etemp: 1500 }], ['li2', 'o2', null],
  ['bh3', 'o2', ['B', 'O', 1.5], { etemp: 1500 }], ['lih', 'o2', ['Li', 'O', 1.8], { etemp: 1500 }], ['co', 'o2', ['C', 'O', 1.4], { etemp: 1500 }],
  ['ch4', 'o2', null, { T: 4000, etemp: 3000, ps: 5 }], ['beh2', 'o2', ['Be', 'O', 1.5], { etemp: 1500 }], ['nh3', 'o2', null, { T: 3000, etemp: 3000, ps: 5 }],
  // radicals and spin (tag r). Unpaired electrons: the scene's own rule unless given.
  ['h', 'h', ['H', 'H', 1.0], { tag: 'r' }], ['h', 'h', ['H', 'H', 1.0], { tag: 'r', unpaired: 2 }], ['h', 'h', null, { tag: 'r' }],
  ['oh', 'h', ['O', 'H', 1.1], { tag: 'r' }], ['ch3', 'h', ['C', 'H', 1.2], { tag: 'r' }], ['ch3', 'ch3', ['C', 'C', 1.7], { tag: 'r' }], ['nh2', 'h', ['N', 'H', 1.1], { tag: 'r' }],
  ['oh', 'oh', ['O', 'O', 1.6], { tag: 'r' }], ['li', 'f', ['Li', 'F', 1.7], { tag: 'r' }], ['na', 'f', ['Na', 'F', 2.0], { tag: 'r' }], ['h', 'f', ['H', 'F', 1.0], { tag: 'r' }],
  ['li', 'f', null, { tag: 'r' }], ['na', 'f', null, { tag: 'r' }],
  ['h2', 'f', ['H', 'F', 1.1], { tag: 'r' }], ['f2', 'h', ['F', 'H', 1.1], { tag: 'r' }], ['o2', 'h', ['O', 'H', 1.1], { tag: 'r' }], ['h2', 'oh', ['H', 'O', 1.1], { tag: 'r' }],
  ['ch4', 'f', ['H', 'F', 1.1], { tag: 'r' }], ['ch4', 'oh', ['H', 'O', 1.2], { tag: 'r' }], ['h2o', 'f', ['H', 'F', 1.1], { tag: 'r' }], ['h2', 'o', ['H', 'O', 1.1], { tag: 'r' }],
  ['o3', 'no', ['O', 'N', 1.5], { tag: 'r' }], ['o2', 'o', ['O', 'O', 1.4], { tag: 'r' }], ['o', 'o', ['O', 'O', 1.3], { tag: 'r' }], ['o', 'o', ['O', 'O', 1.3], { tag: 'r', unpaired: 2 }], ['no', 'o', ['N', 'O', 1.3], { tag: 'r' }],
  ['h2o', 'na', ['O', 'Na', 2.2], { tag: 'r' }], ['h2o', 'li', ['O', 'Li', 1.9], { tag: 'r' }], ['hf', 'li', ['F', 'Li', 1.8], { tag: 'r' }], ['hf', 'na', ['F', 'Na', 2.0], { tag: 'r' }],
  // oxygen as a triplet, default electronic temperature
  ['h2', 'o2', null, { tag: 'r', ps: 4 }], ['h2', 'o2', ['H', 'O', 1.0], { tag: 'r' }], ['li2', 'o2', null, { tag: 'r' }], ['bh3', 'o2', ['B', 'O', 1.5], { tag: 'r' }], ['co', 'o2', ['C', 'O', 1.4], { tag: 'r' }],
  ['ch4', 'o2', ['H', 'O', 1.1], { tag: 'r' }], ['li2', 'o2', null, { tag: 'r', unpaired: 0 }],
  // one molecule, heated
  ['h2co', null, null, { T: 4000, ps: 5 }], ['o3', null, null, { T: 3000, ps: 5, etemp: 3000 }], ['h2o2', null, null, { T: 3000, ps: 5, etemp: 3000 }], ['h2co3', null, null, { T: 2500, ps: 5 }],
  ['hcooh', null, null, { T: 3500, ps: 5 }], ['n2h4', null, null, { T: 3000, ps: 5, etemp: 3000 }], ['ch3oh', null, null, { T: 4000, ps: 5 }], ['b2h6', null, null, { T: 2000, ps: 5 }],
];

const filter = process.argv[2];
for (const [a, b, pull, opt = {}] of CASES) {
  if (filter && !(filter.startsWith('#') ? opt.tag === filter.slice(1) : `${a}+${b}`.includes(filter))) continue;
  const { T = 300, etemp = 300, ps = 2 } = opt;
  const scene = buildScene(L(a), b && L(b)), n1 = L(a).elements.length, unpaired = opt.unpaired ?? scene.unpaired;
  const ia = pull ? scene.elements.indexOf(pull[0]) : -1, ib = pull ? n1 + scene.elements.slice(n1).indexOf(pull[1]) : -1;
  const out = [];
  for (let seed = 1; seed <= 2; seed++) {
    let wall = 0; for (let k = 0; k < scene.positions.length; k += 3) wall = Math.max(wall, Math.hypot(...scene.positions.subarray(k, k + 3)));
    const md = new Dynamics({ numbers: scene.numbers, positions: scene.positions, seed, wallRadius: wall + 2.5,
      force: (x) => xtb.compute({ numbers: scene.numbers, positions: x, charge: scene.charge, etemp, unpaired, forces: true, orbitals: false, accuracy: 0.1 }) });
    md.params.temperature = T;
    if (!md.init(false)) { out.push('SCF failed at start'); continue; }
    const e0 = md.epot;
    md.thermalize(T);
    const d = () => Math.hypot(...[0, 1, 2].map((k) => md.x[3 * ia + k] - md.x[3 * ib + k])) * A;
    let failed = null, carried = null;
    if (pull) {
      md.grab = { atom: ib, target: [0, 0, 0] };
      while (md.time < 1500 && d() > pull[2]) {
        const pa = [0, 1, 2].map((k) => md.x[3 * ia + k]), pb = [0, 1, 2].map((k) => md.x[3 * ib + k]), len = Math.hypot(...pb.map((v, k) => v - pa[k]));
        md.grab.target = pa.map((v, k) => v + (pb[k] - v) / len * (pull[2] - 0.1) / A);
        if (!md.step()) { failed = md.error.error; break; }
      }
      carried = md.time;
      md.grab = null;
    }
    const t0 = md.time;
    let es = 0, n = 0;
    while (!failed && md.time < t0 + ps * 1000) { if (!md.step()) { failed = md.error.error; break; } if (md.time > t0 + ps * 1000 - 500) { es += md.epot; n++; } }
    out.push(`${composition({ elements: scene.elements, positions: md.x })}${n ? ` (${((es / n - e0) * HARTREE_EV).toFixed(2)} eV)` : ''}${pull ? ` [pull ${carried.toFixed(0)} fs, ${pull[0]}-${pull[1]} now ${d().toFixed(2)}]` : ''}${failed ? ` SCF FAILED at ${md.time.toFixed(0)} fs` : ''}`);
  }
  console.log(`${`${a}${b ? ' + ' + b : ''}`.padEnd(18)} ${`T ${T}${etemp !== 300 ? ' Te ' + etemp : ''} u${unpaired}`.padEnd(18)} ${out.join('  |  ')}`);
}
