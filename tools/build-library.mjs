// Builds app/src/library.json: every molecule of the library with its geometry optimized
// by the GFN2-xTB engine itself.
//   node tools/build-library.mjs
import { readFileSync, writeFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { BOHR_ANGSTROM, NUMBERS } from '../app/src/units.js';
import { pieces } from '../app/src/scene.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const R = Math.PI / 180, cos = (a) => Math.cos(a * R), sin = (a) => Math.sin(a * R);

// starting geometries in angstrom: [element, x, y, z]
// the geometry this molecule already has in the library (the first few came from the reference
// demo, which is gone; 'pair' is set by hand and lives only there)
const current = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
const kept = (key) => { const m = current.find((x) => x.key === key); return m.elements.map((el, i) => [el, ...m.xyz.slice(3 * i, 3 * i + 3)]); };
const ring = (n, r, el, phase = 0) => Array.from({ length: n }, (_, i) => [typeof el === 'function' ? el(i) : el, r * cos(phase + 360 * i / n), r * sin(phase + 360 * i / n), 0]);
const T = 1.09 / Math.sqrt(3);
const tetra = (centre, r, h = 'H') => { const t = r / Math.sqrt(3); return [[centre, 0, 0, 0], [h, t, t, t], [h, t, -t, -t], [h, -t, t, -t], [h, -t, -t, t]]; };

// Internal coordinates to Cartesian. Rows: [el], [el, a, r], [el, a, r, b, angle],
// [el, a, r, b, angle, c, dihedral]; a, b, c are 1-based numbers of earlier atoms.
function zmat(rows) {
  const out = [];
  const sub = (u, v) => u.map((x, i) => x - v[i]), cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  const unit = (u) => { const l = Math.hypot(...u); return u.map((x) => x / l); };
  rows.forEach(([el, a, r, b, ang, c, dih], i) => {
    if (i === 0) return out.push([el, 0, 0, 0]);
    const A = out[a - 1].slice(1);
    if (i === 1) return out.push([el, A[0] + r, A[1], A[2]]);
    const B = out[b - 1].slice(1), bc = unit(sub(A, B));
    // reference direction for the dihedral: atom c, or any direction off the a-b line for the third atom
    const ref = c ? sub(B, out[c - 1].slice(1)) : (Math.abs(bc[2]) < 0.9 ? [0, 0, 1] : [0, 1, 0]);
    const n = unit(cross(ref, bc)), m = cross(n, bc), t = (dih ?? 0) * R, th = ang * R;
    const d = [-r * Math.cos(th), r * Math.sin(th) * Math.cos(t), r * Math.sin(th) * Math.sin(t)];
    out.push([el, ...A.map((x, k) => x + d[0] * bc[k] + d[1] * m[k] + d[2] * n[k])]);
  });
  return out;
}
const octa = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]];
const methyl = (c, x, ref, first = 180) => [first, first + 120, first - 120].map((d) => ['H', c, 1.09, x, 110, ref, d]);

// group: heading in the list of molecules
const molecules = [
  { key: 'h2o', name: 'H₂O, water', group: 'Simple', start: kept('h2o') },
  { key: 'nh3', name: 'NH₃, ammonia', group: 'Simple', start: kept('nh3') },
  { key: 'ch4', name: 'CH₄, methane', group: 'Simple', start: tetra('C', 1.09) },
  { key: 'h2', name: 'H₂, hydrogen', group: 'Simple', start: zmat([['H'], ['H', 1, 0.74]]) },
  { key: 'n2', name: 'N₂, nitrogen', group: 'Simple', start: [['N', 0.55, 0, 0], ['N', -0.55, 0, 0]] },
  { key: 'o2', name: 'O₂, oxygen', group: 'Simple', unpaired: 2, start: zmat([['O'], ['O', 1, 1.21]]),
    note: 'Oxygen has two unpaired electrons: one in each of the two highest occupied orbitals, with parallel spins (the ↑○ mark). The number of unpaired electrons can be changed under “Settings” in the dynamics panel.' },
  { key: 'co', name: 'CO, carbon monoxide', group: 'Simple', start: zmat([['C'], ['O', 1, 1.13]]) },
  { key: 'co2', name: 'CO₂, carbon dioxide', group: 'Simple', start: [['C', 0, 0, 0], ['O', 1.16, 0, 0], ['O', -1.16, 0, 0]] },
  { key: 'o3', name: 'O₃, ozone', group: 'Simple', start: zmat([['O'], ['O', 1, 1.27], ['O', 1, 1.27, 2, 117]]) },
  { key: 'n2o', name: 'N₂O, nitrous oxide', group: 'Simple', start: [['N', -1.13, 0, 0], ['N', 0, 0, 0], ['O', 1.19, 0, 0]] },
  { key: 'h2o2', name: 'H₂O₂, hydrogen peroxide', group: 'Simple', start: zmat([['O'], ['O', 1, 1.45], ['H', 1, 0.97, 2, 100], ['H', 2, 0.97, 1, 100, 3, 115]]) },
  { key: 'n2h4', name: 'N₂H₄, hydrazine', group: 'Simple', start: zmat([['N'], ['N', 1, 1.45], ['H', 1, 1.02, 2, 108], ['H', 1, 1.02, 2, 108, 3, 115], ['H', 2, 1.02, 1, 108, 3, 90], ['H', 2, 1.02, 1, 108, 5, 115]]) },
  { key: 'hno3', name: 'HNO₃, nitric acid', group: 'Simple', start: zmat([['N'], ['O', 1, 1.21], ['O', 1, 1.21, 2, 130], ['O', 1, 1.41, 2, 115, 3, 180], ['H', 4, 0.97, 1, 102, 2, 0]]) },

  { key: 'bh3', name: 'BH₃, borane', group: 'Boron', start: kept('bh3') },
  { key: 'bh3nh3', name: 'BH₃NH₃, adduct', group: 'Boron', start: kept('bh3nh3') },
  { key: 'pair', name: 'BH₃ + NH₃ apart, 3 Å', group: 'Boron', start: kept('pair'), fixed: true,
    note: 'The geometry is set by hand: two fragments at a B–N distance of 3 Å, not optimized.' },
  { key: 'b2h6', name: 'B₂H₆, diborane', group: 'Boron', start: [['B', 0.885, 0, 0], ['B', -0.885, 0, 0], ['H', 0, 0, 0.99], ['H', 0, 0, -0.99],
    ...[[1, 1], [1, -1], [-1, 1], [-1, -1]].map(([a, b]) => ['H', a * (0.885 + 1.19 * cos(60.5)), b * 1.19 * sin(60.5), 0])] },
  { key: 'h2bnh2', name: 'H₂BNH₂, aminoborane', group: 'Boron', start: [['B', -0.70, 0, 0], ['N', 0.69, 0, 0],
    ...[1, -1].map((q) => ['H', -0.70 - 1.19 * cos(60), q * 1.19 * sin(60), 0]), ...[1, -1].map((q) => ['H', 0.69 + 1.01 * cos(57), q * 1.01 * sin(57), 0])] },
  { key: 'boh3', name: 'B(OH)₃, boric acid', group: 'Boron', start: zmat([['B'], ['O', 1, 1.37], ['O', 1, 1.37, 2, 120], ['O', 1, 1.37, 2, 120, 3, 180],
    ['H', 2, 0.96, 1, 113, 3, 180], ['H', 3, 0.96, 1, 113, 4, 180], ['H', 4, 0.96, 1, 113, 2, 180]]) },
  { key: 'b3n3h6', name: 'B₃N₃H₆, borazine', group: 'Boron', start: [...ring(6, 1.43, (i) => 'BN'[i % 2]), ...ring(6, 2.55, 'H')] },

  { key: 'c2h4', name: 'C₂H₄, ethylene', group: 'Organic', start: kept('c2h4') },
  { key: 'c2h2', name: 'C₂H₂, acetylene', group: 'Organic', start: [['H', -1.66, 0, 0], ['C', -0.60, 0, 0], ['C', 0.60, 0, 0], ['H', 1.66, 0, 0]] },
  { key: 'c2h6', name: 'C₂H₆, ethane', group: 'Organic', start: zmat([['C'], ['C', 1, 1.53], ['H', 1, 1.09, 2, 110.5], ['H', 1, 1.09, 2, 110.5, 3, 120], ['H', 1, 1.09, 2, 110.5, 3, -120],
    ['H', 2, 1.09, 1, 110.5, 3, 180], ['H', 2, 1.09, 1, 110.5, 3, 60], ['H', 2, 1.09, 1, 110.5, 3, -60]]) },
  { key: 'hcn', name: 'HCN, hydrogen cyanide', group: 'Organic', start: [['H', -1.07, 0, 0], ['C', 0, 0, 0], ['N', 1.16, 0, 0]] },
  { key: 'ch3cn', name: 'CH₃CN, acetonitrile', group: 'Organic', start: [['C', 0, 0, 0], ['C', 1.46, 0, 0], ['N', 2.62, 0, 0], ...[0, 120, 240].map((a) => ['H', -0.37, 1.03 * cos(a), 1.03 * sin(a)])] },
  { key: 'hnco', name: 'HNCO, isocyanic acid', group: 'Organic', start: zmat([['N'], ['C', 1, 1.21], ['O', 2, 1.17, 1, 172], ['H', 1, 1.00, 2, 124, 3, 180]]) },
  { key: 'h2co', name: 'H₂CO, formaldehyde', group: 'Organic', start: [['C', 0, 0, 0], ['O', 1.21, 0, 0], ['H', -1.11 * cos(58), 1.11 * sin(58), 0], ['H', -1.11 * cos(58), -1.11 * sin(58), 0]] },
  { key: 'ch3oh', name: 'CH₃OH, methanol', group: 'Organic', start: zmat([['C'], ['O', 1, 1.43], ['H', 2, 0.96, 1, 108], ...methyl(1, 2, 3)]) },
  { key: 'c2h5oh', name: 'C₂H₅OH, ethanol', group: 'Organic', start: zmat([['C'], ['C', 1, 1.52], ['O', 2, 1.43, 1, 108], ['H', 3, 0.96, 2, 108, 1, 180], ...methyl(1, 2, 3),
    ['H', 2, 1.09, 1, 110, 3, 120], ['H', 2, 1.09, 1, 110, 3, -120]]) },
  { key: 'ch3nh2', name: 'CH₃NH₂, methylamine', group: 'Organic', start: zmat([['C'], ['N', 1, 1.47], ['H', 2, 1.01, 1, 110], ['H', 2, 1.01, 1, 110, 3, 115], ...methyl(1, 2, 3, 183)]) },
  { key: 'h2co3', name: 'H₂CO₃, carbonic acid', group: 'Organic', start: zmat([['C'], ['O', 1, 1.21], ['O', 1, 1.34, 2, 125], ['O', 1, 1.34, 2, 125, 3, 180], ['H', 3, 0.97, 1, 107, 2, 0], ['H', 4, 0.97, 1, 107, 2, 0]]) },
  { key: 'hcooh', name: 'HCOOH, formic acid', group: 'Organic', start: zmat([['C'], ['O', 1, 1.20], ['O', 1, 1.34, 2, 125], ['H', 1, 1.10, 2, 125, 3, 180], ['H', 3, 0.97, 1, 107, 2, 0]]) },
  { key: 'ch3cho', name: 'CH₃CHO, acetaldehyde', group: 'Organic', start: zmat([['C'], ['O', 1, 1.21], ['C', 1, 1.50, 2, 124], ['H', 1, 1.11, 2, 120, 3, 180], ...methyl(3, 1, 2, 0)]) },
  { key: 'ch3cooh', name: 'CH₃COOH, acetic acid', group: 'Organic', start: zmat([['C'], ['O', 1, 1.21], ['O', 1, 1.36, 2, 122], ['C', 1, 1.50, 2, 126, 3, 180], ['H', 3, 0.97, 1, 106, 2, 0], ...methyl(4, 1, 2, 0)]) },
  { key: 'urea', name: 'CO(NH₂)₂, urea', group: 'Organic', start: zmat([['C'], ['O', 1, 1.22], ['N', 1, 1.38, 2, 122], ['N', 1, 1.38, 2, 122, 3, 180],
    ['H', 3, 1.01, 1, 119, 2, 0], ['H', 3, 1.01, 1, 119, 2, 180], ['H', 4, 1.01, 1, 119, 2, 0], ['H', 4, 1.01, 1, 119, 2, 180]]) },
  { key: 'glycine', name: 'NH₂CH₂COOH, glycine', group: 'Organic', start: zmat([['N'], ['C', 1, 1.45], ['C', 2, 1.52, 1, 113], ['O', 3, 1.21, 2, 125, 1, 0], ['O', 3, 1.35, 2, 111, 1, 180],
    ['H', 5, 0.97, 3, 107, 2, 180], ['H', 1, 1.01, 2, 110, 3, 60], ['H', 1, 1.01, 2, 110, 3, -60], ['H', 2, 1.09, 1, 109, 3, 120], ['H', 2, 1.09, 1, 109, 3, -120]]) },

  { key: 'c3h6', name: 'C₃H₆, cyclopropane', group: 'Rings', start: [...ring(3, 0.872, 'C'), ...[0.91, -0.91].flatMap((z) => ring(3, 1.452, 'H').map(([e, x, y]) => [e, x, y, z]))] },
  { key: 'c4h5n', name: 'C₄H₅N, pyrrole', group: 'Rings', start: [...ring(5, 1.17, (i) => 'NCCCC'[i]), ...ring(5, 2.22, 'H')] },
  { key: 'c4h4o', name: 'C₄H₄O, furan', group: 'Rings', start: [...ring(5, 1.16, (i) => 'OCCCC'[i]), ...ring(5, 2.22, 'H').slice(1)] },
  { key: 'c6h6', name: 'C₆H₆, benzene', group: 'Rings', start: [...ring(6, 1.39, 'C'), ...ring(6, 2.48, 'H')] },
  { key: 'c5h5n', name: 'C₅H₅N, pyridine', group: 'Rings', start: [...ring(6, 1.39, (i) => 'NCCCCC'[i]), ...ring(6, 2.47, 'H').slice(1)] },

  { key: 'hf', name: 'HF, hydrogen fluoride', group: 'Fluorine and metals', start: zmat([['H'], ['F', 1, 0.92]]) },
  { key: 'f2', name: 'F₂, fluorine', group: 'Fluorine and metals', start: zmat([['F'], ['F', 1, 1.41]]) },
  // BF3 in the plane perpendicular to x, like BH3: a second molecule arrives along its empty orbital
  { key: 'bf3', name: 'BF₃, boron trifluoride', group: 'Fluorine and metals', start: [['B', 0, 0, 0], ...[90, 210, 330].map((a) => ['F', 0, 1.31 * cos(a), 1.31 * sin(a)])] },
  { key: 'cf4', name: 'CF₄, tetrafluoromethane', group: 'Fluorine and metals', start: tetra('C', 1.32, 'F') },
  { key: 'ch3f', name: 'CH₃F, fluoromethane', group: 'Fluorine and metals', start: zmat([['C'], ['F', 1, 1.38], ['H', 1, 1.09, 2, 109], ['H', 1, 1.09, 2, 109, 3, 120], ['H', 1, 1.09, 2, 109, 3, -120]]) },
  { key: 'lih', name: 'LiH, lithium hydride', group: 'Fluorine and metals', start: zmat([['Li'], ['H', 1, 1.60]]), note: 'The method noticeably underestimates the length of this bond: 1.37 Å calculated, about 1.60 Å in experiment.' },
  { key: 'lif', name: 'LiF, lithium fluoride', group: 'Fluorine and metals', start: zmat([['Li'], ['F', 1, 1.56]]) },
  { key: 'li2', name: 'Li₂, dilithium', group: 'Fluorine and metals', start: zmat([['Li'], ['Li', 1, 2.67]]), note: 'The method noticeably underestimates the length of this bond: 2.30 Å calculated, about 2.67 Å in experiment.' },
  { key: 'beh2', name: 'BeH₂, beryllium hydride', group: 'Fluorine and metals', start: [['H', -1.33, 0, 0], ['Be', 0, 0, 0], ['H', 1.33, 0, 0]] },
  { key: 'bef2', name: 'BeF₂, beryllium fluoride', group: 'Fluorine and metals', start: [['F', -1.37, 0, 0], ['Be', 0, 0, 0], ['F', 1.37, 0, 0]] },
  { key: 'nah', name: 'NaH, sodium hydride', group: 'Fluorine and metals', start: zmat([['Na'], ['H', 1, 1.89]]), note: 'The method noticeably underestimates the length of this bond: 1.56 Å calculated, about 1.89 Å in experiment.' },
  { key: 'naf', name: 'NaF, sodium fluoride', group: 'Fluorine and metals', start: zmat([['Na'], ['F', 1, 1.93]]), note: 'The method noticeably underestimates the length of this bond: 1.75 Å calculated, about 1.93 Å in experiment.' },
  { key: 'kf', name: 'KF, potassium fluoride', group: 'Fluorine and metals', start: zmat([['K'], ['F', 1, 2.17]]), note: 'The method noticeably underestimates the length of this bond: 2.02 Å calculated, about 2.17 Å in experiment.' },
  { key: 'lioh', name: 'LiOH, lithium hydroxide', group: 'Fluorine and metals', start: zmat([['Li'], ['O', 1, 1.60], ['H', 2, 0.96, 1, 150]]) },
  { key: 'naoh', name: 'NaOH, caustic soda', group: 'Fluorine and metals', start: zmat([['Na'], ['O', 1, 1.95], ['H', 2, 0.96, 1, 150]]) },
  { key: 'koh', name: 'KOH, caustic potash', group: 'Fluorine and metals', start: zmat([['K'], ['O', 1, 2.20], ['H', 2, 0.96, 1, 150]]) },
  { key: 'he', name: 'He, helium atom', group: 'Fluorine and metals', start: [['He', 0, 0, 0]] },

  // from here on the atoms carry d functions
  { key: 'h2s', name: 'H₂S, hydrogen sulfide', group: 'Third period', start: zmat([['S'], ['H', 1, 1.34], ['H', 1, 1.34, 2, 92]]) },
  { key: 'hcl', name: 'HCl, hydrogen chloride', group: 'Third period', start: zmat([['H'], ['Cl', 1, 1.27]]) },
  { key: 'cl2', name: 'Cl₂, chlorine', group: 'Third period', start: zmat([['Cl'], ['Cl', 1, 1.99]]) },
  { key: 'ph3', name: 'PH₃, phosphine', group: 'Third period', start: zmat([['P'], ['H', 1, 1.42], ['H', 1, 1.42, 2, 93.5], ['H', 1, 1.42, 2, 93.5, 3, 94]]) },
  { key: 'sih4', name: 'SiH₄, silane', group: 'Third period', start: tetra('Si', 1.48) },
  { key: 'ch3cl', name: 'CH₃Cl, chloromethane', group: 'Third period', start: zmat([['C'], ['Cl', 1, 1.78], ['H', 1, 1.09, 2, 108], ['H', 1, 1.09, 2, 108, 3, 120], ['H', 1, 1.09, 2, 108, 3, -120]]) },
  { key: 'ccl4', name: 'CCl₄, carbon tetrachloride', group: 'Third period', start: tetra('C', 1.77, 'Cl') },
  { key: 'so2', name: 'SO₂, sulfur dioxide', group: 'Third period', start: zmat([['S'], ['O', 1, 1.43], ['O', 1, 1.43, 2, 119]]) },
  { key: 'so3', name: 'SO₃, sulfur trioxide', group: 'Third period', start: zmat([['S'], ['O', 1, 1.42], ['O', 1, 1.42, 2, 120], ['O', 1, 1.42, 2, 120, 3, 180]]) },
  { key: 'h2so4', name: 'H₂SO₄, sulfuric acid', group: 'Third period', start: zmat([['S'], ['O', 1, 1.42], ['O', 1, 1.42, 2, 123], ['O', 1, 1.57, 2, 108, 3, 115], ['O', 1, 1.57, 2, 108, 3, -115],
    ['H', 4, 0.97, 1, 108, 2, 30], ['H', 5, 0.97, 1, 108, 3, 30]]) },
  { key: 'h3po4', name: 'H₃PO₄, phosphoric acid', group: 'Third period', start: zmat([['P'], ['O', 1, 1.48], ['O', 1, 1.58, 2, 115], ['O', 1, 1.58, 2, 115, 3, 120], ['O', 1, 1.58, 2, 115, 3, -120],
    ['H', 3, 0.97, 1, 112, 2, 40], ['H', 4, 0.97, 1, 112, 2, 40], ['H', 5, 0.97, 1, 112, 2, 40]]) },
  { key: 'pcl5', name: 'PCl₅, phosphorus pentachloride', group: 'Third period', start: [['P', 0, 0, 0], ['Cl', 0, 0, 2.12], ['Cl', 0, 0, -2.12], ...[0, 120, 240].map((a) => ['Cl', 2.02 * cos(a), 2.02 * sin(a), 0])] },
  { key: 'sf6', name: 'SF₆, sulfur hexafluoride', group: 'Third period', start: [['S', 0, 0, 0], ...octa.map((v) => ['F', ...v.map((x) => 1.56 * x)])] },
  { key: 'nacl', name: 'NaCl, sodium chloride', group: 'Third period', start: zmat([['Na'], ['Cl', 1, 2.36]]) },
  { key: 'mgcl2', name: 'MgCl₂, magnesium chloride', group: 'Third period', start: [['Cl', -2.18, 0, 0], ['Mg', 0, 0, 0], ['Cl', 2.18, 0, 0]] },
  { key: 'alcl3', name: 'AlCl₃, aluminium chloride', group: 'Third period', start: [['Al', 0, 0, 0], ...[90, 210, 330].map((a) => ['Cl', 0, 2.06 * cos(a), 2.06 * sin(a)])] },
  { key: 'ne', name: 'Ne, neon atom', group: 'Third period', start: [['Ne', 0, 0, 0]] },
  { key: 'ar', name: 'Ar, argon atom', group: 'Third period', start: [['Ar', 0, 0, 0]] },

  // d metals: here the d functions are the valence orbitals themselves
  { key: 'ticl4', name: 'TiCl₄, titanium tetrachloride', group: 'd metals', start: tetra('Ti', 2.17, 'Cl'),
    note: 'The five lowest empty orbitals are the d orbitals of titanium: the tetrahedron of chlorines splits them into a pair below and a triple above.' },
  { key: 'nico4', name: 'Ni(CO)₄, nickel carbonyl', group: 'd metals', start: [['Ni', 0, 0, 0], ...tetra('C', 1.84, 'C').slice(1), ...tetra('O', 2.98, 'O').slice(1)],
    note: 'Nickel has all five d orbitals filled: they are the highest occupied orbitals of the molecule.' },
  { key: 'feco5', name: 'Fe(CO)₅, iron pentacarbonyl', group: 'd metals', start: [['Fe', 0, 0, 0],
    ...[1.81, 2.96].flatMap((r, i) => [[i ? 'O' : 'C', 0, 0, r], [i ? 'O' : 'C', 0, 0, -r], ...[0, 120, 240].map((a) => [i ? 'O' : 'C', r * cos(a), r * sin(a), 0])])] },
  { key: 'crco6', name: 'Cr(CO)₆, chromium hexacarbonyl', group: 'd metals', start: [['Cr', 0, 0, 0], ...octa.map((v) => ['C', ...v.map((x) => 1.91 * x)]), ...octa.map((v) => ['O', ...v.map((x) => 3.05 * x)])],
    note: 'The three highest occupied orbitals are the chromium d orbitals that point between the CO groups (dxy, dyz, dxz); the two that point at them (dz², dx²−y²) are empty.' },
  { key: 'ferrocene', name: 'Fe(C₅H₅)₂, ferrocene', group: 'd metals', start: [['Fe', 0, 0, 0],
    ...[1.66, -1.66].flatMap((z) => [...ring(5, 1.21, 'C'), ...ring(5, 2.29, 'H')].map(([el, x, y]) => [el, x, y, z]))],
    note: 'The frontier orbitals are the five d orbitals of iron: three filled (dz², dxy, dx²−y²) and, above them, two empty (dxz, dyz).' },
  { key: 'mno4-', name: 'MnO₄⁻, permanganate', group: 'd metals', charge: -1, start: tetra('Mn', 1.63, 'O') },

  { key: 'cl-', name: 'Cl⁻, chloride', group: 'Ions', charge: -1, start: [['Cl', 0, 0, 0]] },
  { key: 'cl', name: 'Cl·, chlorine atom', group: 'Radicals', start: [['Cl', 0, 0, 0]] },
  { key: 'f-', name: 'F⁻, fluoride', group: 'Ions', charge: -1, start: [['F', 0, 0, 0]] },
  { key: 'bf4-', name: 'BF₄⁻, tetrafluoroborate', group: 'Ions', charge: -1, start: tetra('B', 1.40, 'F') },
  { key: 'h', name: 'H·, hydrogen atom', group: 'Radicals', start: [['H', 0, 0, 0]] },
  { key: 'o', name: 'O, oxygen atom', group: 'Radicals', unpaired: 2, start: [['O', 0, 0, 0]] },
  { key: 'f', name: 'F·, fluorine atom', group: 'Radicals', start: [['F', 0, 0, 0]] },
  { key: 'li', name: 'Li·, lithium atom', group: 'Radicals', start: [['Li', 0, 0, 0]] },
  { key: 'na', name: 'Na·, sodium atom', group: 'Radicals', start: [['Na', 0, 0, 0]] },
  { key: 'oh', name: 'OH·, hydroxyl', group: 'Radicals', start: zmat([['O'], ['H', 1, 0.97]]) },
  { key: 'ho2', name: 'HO₂·, hydroperoxyl', group: 'Radicals', start: zmat([['O'], ['O', 1, 1.33], ['H', 1, 0.97, 2, 104]]) },
  { key: 'ch3', name: 'CH₃·, methyl', group: 'Radicals', start: [['C', 0, 0, 0], ...[90, 210, 330].map((a) => ['H', 0, 1.08 * cos(a), 1.08 * sin(a)])] },
  { key: 'nh2', name: 'NH₂·, amino radical', group: 'Radicals', start: zmat([['N'], ['H', 1, 1.02], ['H', 1, 1.02, 2, 103]]) },
  { key: 'no', name: 'NO·, nitric oxide', group: 'Radicals', start: zmat([['N'], ['O', 1, 1.15]]) },
  { key: 'no2', name: 'NO₂·, nitrogen dioxide', group: 'Radicals', start: zmat([['N'], ['O', 1, 1.20], ['O', 1, 1.20, 2, 134]]) },

  { key: 'hco3-', name: 'HCO₃⁻, bicarbonate', group: 'Ions', charge: -1, start: zmat([['C'], ['O', 1, 1.25], ['O', 1, 1.25, 2, 130], ['O', 1, 1.42, 2, 115, 3, 180], ['H', 4, 0.97, 1, 104, 2, 0]]) },
  { key: 'hcoo-', name: 'HCOO⁻, formate', group: 'Ions', charge: -1, start: zmat([['C'], ['O', 1, 1.25], ['O', 1, 1.25, 2, 130], ['H', 1, 1.12, 2, 115, 3, 180]]) },
  { key: 'ch3coo-', name: 'CH₃COO⁻, acetate', group: 'Ions', charge: -1, start: zmat([['C'], ['O', 1, 1.26], ['O', 1, 1.26, 2, 128], ['C', 1, 1.54, 2, 116, 3, 180], ...methyl(4, 1, 2, 0)]) },
  { key: 'no3-', name: 'NO₃⁻, nitrate', group: 'Ions', charge: -1, start: zmat([['N'], ['O', 1, 1.25], ['O', 1, 1.25, 2, 120], ['O', 1, 1.25, 2, 120, 3, 180]]) },
  { key: 'h3o+', name: 'H₃O⁺, hydronium', group: 'Ions', charge: 1, start: zmat([['O'], ['H', 1, 0.98], ['H', 1, 0.98, 2, 112], ['H', 1, 0.98, 2, 112, 3, 125]]) },
  { key: 'nh4+', name: 'NH₄⁺, ammonium', group: 'Ions', charge: 1, start: tetra('N', 1.02) },
  { key: 'oh-', name: 'OH⁻, hydroxide', group: 'Ions', charge: -1, start: zmat([['O'], ['H', 1, 0.97]]) },
  { key: 'cn-', name: 'CN⁻, cyanide', group: 'Ions', charge: -1, start: zmat([['C'], ['N', 1, 1.18]]) },
  { key: 'bh4-', name: 'BH₄⁻, borohydride', group: 'Ions', charge: -1, start: tetra('B', 1.24) },
];

let unpairedNow; // number of unpaired electrons of the molecule being built (undefined: lowest possible)
function energyAndGradient(numbers, x, charge = 0) {
  const r = xtb.compute({ numbers, positions: x, charge, etemp: 300, forces: true, unpaired: unpairedNow });
  if (!r.ok) throw new Error(r.error);
  return { e: r.energy, g: r.forces.map((f) => -f), r };
}

// BFGS in Cartesian coordinates with a step cap and backtracking; x in bohr
function optimize(numbers, x0, charge = 0, gtol = 2e-5, maxIter = 800) {
  const n = x0.length;
  let x = Float64Array.from(x0), { e, g } = energyAndGradient(numbers, x, charge);
  let H = Array.from({ length: n }, (_, i) => { const row = new Float64Array(n); row[i] = 2; return row; }); // inverse Hessian guess
  for (let it = 0; it < maxIter; it++) {
    const gmax = g.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    if (gmax < gtol) return { x, e, gmax, it };
    let d = new Float64Array(n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) d[i] -= H[i][j] * g[j];
    if (d.reduce((s, v, i) => s + v * g[i], 0) >= 0) d = g.map((v) => -v); // not a descent direction
    const dmax = d.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    let scale = Math.min(1, 0.3 / dmax), next, xn;
    for (;;) {
      xn = x.map((v, i) => v + scale * d[i]);
      next = energyAndGradient(numbers, xn, charge);
      if (next.e <= e + 1e-10 || scale < 1e-4) break;
      scale *= 0.4;
    }
    const s = xn.map((v, i) => v - x[i]), y = next.g.map((v, i) => v - g[i]);
    const sy = s.reduce((t, v, i) => t + v * y[i], 0);
    if (sy > 1e-12) {
      const Hy = new Float64Array(n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Hy[i] += H[i][j] * y[j];
      const yHy = y.reduce((t, v, i) => t + v * Hy[i], 0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++)
        H[i][j] += ((sy + yHy) * s[i] * s[j]) / (sy * sy) - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
    }
    x = xn; e = next.e; g = next.g;
  }
  throw new Error('geometry optimization did not converge');
}

const dist = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]) * BOHR_ANGSTROM;
const out = [];
for (const m of molecules) {
  const elements = m.start.map((a) => a[0]), numbers = elements.map((el) => NUMBERS[el]);
  const x0 = Float64Array.from(m.start.flatMap((a) => a.slice(1).map((v) => v / BOHR_ANGSTROM)));
  let x = x0, info;
  unpairedNow = m.unpaired;
  if (m.fixed) {
    const { e, g } = energyAndGradient(numbers, x0, m.charge || 0);
    info = { e, gmax: g.reduce((mx, v) => Math.max(mx, Math.abs(v)), 0), it: 0 };
  } else {
    info = optimize(numbers, x0, m.charge || 0);
    x = info.x;
  }
  // centre on the centroid, as the reference demo does
  const c = [0, 1, 2].map((k) => x.reduce((s, v, i) => s + (i % 3 === k ? v : 0), 0) / numbers.length);
  const xyz = Array.from(x, (v, i) => +((v - c[i % 3]) * BOHR_ANGSTROM).toFixed(6));
  out.push({ key: m.key, name: m.name, group: m.group, charge: m.charge || 0, unpaired: m.unpaired, elements, xyz, optimized: !m.fixed, energy: info.e, note: m.note });
  // sanity: one piece (unless meant otherwise), and how far the optimizer moved the guess
  const parts = pieces({ elements, positions: x }).count;
  let moved = 0;
  for (let i = 0; i < x.length; i += 3) moved = Math.max(moved, Math.hypot(x[i] - x0[i], x[i + 1] - x0[i + 1], x[i + 2] - x0[i + 2]) * BOHR_ANGSTROM);
  const bonds = [];
  for (let i = 0; i < numbers.length; i++) for (let j = i + 1; j < numbers.length; j++) { const d = dist(x, i, j); if (d < 1.75) bonds.push(`${elements[i]}${elements[j]} ${d.toFixed(2)}`); }
  console.log(`${m.key.padEnd(8)} E = ${info.e.toFixed(6)}  max|g| ${info.gmax.toExponential(1)}  steps ${String(info.it).padStart(3)}  pieces ${parts}${parts !== 1 && !m.fixed ? ' !!!' : ''}  moved ${moved.toFixed(2)} A  | ${[...new Set(bonds)].join(', ')}`);
}
writeFileSync(new URL('app/src/library.json', root), JSON.stringify({ comment: 'generated by tools/build-library.mjs; xyz in angstrom', molecules: out }, null, 1) + '\n');
