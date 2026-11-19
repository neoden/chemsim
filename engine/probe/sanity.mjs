// Checks that the basis description coming out of the engine is interpreted correctly:
// the density built from shells + coefficients must integrate to the valence electron count.
import { readFileSync } from 'node:fs';
import { loadXtb, parseXyz, unpack } from './xtb.mjs';
const root = new URL('../', import.meta.url);
const xtb = await loadXtb(readFileSync(new URL('xtb-wasm/target/wasm32-wasip1/release/xtb_wasm.wasm', root)));
const HARTREE_EV = 27.211386245988;
for (const [name, nval] of [['h2o', 8], ['bh3_nh3', 14]]) {
  const mol = parseXyz(readFileSync(new URL(`cases/${name}.xyz`, root), 'utf8'));
  const r = unpack(xtb.create(mol.numbers, mol.pos).singlepoint(mol.pos).raw);
  // AO list: p order is y, z, x
  const aos = [];
  for (const sh of r.shells) {
    const c = mol.pos.subarray(3 * sh.atom, 3 * sh.atom + 3);
    if (sh.l === 0) aos.push({ c, sh, k: -1 });
    else for (const k of [1, 2, 0]) aos.push({ c, sh, k });
  }
  const labels = aos.map(({ sh, k }) => `${'?H???BCNO'[mol.numbers[sh.atom]]}${sh.atom}${k < 0 ? 's' : 'p' + 'xyz'[k]}`);
  const h = 0.12, L = 9, n = Math.round(2 * L / h);
  let total = 0;
  const ao = new Float64Array(aos.length);
  const cx = mol.pos.filter((_, i) => i % 3 === 0).reduce((a, b) => a + b) / mol.numbers.length;
  const cz = mol.pos.filter((_, i) => i % 3 === 2).reduce((a, b) => a + b) / mol.numbers.length;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let k = 0; k < n; k++) {
    const p = [cx - L + (i + 0.5) * h, -L + (j + 0.5) * h, cz - L + (k + 0.5) * h];
    aos.forEach(({ c, sh, k }, m) => {
      const d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]], r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      let v = 0;
      for (let q = 0; q < sh.alpha.length; q++) v += sh.coeff[q] * Math.exp(-sh.alpha[q] * r2);
      ao[m] = k < 0 ? v : v * d[k];
    });
    for (let o = 0; o < r.nmo; o++) {
      if (r.occ[o] < 1e-12) continue;
      let psi = 0;
      for (let m = 0; m < r.nao; m++) psi += r.coeff[o * r.nao + m] * ao[m];
      total += r.occ[o] * psi * psi;
    }
  }
  console.log(`${name}: electrons on grid = ${(total * h ** 3).toFixed(4)} (expected ${nval}); primitives per shell: ${[...new Set(r.shells.map((s) => s.alpha.length))]}`);
  r.eps.forEach((e, o) => {
    const top = [...r.coeff.subarray(o * r.nao, (o + 1) * r.nao)].map((c, m) => [c, labels[m]]).filter(([c]) => Math.abs(c) > 0.25)
      .map(([c, l]) => `${l}:${c.toFixed(2)}`).join(' ');
    console.log(`  ${(e * HARTREE_EV).toFixed(3).padStart(8)} eV  occ ${r.occ[o].toFixed(3)}  ${top}`);
  });
}
