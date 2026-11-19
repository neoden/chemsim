// Extended Hückel engine. The numerical core is the one from reference/eht-core.js
// (STO-3G-like minimal basis, YAeHMOP parameters, K = 1.75), wrapped in the common
// engine interface. No forces: it stays as the reference for checking the renderer.
import { HARTREE_EV, SYMBOLS } from '../units.js';

const G1 = { a: [2.227660584, 0.405771156, 0.109818], s: [0.154328967, 0.535328142, 0.444634542] };
const G2 = { a: [0.994203, 0.231031, 0.0751386], s: [-0.09996723, 0.39951283, 0.70011547], p: [0.15591627, 0.60768372, 0.39195739] };
const EL = {
  H: { z: 1.3, Hs: -13.6, ve: 1 }, B: { z: 1.3, Hs: -15.2, Hp: -8.5, ve: 3 }, C: { z: 1.625, Hs: -21.4, Hp: -11.4, ve: 4 },
  N: { z: 1.95, Hs: -26.0, Hp: -13.4, ve: 5 }, O: { z: 2.275, Hs: -32.3, Hp: -14.8, ve: 6 },
};

// atomic numbers this engine has parameters for
export const EHT_NUMBERS = [1, 5, 6, 7, 8];

// r in bohr
function mkAtom(el, r) {
  const E = EL[el], g = el === 'H' ? G1 : G2, al = g.a.map((v) => v * E.z * E.z);
  let ds = g.s.map((d, i) => d * Math.pow(2 * al[i] / Math.PI, 0.75));
  let dp = g.p ? g.p.map((d, i) => d * Math.pow(2 * al[i] / Math.PI, 0.75) * 2 * Math.sqrt(al[i])) : [0, 0, 0];
  let ns = 0, np = 0;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    const p = al[i] + al[j], b = Math.pow(Math.PI / p, 1.5);
    ns += ds[i] * ds[j] * b; np += dp[i] * dp[j] * b / (2 * p);
  }
  ds = ds.map((v) => v / Math.sqrt(ns)); if (np > 0) dp = dp.map((v) => v / Math.sqrt(np));
  return { el, r, al, ds, dp, E };
}
function basis(atoms) {
  const b = [];
  atoms.forEach((a, k) => { b.push({ k, t: 0 }); if (a.el !== 'H') { b.push({ k, t: 1 }); b.push({ k, t: 2 }); b.push({ k, t: 3 }); } });
  return b;
}
function ovl(A, ta, B, tb) {
  let s = 0; const ca = ta ? A.dp : A.ds, cb = tb ? B.dp : B.ds;
  let ab2 = 0; for (let d = 0; d < 3; d++) ab2 += (A.r[d] - B.r[d]) ** 2;
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    const a = A.al[i], b = B.al[j], p = a + b;
    let f = Math.pow(Math.PI / p, 1.5) * Math.exp(-a * b / p * ab2);
    for (let d = 0; d < 3; d++) {
      const P = (a * A.r[d] + b * B.r[d]) / p, pa = P - A.r[d], pb = P - B.r[d], la = ta === d + 1, lb = tb === d + 1;
      if (la && lb) f *= pa * pb + 0.5 / p; else if (la) f *= pa; else if (lb) f *= pb;
    }
    s += ca[i] * cb[j] * f;
  }
  return s;
}
function jacobi(M, n) {
  const a = M.map((r) => r.slice()), v = [];
  for (let i = 0; i < n; i++) { v.push(new Array(n).fill(0)); v[i][i] = 1; }
  for (let sw = 0; sw < 80; sw++) {
    let off = 0; for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += a[p][q] * a[p][q];
    if (off < 1e-22) break;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
      for (let k = 0; k < n; k++) { const x = v[k][p], y = v[k][q]; v[k][p] = c * x - s * y; v[k][q] = s * x + c * y; }
    }
  }
  return { w: a.map((r, i) => r[i]), v };
}
function mm(A, B, n) {
  const C = [];
  for (let i = 0; i < n; i++) { C.push([]); for (let j = 0; j < n; j++) { let t = 0; for (let k = 0; k < n; k++) t += A[i][k] * B[k][j]; C[i][j] = t; } }
  return C;
}
export function solve(atoms) {
  const bs = basis(atoms), n = bs.length, S = [], H = [];
  for (let i = 0; i < n; i++) { S.push([]); H.push([]); }
  const hd = bs.map((b) => (b.t ? atoms[b.k].E.Hp : atoms[b.k].E.Hs));
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    S[i][j] = i === j ? 1 : ovl(atoms[bs[i].k], bs[i].t, atoms[bs[j].k], bs[j].t);
    H[i][j] = i === j ? hd[i] : 1.75 * S[i][j] * (hd[i] + hd[j]) / 2;
  }
  const es = jacobi(S, n), X = [];
  for (let i = 0; i < n; i++) { X.push([]); for (let j = 0; j < n; j++) { let t = 0; for (let k = 0; k < n; k++) t += es.v[i][k] * es.v[j][k] / Math.sqrt(es.w[k]); X[i][j] = t; } }
  const eh = jacobi(mm(mm(X, H, n), X, n), n), C = mm(X, eh.v, n);
  const ord = eh.w.map((e, i) => i).sort((a, b) => eh.w[a] - eh.w[b]);
  return { bs, n, S, H, E: ord.map((i) => eh.w[i]), C: ord.map((i) => C.map((r) => r[i])), nocc: atoms.reduce((s, a) => s + a.E.ve, 0) / 2 };
}

export function createEhtEngine() {
  return {
    name: 'eht',
    hasForces: false,
    compute(req) {
      const other = req.numbers.find((z) => !EHT_NUMBERS.includes(z));
      if (other !== undefined) return { ok: false, converged: true, error: `the Hückel method has no parameters here for the element ${SYMBOLS[other] ?? other}` };
      const nel = req.numbers.reduce((s, z) => s + EL[SYMBOLS[z]].ve, 0) - (req.charge || 0);
      if (nel % 2) return { ok: false, converged: true, error: 'The Hückel method here handles only molecules with an even number of electrons' };
      const atoms = req.numbers.map((z, k) => mkAtom(SYMBOLS[z], [req.positions[3 * k], req.positions[3 * k + 1], req.positions[3 * k + 2]]));
      const s = solve(atoms), n = s.n, nocc = nel / 2;
      const coeff = new Float64Array(n * n), eps = new Float64Array(n), occ = new Float64Array(n);
      let energy = 0;
      for (let m = 0; m < n; m++) {
        eps[m] = s.E[m] / HARTREE_EV;
        occ[m] = m < nocc ? 2 : 0;
        energy += occ[m] * eps[m];
        coeff.set(s.C[m], m * n);
      }
      const shells = [];
      atoms.forEach((a, k) => {
        shells.push({ atom: k, l: 0, alpha: Float64Array.from(a.al), coeff: Float64Array.from(a.ds) });
        if (a.el !== 'H') shells.push({ atom: k, l: 1, alpha: Float64Array.from(a.al), coeff: Float64Array.from(a.dp) });
      });
      // energy: sum of orbital energies, the only "total energy" the method has
      return { ok: true, converged: true, engine: 'eht', energy, forces: null, eps, occ, coeff, nao: n, nmo: n, shells, iterations: 0 };
    },
  };
}
