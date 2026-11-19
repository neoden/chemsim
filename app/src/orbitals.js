// Orbital bookkeeping on top of an engine result.

// The sign of an eigenvector is arbitrary. For a freshly loaded molecule fix it by making
// the largest coefficient positive (the rule of the reference demo).
export function fixSigns(result) {
  const { nao, nmo, coeff } = result;
  for (let m = 0; m < nmo; m++) {
    const c = coeff.subarray(m * nao, (m + 1) * nao);
    let big = 0;
    for (let i = 1; i < nao; i++) if (Math.abs(c[i]) > Math.abs(c[big]) + 1e-9) big = i;
    if (c[big] < 0) for (let i = 0; i < nao; i++) c[i] = -c[i];
  }
}

// Index of the highest orbital that is at least half filled, and of the one above it.
export function frontier(result) {
  let homo = -1;
  for (let m = 0; m < result.nmo; m++) if (result.occ[m] > 0.999) homo = m; // 0.999: a single electron comes out as 0.99999…
  return { homo, lumo: homo + 1 < result.nmo ? homo + 1 : -1 };
}

// Human-readable makeup of orbital m: the largest AO coefficients.
export function describe(result, elements, m, threshold = 0.25) {
  const labels = [];
  for (const sh of result.shells)
    for (const t of sh.l ? ['px', 'py', 'pz'] : ['s']) labels.push(`${elements[sh.atom]}${sh.atom}${t}`);
  const c = result.coeff.subarray(m * result.nao, (m + 1) * result.nao);
  return Array.from(c, (v, i) => [v, labels[i]]).filter(([v]) => Math.abs(v) > threshold);
}

// Overlap matrix of the basis functions (nao x nao, row-major) from the shells and the atom
// positions (bohr): analytic integrals of contracted s and p Gaussians.
export function overlapMatrix(shells, positions) {
  const fns = []; // one entry per AO: shell, centre, component (-1 for s, 0..2 for p)
  for (const sh of shells) {
    const c = [positions[3 * sh.atom], positions[3 * sh.atom + 1], positions[3 * sh.atom + 2]];
    for (const k of sh.l ? [0, 1, 2] : [-1]) fns.push({ sh, c, k });
  }
  const n = fns.length, S = new Float64Array(n * n);
  for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
    const A = fns[i], B = fns[j];
    const ab2 = (A.c[0] - B.c[0]) ** 2 + (A.c[1] - B.c[1]) ** 2 + (A.c[2] - B.c[2]) ** 2;
    let s = 0;
    for (let p = 0; p < A.sh.alpha.length; p++) for (let q = 0; q < B.sh.alpha.length; q++) {
      const a = A.sh.alpha[p], b = B.sh.alpha[q], g = a + b;
      let f = Math.pow(Math.PI / g, 1.5) * Math.exp(-a * b / g * ab2);
      if (A.k >= 0 || B.k >= 0) {
        const pa = A.k >= 0 ? (a * A.c[A.k] + b * B.c[A.k]) / g - A.c[A.k] : 1;
        const pb = B.k >= 0 ? (a * A.c[B.k] + b * B.c[B.k]) / g - B.c[B.k] : 1;
        f *= pa * pb + (A.k >= 0 && A.k === B.k ? 0.5 / g : 0);
      }
      s += A.sh.coeff[p] * B.sh.coeff[q] * f;
    }
    S[i * n + j] = S[j * n + i] = s;
  }
  return S;
}

// Continuity of orbitals between two consecutive results of the same molecule. The sign of
// an eigenvector is arbitrary and levels cross as the nuclei move, so each new orbital is
// matched to the previous one it overlaps most (prev . S . next), its sign is aligned with
// it, and the selected orbital follows its match rather than its index.
// Changes next.coeff in place; returns the new index of `selected`.
export function followOrbitals(prev, next, S, selected) {
  const { nao, nmo } = next;
  if (!prev || prev.nao !== nao || prev.nmo !== nmo) { fixSigns(next); return Math.min(selected, nmo - 1); }
  // SC[j] = S . next_j
  const SC = new Float64Array(nmo * nao);
  for (let j = 0; j < nmo; j++) for (let a = 0; a < nao; a++) {
    let t = 0;
    for (let b = 0; b < nao; b++) t += S[a * nao + b] * next.coeff[j * nao + b];
    SC[j * nao + a] = t;
  }
  const O = new Float64Array(nmo * nmo); // O[i, j] = <prev_i | next_j>
  const pairs = [];
  for (let i = 0; i < nmo; i++) for (let j = 0; j < nmo; j++) {
    let t = 0;
    for (let a = 0; a < nao; a++) t += prev.coeff[i * nao + a] * SC[j * nao + a];
    O[i * nmo + j] = t;
    pairs.push(i * nmo + j);
  }
  // one-to-one matching, strongest overlaps first
  pairs.sort((p, q) => Math.abs(O[q]) - Math.abs(O[p]));
  const partner = new Int32Array(nmo).fill(-1), taken = new Uint8Array(nmo);
  for (const p of pairs) {
    const i = (p / nmo) | 0, j = p % nmo;
    if (partner[i] >= 0 || taken[j]) continue;
    partner[i] = j; taken[j] = 1;
    if (O[p] < 0) for (let a = 0; a < nao; a++) next.coeff[j * nao + a] = -next.coeff[j * nao + a];
  }
  const follow = partner[selected];
  return follow;
}
