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

// The orbital a scene is first shown with: the highest filled one; with no valence electrons
// at all (a bare K+ or Na+), the lowest empty one.
export function firstShown(result) {
  const { homo, lumo } = frontier(result);
  return homo >= 0 ? homo : lumo;
}

// Functions of a shell, in the order of the engine interface. Each is the radial part of the
// shell times a polynomial of the displacement from the nucleus:
//   s;  p: x, y, z;  d: sqrt(3) xy, sqrt(3) yz, z^2 - (x^2 + y^2)/2, sqrt(3) xz, sqrt(3)/2 (x^2 - y^2).
// The contraction coefficients of a shell carry the normalization of its x^l function, and
// these factors make the other functions of the shell normalized too.
const S3 = Math.sqrt(3);
export const AO_LABELS = [['s'], ['px', 'py', 'pz'], ['dxy', 'dyz', 'dz²', 'dxz', 'dx²−y²']];
// Cartesian monomials x^i y^j z^k of a shell, and the functions above as their combinations
export const CART = [[[0, 0, 0]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]];
export const COMBO = [[[1]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  [[0, S3, 0, 0, 0, 0], [0, 0, 0, 0, S3, 0], [-0.5, 0, 0, -0.5, 0, 1], [0, 0, S3, 0, 0, 0], [S3 / 2, 0, 0, -S3 / 2, 0, 0]]];

// Values of the 2l + 1 functions of a shell at displacement (x, y, z) from its nucleus,
// written to out[at...]. The reference the shader is checked against.
export function shellValues(sh, x, y, z, out, at = 0) {
  const r2 = x * x + y * y + z * z;
  let R = 0;
  for (let k = 0; k < sh.alpha.length; k++) R += sh.coeff[k] * Math.exp(-sh.alpha[k] * r2);
  if (sh.l === 0) out[at] = R;
  else if (sh.l === 1) { out[at] = R * x; out[at + 1] = R * y; out[at + 2] = R * z; }
  else {
    out[at] = S3 * R * x * y; out[at + 1] = S3 * R * y * z; out[at + 2] = R * (z * z - 0.5 * (x * x + y * y));
    out[at + 3] = S3 * R * x * z; out[at + 4] = 0.5 * S3 * R * (x * x - y * y);
  }
  return 2 * sh.l + 1;
}

// Human-readable makeup of orbital m: the largest AO coefficients.
export function describe(result, elements, m, threshold = 0.25) {
  const labels = [];
  for (const sh of result.shells)
    for (const t of AO_LABELS[sh.l]) labels.push(`${elements[sh.atom]}${sh.atom}${t}`);
  const c = result.coeff.subarray(m * result.nao, (m + 1) * result.nao);
  return Array.from(c, (v, i) => [v, labels[i]]).filter(([v]) => Math.abs(v) > threshold);
}

// Overlap matrix of the basis functions (nao x nao, row-major) from the shells and the atom
// positions (bohr): analytic integrals of contracted s, p and d Gaussians, shell pair by
// shell pair. A primitive integral is a product of three one-dimensional ones, which come
// from the Obara–Saika recurrence.
export function overlapMatrix(shells, positions) {
  const first = [];
  let n = 0;
  for (const sh of shells) { first.push(n); n += 2 * sh.l + 1; }
  const S = new Float64Array(n * n), cart = new Float64Array(36);
  const X = [new Float64Array(9), new Float64Array(9), new Float64Array(9)]; // X[axis][3 i + j]: powers i on A, j on B
  for (let sa = 0; sa < shells.length; sa++) for (let sb = 0; sb <= sa; sb++) {
    const A = shells[sa], B = shells[sb], la = A.l, lb = B.l, ca = CART[la], cb = CART[lb];
    const pa = positions.subarray(3 * A.atom, 3 * A.atom + 3), pb = positions.subarray(3 * B.atom, 3 * B.atom + 3);
    const ab2 = (pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2 + (pa[2] - pb[2]) ** 2;
    cart.fill(0);
    for (let p = 0; p < A.alpha.length; p++) for (let q = 0; q < B.alpha.length; q++) {
      const a = A.alpha[p], b = B.alpha[q], g = a + b, h = 0.5 / g;
      const f = A.coeff[p] * B.coeff[q] * Math.pow(Math.PI / g, 1.5) * Math.exp(-a * b / g * ab2);
      for (let d = 0; d < 3; d++) {
        const c = (a * pa[d] + b * pb[d]) / g, PA = c - pa[d], PB = c - pb[d], t = X[d];
        t[0] = 1;
        for (let i = 1; i <= la; i++) t[3 * i] = PA * t[3 * i - 3] + (i > 1 ? (i - 1) * h * t[3 * i - 6] : 0);
        for (let i = 0; i <= la; i++) for (let j = 1; j <= lb; j++)
          t[3 * i + j] = PB * t[3 * i + j - 1] + (i ? i * h * t[3 * i - 3 + j - 1] : 0) + (j > 1 ? (j - 1) * h * t[3 * i + j - 2] : 0);
      }
      for (let u = 0; u < ca.length; u++) for (let v = 0; v < cb.length; v++)
        cart[6 * u + v] += f * X[0][3 * ca[u][0] + cb[v][0]] * X[1][3 * ca[u][1] + cb[v][1]] * X[2][3 * ca[u][2] + cb[v][2]];
    }
    const ka = COMBO[la], kb = COMBO[lb];
    for (let i = 0; i < ka.length; i++) for (let j = 0; j < kb.length; j++) {
      let s = 0;
      for (let u = 0; u < ca.length; u++) if (ka[i][u]) for (let v = 0; v < cb.length; v++) s += ka[i][u] * kb[j][v] * cart[6 * u + v];
      S[(first[sa] + i) * n + first[sb] + j] = S[(first[sb] + j) * n + first[sa] + i] = s;
    }
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
