// The electrostatic potential around a molecule, from its valence electrons and the charges
// of its cores (the nuclei with the inner shells, which GFN2-xTB does not have). No DOM.
//
// The electron density is a sum of products of basis functions. Every product is replaced by
// the point charge, point dipole and point quadrupole it has, placed
//   on the nucleus, for two functions of one atom (together with the core's charge);
//   halfway between the nuclei, for functions of two atoms that share a noticeable charge
//   (a bond); for the other pairs of atoms the little they share is split between the two.
// Point charges forget that a cloud has a size: seen from within it, part of it is already
// behind, and the potential is higher than the points give. That is put back for every atom
// as if its shells were round, each with the electrons it holds (`pen` below).
// On the surface the app colours, this is the potential of the density itself to within some
// 0.003 Eh; see tests/esp.test.mjs. The engine's own multipoles (result.charges, .dipoles,
// .quadrupoles) are the same thing with the bonds' share put on the atoms as well, and that is
// too coarse: it paints the rim of benzene, between the hydrogens, as negative.
import { CART, COMBO } from './orbitals.js';
import { VALENCE, COVALENT, BOHR_ANGSTROM } from './units.js';

// The potential at which the colours are full on the scale that is the same for every molecule,
// Eh per unit charge (3 V). Neutral molecules with covalent bonds stay within it (benzene
// 0.5 V, water 2.2 V, HF 3 V); salts such as LiF reach 4 V and ions 6 to 10 V, off the scale.
export const ESP_FIXED_SCALE = 3 / 27.211386245988;
// Density of the surface the potential is shown on, e/bohr^3: about where molecules touch.
export const ESP_DENSITY = 0.002;
// Closer than this to a site (bohr) the expansion means nothing; the shader stops there too.
export const ESP_CORE = 0.7;
// Two atoms get a site of their own between them if they share at least this many electrons.
const BOND_SHARE = 0.02;
// The correction for the size of an atom's cloud is tabulated at PEN_N distances from 0 to PEN_R bohr.
export const PEN_N = 64, PEN_R = 12;

// For one shell, at the PEN_N distances: g, r times (1/r minus the potential of one electron
// spread like the shell, evenly over the directions), which falls from 1 at the nucleus to 0;
// and rho, the density of that electron.
const shellTables = new Map();
function shellTable(sh) {
  const key = `${sh.l}|${sh.alpha.join()}|${sh.coeff.join()}`;
  let g = shellTables.get(key);
  if (g) return g;
  const each = 40, n = (PEN_N - 1) * each, h = PEN_R / n, rho = new Float64Array(n + 1);
  for (let i = 0; i <= n; i++) {
    const r = i * h;
    let R = 0;
    for (let k = 0; k < sh.alpha.length; k++) R += sh.coeff[k] * Math.exp(-sh.alpha[k] * r * r);
    rho[i] = R * R * r ** (2 * sh.l);
  }
  // inside[i]: the charge within r_i; outer[i]: the potential at r_i of the charge beyond it
  const inside = new Float64Array(n + 1), outer = new Float64Array(n + 1);
  for (let i = 1; i <= n; i++) inside[i] = inside[i - 1] + 0.5 * h * (rho[i] * (i * h) ** 2 + rho[i - 1] * ((i - 1) * h) ** 2);
  for (let i = n - 1; i >= 0; i--) outer[i] = outer[i + 1] + 0.5 * h * (rho[i] * i * h + rho[i + 1] * (i + 1) * h);
  g = { g: new Float64Array(PEN_N), rho: new Float64Array(PEN_N) };
  for (let k = 0; k < PEN_N; k++) {
    g.g[k] = 1 - (inside[k * each] + k * each * h * outer[k * each]) / inside[n];
    g.rho[k] = rho[k * each] / (4 * Math.PI * inside[n]);
  }
  shellTables.set(key, g);
  return g;
}

// the monomials whose integrals make a charge, a dipole and the second moments
const MOMENTS = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0], [1, 1, 0], [0, 2, 0], [1, 0, 1], [0, 1, 1], [0, 0, 2]];

// The sites of a molecule: { n, atoms, positions (bohr, 3n), multipoles (10n), pen, reach }. The
// first `atoms` sites are the nuclei, in their order. A site's ten numbers are its charge, its
// dipole (x, y, z) and its traceless quadrupole (xx, xy, yy, xz, yz, zz; sum over the charges
// of q (3 x_i x_j - r^2 delta_ij) / 2), in atomic units. pen: for every atom PEN_N numbers,
// r times what the size of its cloud adds to the potential at distance r. reach: the distance
// from the origin within which the surface the potential is shown on lies.
// bondShare: Infinity leaves no sites between the atoms, which makes the multipoles those of
// the engine (the checks compare the two).
export function chargeSites(numbers, positions, result, { bondShare = BOND_SHARE } = {}) {
  const { shells, nao, nmo, coeff, occ } = result, nat = numbers.length;
  // occupied orbitals times sqrt(occupation), function by function: P_ij = sum_k W[i][k] W[j][k]
  const occupied = [];
  for (let m = 0; m < nmo; m++) if (occ[m] > 1e-9) occupied.push(m);
  const no = occupied.length, W = new Float64Array(nao * no);
  occupied.forEach((m, k) => { const w = Math.sqrt(occ[m]); for (let i = 0; i < nao; i++) W[i * no + k] = w * coeff[m * nao + i]; });
  const first = [], byAtom = Array.from({ length: nat }, () => []);
  let n = 0;
  // softest: the smallest exponent of an atom, to tell that two atoms are too far apart to share anything
  const softest = new Float64Array(nat).fill(Infinity), held = new Float64Array(shells.length); // held: electrons of a shell
  shells.forEach((sh, s) => {
    first.push(n); n += 2 * sh.l + 1; byAtom[sh.atom].push(s);
    softest[sh.atom] = Math.min(softest[sh.atom], ...sh.alpha);
  });

  const X = [new Float64Array(15), new Float64Array(15), new Float64Array(15)]; // X[axis][3 i + j]: powers i on A (two more than the shell has, for the moments), j on B
  const D = new Float64Array(36), m = new Float64Array(10);
  // adds to m the moments about nucleus A of the density that shells sa (on A) and sb share
  const shellPair = (sa, sb, weight) => {
    const A = shells[sa], B = shells[sb], la = A.l, lb = B.l, ca = CART[la], cb = CART[lb], ka = COMBO[la], kb = COMBO[lb];
    const pa = positions.subarray(3 * A.atom, 3 * A.atom + 3), pb = positions.subarray(3 * B.atom, 3 * B.atom + 3);
    const ab2 = (pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2 + (pa[2] - pb[2]) ** 2;
    // the density matrix of the pair, in Cartesian monomials
    D.fill(0);
    for (let i = 0; i < ka.length; i++) for (let j = 0; j < kb.length; j++) {
      let p = 0;
      for (let k = 0, wi = (first[sa] + i) * no, wj = (first[sb] + j) * no; k < no; k++) p += W[wi + k] * W[wj + k];
      for (let u = 0; u < ca.length; u++) if (ka[i][u]) for (let v = 0; v < cb.length; v++) D[6 * u + v] += p * ka[i][u] * kb[j][v];
    }
    for (let p = 0; p < A.alpha.length; p++) for (let q = 0; q < B.alpha.length; q++) {
      const a = A.alpha[p], b = B.alpha[q], g = a + b, h = 0.5 / g, damp = a * b / g * ab2;
      if (damp > 32) continue;
      const f = weight * A.coeff[p] * B.coeff[q] * Math.pow(Math.PI / g, 1.5) * Math.exp(-damp);
      for (let d = 0; d < 3; d++) {
        const c = (a * pa[d] + b * pb[d]) / g, PA = c - pa[d], PB = c - pb[d], t = X[d];
        t[0] = 1;
        for (let i = 1; i <= la + 2; i++) t[3 * i] = PA * t[3 * i - 3] + (i > 1 ? (i - 1) * h * t[3 * i - 6] : 0);
        for (let i = 0; i <= la + 2; i++) for (let j = 1; j <= lb; j++)
          t[3 * i + j] = PB * t[3 * i + j - 1] + (i ? i * h * t[3 * i - 3 + j - 1] : 0) + (j > 1 ? (j - 1) * h * t[3 * i + j - 2] : 0);
      }
      for (let u = 0; u < ca.length; u++) for (let v = 0; v < cb.length; v++) {
        const d = f * D[6 * u + v];
        if (d === 0) continue;
        const [ax, ay, az] = ca[u], [bx, by, bz] = cb[v];
        for (let e = 0; e < 10; e++) {
          const [ex, ey, ez] = MOMENTS[e];
          m[e] += d * X[0][3 * (ax + ex) + bx] * X[1][3 * (ay + ey) + by] * X[2][3 * (az + ez) + bz];
        }
      }
    }
  };

  const sitePositions = [], siteMoments = [];
  // a site at point c with the electrons of m (moments about A; s = A - c), times `share`
  const moved = new Float64Array(10);
  const shifted = (s, share) => {
    const [m0, mx, my, mz] = m, d = [mx, my, mz];
    moved[0] = m0;
    for (let i = 0; i < 3; i++) moved[1 + i] = d[i] + s[i] * m0;
    const at = (i, j) => m[[4, 5, 7, 5, 6, 8, 7, 8, 9][3 * i + j]] + s[i] * d[j] + s[j] * d[i] + s[i] * s[j] * m0;
    moved.set([at(0, 0), at(0, 1), at(1, 1), at(0, 2), at(1, 2), at(2, 2)], 4);
    for (let e = 0; e < 10; e++) moved[e] *= -share; // electrons are negative
    return moved;
  };
  const add = (site, values) => { for (let e = 0; e < 10; e++) siteMoments[site][e] += values[e]; };
  for (let a = 0; a < nat; a++) {
    sitePositions.push(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]);
    siteMoments.push(new Float64Array(10));
    siteMoments[a][0] = VALENCE[numbers[a]];
  }
  const zero = [0, 0, 0];
  for (let a = 0; a < nat; a++) {
    const pa = positions.subarray(3 * a, 3 * a + 3);
    m.fill(0);
    for (const sa of byAtom[a]) for (const sb of byAtom[a]) { const before = m[0]; shellPair(sa, sb, 1); held[sa] += m[0] - before; }
    add(a, shifted(zero, 1));
    for (let b = 0; b < a; b++) {
      const pb = positions.subarray(3 * b, 3 * b + 3), ab = [pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]];
      if (softest[a] * softest[b] / (softest[a] + softest[b]) * (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]) > 32) continue; // nothing is shared that far
      m.fill(0);
      for (const sa of byAtom[a]) for (const sb of byAtom[b]) {
        const before = m[0];
        shellPair(sa, sb, 2);
        held[sa] += 0.5 * (m[0] - before); held[sb] += 0.5 * (m[0] - before);
      }
      if (Math.abs(m[0]) >= bondShare) {
        sitePositions.push(0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2]));
        siteMoments.push(Float64Array.from(shifted(ab.map((v) => 0.5 * v), 1)));
      } else {
        add(a, shifted(zero, 0.5));
        add(b, shifted(ab, 0.5));
      }
    }
  }
  const count = siteMoments.length, multipoles = new Float64Array(10 * count);
  siteMoments.forEach((v, k) => {
    const tr = v[4] + v[6] + v[9]; // second moments xx, xy, yy, xz, yz, zz -> traceless
    multipoles.set([v[0], v[1], v[2], v[3], 1.5 * v[4] - 0.5 * tr, 1.5 * v[5], 1.5 * v[6] - 0.5 * tr, 1.5 * v[7], 1.5 * v[8], 1.5 * v[9] - 0.5 * tr], 10 * k);
  });
  // every atom as a round cloud of the electrons it holds: its part of the potential, and its density
  const pen = new Float32Array(nat * PEN_N), round = new Float64Array(nat * PEN_N);
  shells.forEach((sh, s) => {
    const t = shellTable(sh);
    for (let k = 0; k < PEN_N; k++) { pen[sh.atom * PEN_N + k] += held[s] * t.g[k]; round[sh.atom * PEN_N + k] += held[s] * t.rho[k]; }
  });
  // where an atom's own cloud has thinned to half the density of the surface: neighbours add
  // to it, but not twice as much that far out
  let reach = 0;
  for (let a = 0; a < nat; a++) {
    let k = PEN_N - 1;
    while (k > 0 && round[a * PEN_N + k - 1] < 0.5 * ESP_DENSITY) k--;
    // an atom with no cloud to speak of (a bare K+) still takes the room of a small atom in the picture
    reach = Math.max(reach, Math.hypot(positions[3 * a], positions[3 * a + 1], positions[3 * a + 2]) + Math.max(3, k * PEN_R / (PEN_N - 1)));
  }
  return { n: count, atoms: nat, positions: Float64Array.from(sitePositions), multipoles, pen, reach };
}

// Potential (Eh per unit charge) of sites at point p (bohr): positive where a positive charge
// is pushed away.
export function potential(sites, p) {
  const { positions, multipoles: M, pen } = sites;
  let v = 0;
  for (let k = 0; k < sites.n; k++) {
    const x = p[0] - positions[3 * k], y = p[1] - positions[3 * k + 1], z = p[2] - positions[3 * k + 2], q = 10 * k;
    const r2 = Math.max(x * x + y * y + z * z, ESP_CORE * ESP_CORE), ir = 1 / Math.sqrt(r2), ir3 = ir * ir * ir;
    if (k < sites.atoms && r2 < PEN_R * PEN_R) {
      const t = Math.sqrt(r2) / PEN_R * (PEN_N - 1), i = Math.min(Math.floor(t), PEN_N - 2), at = k * PEN_N + i;
      v += (pen[at] + (pen[at + 1] - pen[at]) * (t - i)) * ir;
    }
    v += M[q] * ir + (M[q + 1] * x + M[q + 2] * y + M[q + 3] * z) * ir3
      + (M[q + 4] * x * x + M[q + 6] * y * y + M[q + 9] * z * z + 2 * (M[q + 5] * x * y + M[q + 7] * x * z + M[q + 8] * y * z)) * ir3 * ir * ir;
  }
  return v;
}

// The potential at which the colours are full when the scale is stretched to a scene: the largest magnitude found on
// spheres of roughly van der Waals size around the atoms (the parts not inside a neighbour's
// sphere), which stand in for the density surface. Never under 0.02 Eh, so that a molecule
// with no polarity to speak of (methane, hydrogen) comes out pale rather than as noise.
const DIRECTIONS = (() => {
  const d = [], n = 24;
  for (let i = 0; i < n; i++) { // a spiral over the sphere
    const z = 1 - (2 * i + 1) / n, s = Math.sqrt(1 - z * z), ph = i * 2.399963;
    d.push([s * Math.cos(ph), s * Math.sin(ph), z]);
  }
  return d;
})();
export function potentialScale(numbers, positions, sites) {
  const nat = numbers.length, radius = numbers.map((z) => (COVALENT[z] + 0.9) / BOHR_ANGSTROM);
  let top = 0;
  for (let a = 0; a < nat; a++) for (const d of DIRECTIONS) {
    const p = [0, 1, 2].map((c) => positions[3 * a + c] + radius[a] * d[c]);
    let outside = true;
    for (let b = 0; b < nat && outside; b++)
      if (b !== a && Math.hypot(p[0] - positions[3 * b], p[1] - positions[3 * b + 1], p[2] - positions[3 * b + 2]) < radius[b]) outside = false;
    if (outside) top = Math.max(top, Math.abs(potential(sites, p)));
  }
  return Math.min(1, Math.max(0.02, top));
}
