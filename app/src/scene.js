// A scene is one or two molecules of the library put side by side: elements, positions
// (bohr, centred on the centroid) and, for every atom, the fragment it came from.
import { NUMBERS, VALENCE, COVALENT } from './units.js';

// bohr between the facing atoms of two fragments when the second is added. Measured for
// BH3 + NH3 at 300 K: from 6 bohr they bind on their own within 0.1–3 ps, from 8 bohr they
// did not in 4 ps (6 runs), so bringing them together is left to the hand.
const GAP = 8.0;

function centred(bohr) {
  const n = bohr.length / 3, c = [0, 0, 0];
  for (let i = 0; i < bohr.length; i++) c[i % 3] += bohr[i] / n;
  // already centred: leave the numbers untouched (a shift in the last digit is enough to
  // rotate degenerate orbitals among themselves)
  if (Math.hypot(...c) < 1e-9) return Float64Array.from(bohr);
  return Float64Array.from(bohr, (v, i) => v - c[i % 3]);
}

// number of valence electrons, and the unpaired electrons a molecule has on its own:
// what the library says (O2: 2), otherwise one if the count is odd
const electronsOf = (m) => m.elements.reduce((s, el) => s + VALENCE[NUMBERS[el]], 0) - (m.charge || 0);
const unpairedOf = (m) => m.unpaired ?? electronsOf(m) % 2;

// first, second: library entries with .bohr; second may be null.
// scene.unpaired is the number of unpaired electrons the scene starts with: a molecule's
// own; for two molecules with unpaired electrons the difference, i.e. their spins opposite,
// which is the arrangement that lets them form a bond.
export function buildScene(first, second) {
  if (!second) {
    return {
      key: first.key, name: first.name, note: first.note, charge: first.charge || 0, elements: first.elements,
      numbers: first.elements.map((el) => NUMBERS[el]), positions: centred(first.bohr), fragments: new Uint8Array(first.elements.length), nFragments: 1,
      electrons: electronsOf(first), unpaired: unpairedOf(first),
    };
  }
  const u1 = unpairedOf(first), u2 = unpairedOf(second);
  // the second molecule goes to +x of the first, as in the library orientation: for
  // BH3 + NH3 that points the nitrogen lone pair at the empty orbital of boron
  const a = centred(first.bohr), b = centred(second.bohr);
  let right = -Infinity, left = Infinity;
  for (let i = 0; i < a.length; i += 3) right = Math.max(right, a[i]);
  for (let i = 0; i < b.length; i += 3) left = Math.min(left, b[i]);
  const shift = right - left + GAP, all = new Float64Array(a.length + b.length);
  all.set(a);
  for (let i = 0; i < b.length; i++) all[a.length + i] = b[i] + (i % 3 === 0 ? shift : 0);
  const elements = [...first.elements, ...second.elements];
  const fragments = new Uint8Array(elements.length).fill(1, first.elements.length);
  return {
    key: `${first.key}+${second.key}`, name: `${first.name} + ${second.name}`, note: '', charge: (first.charge || 0) + (second.charge || 0),
    elements, numbers: elements.map((el) => NUMBERS[el]), positions: centred(all), fragments, nFragments: 2,
    electrons: electronsOf(first) + electronsOf(second), unpaired: u1 && u2 ? Math.abs(u1 - u2) : u1 + u2,
  };
}

// Which atoms currently hang together. Two atoms count as bonded when they are closer than
// 1.3 times the sum of their covalent radii (and never less than 1.05 A: a vibrating H2
// stretches past 0.81 A, which the radii alone would call broken); a piece is everything connected through such
// bonds. Returns { piece: for every atom the number of its piece, count }.
// (The fragments above record where an atom came from and never change; pieces follow what
// has happened since: after a molecule falls apart, each part is its own piece.)
export function pieces(scene) {
  const n = scene.elements.length, p = scene.positions, piece = new Int32Array(n).fill(-1);
  let count = 0;
  for (let start = 0; start < n; start++) {
    if (piece[start] >= 0) continue;
    const stack = [start];
    piece[start] = count;
    while (stack.length) {
      const i = stack.pop();
      for (let j = 0; j < n; j++) {
        if (piece[j] >= 0) continue;
        const d = Math.hypot(p[3 * i] - p[3 * j], p[3 * i + 1] - p[3 * j + 1], p[3 * i + 2] - p[3 * j + 2]) * 0.529177210903;
        if (d < Math.max(1.05, 1.3 * (COVALENT[NUMBERS[scene.elements[i]]] + COVALENT[NUMBERS[scene.elements[j]]]))) { piece[j] = count; stack.push(j); }
      }
    }
    count++;
  }
  return { piece, count };
}

// What the scene consists of right now, as formulas of its pieces: 'H2O + NH4'.
// Carbon first, then hydrogen, then the rest alphabetically; equal pieces are counted.
export function composition(scene) {
  const { piece, count } = pieces(scene), names = new Map();
  for (let p = 0; p < count; p++) {
    const n = {};
    scene.elements.forEach((el, k) => { if (piece[k] === p) n[el] = (n[el] || 0) + 1; });
    const order = Object.keys(n).sort((a, b) => (a === 'C' ? -2 : a === 'H' && n.C ? -1 : 0) - (b === 'C' ? -2 : b === 'H' && n.C ? -1 : 0) || a.localeCompare(b));
    const f = order.map((el) => el + (n[el] > 1 ? n[el] : '')).join('');
    names.set(f, (names.get(f) || 0) + 1);
  }
  return [...names].map(([f, k]) => (k > 1 ? `${k} ${f}` : f)).join(' + ');
}

// positions with every atom for which member(k) is true moved by `delta`
export function translated(scene, member, delta) {
  const p = Float64Array.from(scene.positions);
  for (let k = 0; k < scene.elements.length; k++) if (member(k)) for (let c = 0; c < 3; c++) p[3 * k + c] += delta[c];
  return p;
}

// positions with the atoms for which member(k) is true rotated about their own centroid:
// by `a` about unit axis u, then by `b` about unit axis w
export function rotated(scene, member, u, a, w, b) {
  const p = Float64Array.from(scene.positions), c = [0, 0, 0];
  let n = 0;
  for (let k = 0; k < scene.elements.length; k++) if (member(k)) { n++; for (let d = 0; d < 3; d++) c[d] += p[3 * k + d]; }
  const rot = (v, axis, angle) => {
    const cs = Math.cos(angle), sn = Math.sin(angle), dt = axis[0] * v[0] + axis[1] * v[1] + axis[2] * v[2];
    const cr = [axis[1] * v[2] - axis[2] * v[1], axis[2] * v[0] - axis[0] * v[2], axis[0] * v[1] - axis[1] * v[0]];
    return v.map((x, i) => x * cs + cr[i] * sn + axis[i] * dt * (1 - cs));
  };
  for (let k = 0; k < scene.elements.length; k++) if (member(k)) {
    const v = rot(rot([0, 1, 2].map((d) => p[3 * k + d] - c[d] / n), u, a), w, b);
    for (let d = 0; d < 3; d++) p[3 * k + d] = c[d] / n + v[d];
  }
  return p;
}
