// The molecule editor's model: atoms, bonds and the operations on them. No DOM and no engine,
// so it runs in Node as well.
//
// A molecule here is { elements: ['C', ...], positions: [x, y, z, ...] (bohr), bonds:
// [[i, j, order], ...], charge }. Every operation takes a molecule and returns a new one, or
// throws an Error whose message is meant for the user.
//
// Bonds and their orders exist only here. The engine knows nuclei and a number of electrons;
// the editor keeps bonds to count hydrogens and to put new atoms in a sensible place. What
// it produces is a starting geometry: the real one comes from relaxing it with the engine.
//
// Hydrogens are counted by difference: an operation adds or takes as many as it changes the
// number of bonds by. Nothing here decides how many hydrogens an atom "should" have, so a
// molecule that came from elsewhere (the library, an XYZ file) can be edited from any point.
import { NUMBERS, VALENCE, COVALENT, BOHR_ANGSTROM } from './units.js';

// the usual number of bonds of an element: how many hydrogens a new atom arrives with.
// Elements that are not listed arrive bare, and nothing is added to or taken from them.
export const BONDS = { H: 1, B: 3, C: 4, N: 3, O: 2, F: 1, Si: 4, P: 3, S: 2, Cl: 1, Br: 1, I: 1 };
export const MAX_ATOMS = 10000;

const SHORTER = [0, 1, 0.87, 0.78]; // a double and a triple bond against a single one
const radius = (el) => COVALENT[NUMBERS[el]] / BOHR_ANGSTROM;
export const bondLength = (a, b, order = 1) => (radius(a) + radius(b)) * SHORTER[order];

const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
const unit = (u) => { const l = Math.hypot(u[0], u[1], u[2]) || 1; return [u[0] / l, u[1] / l, u[2] / l]; };
const at = (m, k) => [m.positions[3 * k], m.positions[3 * k + 1], m.positions[3 * k + 2]];
const distance = (m, i, j) => Math.hypot(...sub(at(m, i), at(m, j)));

export const emptyMolecule = () => ({ elements: [], positions: [], bonds: [], charge: 0 });

// ---- a working copy: atoms are struck out rather than removed, so that indices hold until
// the operation is over
const draft = (m) => ({ elements: [...m.elements], positions: [...m.positions], bonds: m.bonds.map((b) => [...b]), charge: m.charge, dead: new Set() });
function finish(d) {
  if (d.elements.length - d.dead.size > MAX_ATOMS) throw new Error(`The editor takes molecules of up to ${MAX_ATOMS} atoms.`);
  const index = [];
  let n = 0;
  d.elements.forEach((_, k) => index.push(d.dead.has(k) ? -1 : n++));
  return {
    elements: d.elements.filter((_, k) => index[k] >= 0),
    positions: d.positions.filter((_, i) => index[(i / 3) | 0] >= 0),
    bonds: d.bonds.map(([i, j, o]) => [index[i], index[j], o]),
    charge: d.charge,
  };
}
const addAtom = (d, el, p) => { d.elements.push(el); d.positions.push(...p); return d.elements.length - 1; };
function strike(d, k) {
  d.dead.add(k);
  d.bonds = d.bonds.filter(([i, j]) => i !== k && j !== k);
}

// the atoms bonded to k: [{ atom, order }]
export function neighbours(m, k) {
  const out = [];
  for (const [i, j, order] of m.bonds) {
    if (i === k) out.push({ atom: j, order });
    else if (j === k) out.push({ atom: i, order });
  }
  return out;
}
// the hydrogens that belong to atom k: bonded to it and to nothing else
export function hydrogens(m, k) {
  if (m.elements[k] === 'H') return [];
  return neighbours(m, k).map((n) => n.atom).filter((a) => m.elements[a] === 'H' && neighbours(m, a).length === 1);
}
const findBond = (m, i, j) => m.bonds.findIndex(([a, b]) => (a === i && b === j) || (a === j && b === i));

// lone pairs of atom k once it has `used` bonds (orders summed): only for main-group elements
function lonePairs(el, used) {
  const z = NUMBERS[el];
  if ((z > 20 && z < 31) || (z > 38 && z < 49) || (z > 56 && z < 81)) return 0; // the d metals and the lanthanides
  return Math.max(0, Math.floor((VALENCE[z] - used) / 2));
}

// Directions for `count` new atoms on atom k. The bonds it already has stay where they are;
// the new ones and the lone pairs spread over the sphere as far from them and from each
// other as they can (points repelling on a sphere: a tetrahedron for four, a triangle for
// three, a line for two). With a single bond already there, the whole set can still turn
// about it: it is turned to where the new atoms are farthest from the rest of the molecule
// (staggered against the neighbour's own substituents), or, on a double bond, into the plane
// of the other end.
// `used` is the number of bonds the atom will have in the end, when not all of them are in
// the molecule yet (a molecule being built from its graph).
function spread(d, k, count, length, used) {
  const centre = at(d, k), nb = neighbours(d, k);
  const fixed = nb.map((n) => unit(sub(at(d, n.atom), centre)));
  const lone = lonePairs(d.elements[k], used ?? nb.reduce((s, n) => s + n.order, 0) + count);
  // a deterministic start, well spread and in no special position
  let free = Array.from({ length: count + lone }, (_, i) => {
    const z = 1 - 2 * (i + 0.5) / (count + lone), r = Math.sqrt(1 - z * z), phi = 2.399963 * i + 0.7;
    return unit([r * Math.cos(phi) + 0.13, r * Math.sin(phi) - 0.07, z + 0.11]);
  });
  // They start on the side away from the bonds that are there. Between three bonds the
  // repulsion is almost level, and a point that started there could stay there: a hydrogen
  // pointing into a cage.
  const taken = fixed.reduce((s, v) => s.map((x, c) => x + v[c]), [0, 0, 0]);
  if (Math.hypot(...taken) > 0.3) free = free.map((v) => unit(sub(v, taken)));
  for (let iter = 0; iter < 400; iter++) {
    free.forEach((v, i) => {
      const f = [0, 0, 0];
      const push = (from) => { const r = sub(v, from), l = Math.hypot(...r) || 1e-3, w = 1 / (l * l * l); for (let c = 0; c < 3; c++) f[c] += w * r[c]; };
      fixed.forEach(push);
      free.forEach((u, j) => { if (j !== i) push(u); });
      const along = dot(f, v), t = f.map((x, c) => x - along * v[c]), size = Math.max(1, Math.hypot(...t));
      free[i] = unit(v.map((x, c) => x + 0.12 * t[c] / size));
    });
  }
  if (fixed.length === 1 && count > 0) {
    const axis = fixed[0], partner = nb[0].atom;
    // the plane a double-bonded neighbour sets with its own substituents
    let normal = null;
    if (nb[0].order === 2 && free.length === 2) {
      for (const other of neighbours(d, partner)) {
        if (other.atom === k) continue;
        const c = cross(axis, unit(sub(at(d, other.atom), at(d, partner))));
        if (Math.hypot(...c) > 0.2) { normal = unit(c); break; }
      }
    }
    const others = d.elements.map((_, a) => a).filter((a) => !d.dead.has(a) && a !== k && a !== partner).map((a) => at(d, a));
    const turned = (angle) => free.map((v) => {
      const cs = Math.cos(angle), sn = Math.sin(angle), cr = cross(axis, v), al = dot(axis, v);
      return v.map((x, c) => x * cs + cr[c] * sn + axis[c] * al * (1 - cs));
    });
    let best = free, bestScore = Infinity;
    for (let step = 0; step < 72; step++) {
      const set = turned(step * Math.PI / 36);
      let score = 0;
      for (const v of set.slice(0, count)) {
        if (normal) score += 100 * dot(v, normal) ** 2;
        for (const q of others) { const r = sub(v.map((x, c) => centre[c] + length * x), q); score += 1 / (dot(r, r) ** 2 + 1e-6); }
      }
      if (score < bestScore - 1e-12) { bestScore = score; best = set; }
    }
    free = best;
  }
  return free.slice(0, count);
}

// gives atom k exactly n hydrogens, all placed afresh
function rehydrogenate(d, k, n) {
  for (const h of hydrogens(d, k)) strike(d, h);
  const length = bondLength(d.elements[k], 'H'), centre = at(d, k);
  for (const v of spread(d, k, n, length)) d.bonds.push([k, addAtom(d, 'H', v.map((x, c) => centre[c] + length * x)), 1]);
}
// changes the number of hydrogens of atom k by `change`, for the elements that keep count
function adjustHydrogens(d, k, change) {
  if (!(d.elements[k] in BONDS) || d.elements[k] === 'H') return;
  rehydrogenate(d, k, Math.max(0, hydrogens(d, k).length + change));
}
// Can atom k take one more bond? Yes if it has a hydrogen to give up for it, or a free
// valence: fewer bonds than is usual for its element (a hydrogen was erased from it).
const hasFreeValence = (m, k) => neighbours(m, k).reduce((s, n) => s + n.order, 0) < BONDS[m.elements[k]];
const canBond = (m, k) => !(m.elements[k] in BONDS) || hydrogens(m, k).length > 0 || hasFreeValence(m, k);

// ---- a rough geometry from the bond graph alone: springs at the ideal bond lengths, at the
// distances the ideal angles give, and a soft push between atoms that are not neighbours.
// Used when a new bond joins atoms that are far apart (closing a ring).
function tidy(d) {
  const alive = d.elements.map((_, k) => k).filter((k) => !d.dead.has(k));
  const springs = new Map(); // 'i,j' -> length
  for (const [i, j, order] of d.bonds) springs.set(`${Math.min(i, j)},${Math.max(i, j)}`, bondLength(d.elements[i], d.elements[j], order));
  const angles = [];
  for (const k of alive) {
    const nb = neighbours(d, k), points = nb.length + lonePairs(d.elements[k], nb.reduce((s, n) => s + n.order, 0));
    const cos = points <= 2 ? -1 : points === 3 ? -0.5 : -1 / 3;
    for (let a = 0; a < nb.length; a++) for (let b = a + 1; b < nb.length; b++) {
      const i = Math.min(nb[a].atom, nb[b].atom), j = Math.max(nb[a].atom, nb[b].atom);
      if (springs.has(`${i},${j}`)) continue;
      const la = bondLength(d.elements[k], d.elements[nb[a].atom], nb[a].order), lb = bondLength(d.elements[k], d.elements[nb[b].atom], nb[b].order);
      angles.push([i, j, Math.sqrt(la * la + lb * lb - 2 * la * lb * cos)]);
    }
  }
  const held = new Set([...springs.keys(), ...angles.map(([i, j]) => `${i},${j}`)]);
  const p = d.positions;
  for (let iter = 0; iter < 600; iter++) {
    const move = new Float64Array(p.length);
    const pull = (i, j, want, k, onlyApart) => {
      const r = [p[3 * j] - p[3 * i], p[3 * j + 1] - p[3 * i + 1], p[3 * j + 2] - p[3 * i + 2]], l = Math.hypot(...r) || 1e-3;
      if (onlyApart && l >= want) return;
      const s = k * (l - want) / l / 2;
      for (let c = 0; c < 3; c++) { move[3 * i + c] += s * r[c]; move[3 * j + c] -= s * r[c]; }
    };
    for (const [key, want] of springs) { const [i, j] = key.split(',').map(Number); pull(i, j, want, 0.5); }
    for (const [i, j, want] of angles) pull(i, j, want, 0.25);
    for (let a = 0; a < alive.length; a++) for (let b = a + 1; b < alive.length; b++) {
      const i = alive[a], j = alive[b];
      if (!held.has(`${i},${j}`)) pull(i, j, radius(d.elements[i]) + radius(d.elements[j]) + 1.5 / BOHR_ANGSTROM, 0.2, true);
    }
    for (let i = 0; i < p.length; i++) p[i] += move[i];
  }
}

// how far the worst bond is from its ideal length, as a fraction of that length
const strain = (d) => Math.max(...d.bonds.map(([i, j, order]) => { const want = bondLength(d.elements[i], d.elements[j], order); return Math.abs(distance(d, i, j) - want) / want; }));

// are two atoms that are not bonded nearly as close as bonded ones would be?
function crowded(d) {
  const n = d.elements.length;
  for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++)
    if (distance(d, i, j) < radius(d.elements[i]) + radius(d.elements[j]) + 0.3 / BOHR_ANGSTROM && findBond(d, i, j) < 0) return true;
  return false;
}

// ---- operations

// a new atom at `point`, with its hydrogens
export function placeAtom(m, el, point) {
  const d = draft(m), k = addAtom(d, el, point);
  rehydrogenate(d, k, el === 'H' ? 0 : BONDS[el] ?? 0);
  return finish(d);
}

// hydrogen h becomes an atom of `el`, at bonding distance from what h was bonded to, and
// receives hydrogens of its own
export function replaceHydrogen(m, h, el) {
  const nb = neighbours(m, h);
  if (m.elements[h] !== 'H' || nb.length !== 1) throw new Error('Only a hydrogen with one bond can be replaced this way.');
  const d = draft(m), parent = nb[0].atom, centre = at(d, parent), v = unit(sub(at(d, h), centre));
  const length = bondLength(d.elements[parent], el);
  d.elements[h] = el;
  for (let c = 0; c < 3; c++) d.positions[3 * h + c] = centre[c] + length * v[c];
  rehydrogenate(d, h, Math.max(0, (BONDS[el] ?? 1) - 1));
  return finish(d);
}

// atom k becomes an atom of `el`; its hydrogens change by the difference in the usual
// number of bonds
export function setElement(m, k, el) {
  const d = draft(m), old = d.elements[k];
  if (old === el) return m;
  if (neighbours(m, k).length === 0) { d.elements[k] = el; rehydrogenate(d, k, el === 'H' ? 0 : BONDS[el] ?? 0); return finish(d); }
  if (el === 'H') throw new Error('To put a hydrogen in place of an atom, erase the atom.');
  const count = hydrogens(d, k).length;
  for (const h of hydrogens(d, k)) strike(d, h);
  d.elements[k] = el;
  const change = old in BONDS && el in BONDS && old !== 'H' ? BONDS[el] - BONDS[old] : 0;
  rehydrogenate(d, k, el in BONDS ? Math.max(0, count + change) : 0);
  return finish(d);
}

// a bond between i and j, or one more order on the bond they have; each end gives up a
// hydrogen, or uses a free valence if it has no hydrogen
export function connect(m, i, j) {
  if (i === j) return m;
  if (m.elements[i] === 'H' || m.elements[j] === 'H') throw new Error('A bond is drawn between two atoms that are not hydrogens: click the hydrogen to replace it instead.');
  const b = findBond(m, i, j);
  if (b >= 0 && m.bonds[b][2] >= 3) throw new Error('This bond is already triple.');
  for (const k of [i, j]) if (!canBond(m, k)) throw new Error(`${m.elements[k]} has all its bonds already, and no hydrogen to give up for another one.`);
  const d = draft(m);
  // a free valence is used first; only an atom without one gives up a hydrogen
  const pays = [i, j].filter((k) => !hasFreeValence(m, k));
  if (b >= 0) d.bonds[b][2]++; else d.bonds.push([i, j, 1]);
  for (const k of [i, j]) adjustHydrogens(d, k, pays.includes(k) ? -1 : 0);
  if (b < 0 && distance(d, i, j) > 1.3 * bondLength(d.elements[i], d.elements[j])) {
    tidy(d);
    // the hydrogens again, now that the atoms they sit on have moved
    for (const k of new Set(d.bonds.flatMap(([a, c]) => [a, c]))) if (d.elements[k] !== 'H' && hydrogens(d, k).length) rehydrogenate(d, k, hydrogens(d, k).length);
  }
  return finish(d);
}

// one more hydrogen on atom k (the way to NH4+ and back from a radical)
export function addHydrogen(m, k) {
  if (m.elements[k] === 'H') return m;
  const d = draft(m);
  rehydrogenate(d, k, hydrogens(d, k).length + 1);
  return finish(d);
}

// bond number b: single -> double -> triple -> single
export function cycleBond(m, b) {
  const [i, j, order] = m.bonds[b];
  if (m.elements[i] === 'H' || m.elements[j] === 'H') throw new Error('A bond to a hydrogen is always single.');
  if (order < 3 && canBond(m, i) && canBond(m, j)) return connect(m, i, j);
  if (order === 1) throw new Error('For a double bond each of the two atoms needs a hydrogen to give up or a free valence, and one of them has neither.');
  const d = draft(m);
  d.bonds[b][2] = 1;
  adjustHydrogens(d, i, order - 1);
  adjustHydrogens(d, j, order - 1);
  return finish(d);
}

// removes bond number b; each end receives hydrogens in its place
export function removeBond(m, b) {
  const [i, j, order] = m.bonds[b];
  if (m.elements[i] === 'H' || m.elements[j] === 'H') throw new Error('To remove a hydrogen, click the hydrogen itself.');
  const d = draft(m);
  d.bonds.splice(b, 1);
  adjustHydrogens(d, i, order);
  adjustHydrogens(d, j, order);
  return finish(d);
}

// removes atom k. A hydrogen simply goes (that is how a radical is made); any other atom
// goes with its hydrogens, and each neighbour receives hydrogens in place of the lost bond.
export function erase(m, k) {
  const d = draft(m);
  if (m.elements[k] === 'H') strike(d, k);
  else {
    const mine = hydrogens(d, k), rest = neighbours(d, k).filter((n) => !mine.includes(n.atom));
    for (const h of mine) strike(d, h);
    strike(d, k);
    for (const n of rest) adjustHydrogens(d, n.atom, n.order);
  }
  if (d.dead.size === d.elements.length) throw new Error('The last atom stays: change its element instead.');
  return finish(d);
}

// atom k moved by `delta` together with its hydrogens
export function moveAtom(m, k, delta) {
  const positions = [...m.positions];
  for (const a of [k, ...hydrogens(m, k)]) for (let c = 0; c < 3; c++) positions[3 * a + c] += delta[c];
  return { ...m, positions };
}

// the same molecule with its centroid at the origin
export function centred(m) {
  const n = m.elements.length, c = [0, 1, 2].map((k) => m.positions.reduce((s, v, i) => s + (i % 3 === k ? v : 0), 0) / n);
  if (Math.hypot(...c) < 1e-9) return m;
  return { ...m, positions: m.positions.map((v, i) => v - c[i % 3]) };
}

// ---- in and out

// A molecule from its bond graph alone, as SMILES gives it: { atoms: [{ el, h }], bonds:
// [[i, j, order]], charge }, h being the number of hydrogens of an atom. Each piece of the
// graph is grown from its first atom the way the editor grows a molecule: an atom places
// all its neighbours and hydrogens at once, spread over the sphere. The bonds that close
// rings are added at the end, and the springs of tidy() pull the rings shut. Separate pieces
// (the ions of a salt) are set side by side.
export function fromGraph(graph) {
  const n = graph.atoms.length, total = graph.atoms.reduce((s, a) => s + 1 + a.h, 0);
  if (total > MAX_ATOMS) throw new Error(`The editor takes molecules of up to ${MAX_ATOMS} atoms; this one has ${total}.`);
  const links = graph.atoms.map(() => []);
  for (const [i, j, order] of graph.bonds) { links[i].push({ atom: j, order }); links[j].push({ atom: i, order }); }
  const used = links.map((l, k) => l.reduce((s, x) => s + x.order, 0) + graph.atoms[k].h);
  const placed = new Array(n).fill(-1), pieces = [];
  for (let root = 0; root < n; root++) {
    if (placed[root] >= 0) continue;
    const d = draft(emptyMolecule()), queue = [root], closing = [];
    placed[root] = addAtom(d, graph.atoms[root].el, [0, 0, 0]);
    for (let q = 0; q < queue.length; q++) {
      const k = queue[q], me = placed[k], { el, h } = graph.atoms[k], centre = at(d, me);
      const kids = links[k].filter((x) => placed[x.atom] < 0);
      for (const x of links[k]) if (placed[x.atom] >= 0 && x.atom < k && findBond(d, me, placed[x.atom]) < 0) closing.push([placed[x.atom], me, x.order]);
      // The heavy neighbours take the directions with the most room around them and the
      // hydrogens the rest: a chain then goes on straight (anti) instead of coiling into itself.
      const reach = 1.5 / BOHR_ANGSTROM, others = d.elements.map((_, a) => a).filter((a) => a !== me).map((a) => at(d, a));
      const crowding = (v) => others.reduce((sum, q) => { const r = sub(v.map((x, c) => centre[c] + reach * x), q); return sum + 1 / (dot(r, r) ** 2 + 1e-6); }, 0);
      const dirs = spread(d, me, kids.length + h, reach, used[k]).map((v) => [crowding(v), v]).sort((a, b) => a[0] - b[0]).map((x) => x[1]);
      kids.forEach((x, i) => {
        const length = bondLength(el, graph.atoms[x.atom].el, x.order);
        placed[x.atom] = addAtom(d, graph.atoms[x.atom].el, dirs[i].map((v, c) => centre[c] + length * v));
        d.bonds.push([me, placed[x.atom], x.order]);
        queue.push(x.atom);
      });
      for (const v of dirs.slice(kids.length)) d.bonds.push([me, addAtom(d, 'H', v.map((x, c) => centre[c] + bondLength(el, 'H') * x)), 1]);
    }
    // a bond between two atoms of the tree that is not in it closes a ring
    const rings = closing.filter(([i, j]) => findBond(d, i, j) < 0);
    // the springs are also what parts two branches of a large molecule that grew into each other
    if (rings.length || crowded(d)) {
      d.bonds.push(...rings);
      tidy(d);
      // A cage can come out tangled, with bonds left stretched. Then the atoms are shaken and
      // the springs run again, a few times, and the least strained result is kept.
      let best = [...d.positions], least = strain(d), seed = 12345;
      for (let attempt = 1; attempt <= 12 && least > 0.08; attempt++) {
        d.positions = best.map((v) => { seed = (seed * 1103515245 + 12345) % 2147483648; return v + (seed / 2147483648 - 0.5) * (1 + attempt) ; });
        tidy(d);
        const now = strain(d);
        if (now < least) { least = now; best = [...d.positions]; }
      }
      d.positions = best;
      for (const k of queue) if (graph.atoms[k].h && graph.atoms[k].el !== 'H') rehydrogenate(d, placed[k], graph.atoms[k].h);
    }
    pieces.push(centred(finish(d)));
  }
  // the pieces in a row, 2.5 A between the spheres that hold them
  const out = emptyMolecule();
  let edge = 0;
  pieces.forEach((p, c) => {
    const size = Math.sqrt(Math.max(...p.elements.map((_, k) => dot(at(p, k), at(p, k))))), x = c ? edge + 2.5 / BOHR_ANGSTROM + size : 0, first = out.elements.length;
    out.elements.push(...p.elements);
    out.positions.push(...p.positions.map((v, i) => (i % 3 ? v : v + x)));
    out.bonds.push(...p.bonds.map(([i, j, o]) => [i + first, j + first, o]));
    edge = x + size;
  });
  out.charge = graph.charge ?? 0;
  return centred(out);
}

// Bonds of a geometry that came without them: two atoms closer than 1.3 times the sum of
// their covalent radii are bonded (the criterion of pieces() in scene.js). Between C, N, O,
// P and S a markedly short bond is taken for a double or a triple one.
export function withBonds(elements, positions, charge = 0) {
  const m = { elements: [...elements], positions: Array.from(positions), bonds: [], charge };
  const multiple = new Set(['C', 'N', 'O', 'P', 'S']);
  for (let i = 0; i < elements.length; i++) for (let j = i + 1; j < elements.length; j++) {
    const single = bondLength(elements[i], elements[j]), l = distance(m, i, j);
    if (l >= Math.max(1.05 / BOHR_ANGSTROM, 1.3 * single)) continue;
    const ratio = l / single, both = multiple.has(elements[i]) && multiple.has(elements[j]);
    m.bonds.push([i, j, both && ratio < 0.8 ? 3 : both && ratio < 0.905 ? 2 : 1]);
  }
  return m;
}

export function toXYZ(m, comment = '') {
  const rows = m.elements.map((el, k) => `${el.padEnd(2)} ${at(m, k).map((v) => (v * BOHR_ANGSTROM).toFixed(6).padStart(11)).join(' ')}`);
  return `${m.elements.length}\n${comment}\n${rows.join('\n')}\n`;
}

// XYZ text (angstrom), with or without the two header lines: every line that reads as
// "element x y z" is an atom
export function fromXYZ(text, charge = 0) {
  const elements = [], positions = [];
  for (const line of text.split('\n')) {
    const f = line.trim().split(/[\s,;]+/);
    if (f.length < 4 || !/^[A-Za-z]{1,2}$/.test(f[0])) continue;
    const xyz = f.slice(1, 4).map(Number);
    if (xyz.some((v) => !Number.isFinite(v))) continue;
    const el = f[0][0].toUpperCase() + f[0].slice(1).toLowerCase();
    if (!(el in NUMBERS)) throw new Error(`Element ${f[0]} is not among the 86 (H to Rn) the program works with.`);
    elements.push(el);
    positions.push(...xyz.map((v) => v / BOHR_ANGSTROM));
  }
  if (!elements.length) throw new Error('No atoms found. Each line should read: element x y z, in angstrom.');
  if (elements.length > MAX_ATOMS) throw new Error(`The editor takes molecules of up to ${MAX_ATOMS} atoms; this one has ${elements.length}.`);
  return withBonds(elements, positions, charge);
}
