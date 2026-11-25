// SMILES in and out. No DOM and no engine, so it runs in Node as well.
//
// parseSMILES reads the text into a graph: atoms with their hydrogens counted and bonds with
// whole orders (aromatic rings get alternating double bonds). fromSMILES hands that graph to
// the editor, which gives it a starting geometry. toSMILES writes an editor molecule back.
//
// What is not kept: the handedness marks (@, @@, / and \) and isotopes are read and dropped.
// The geometry built from the text is one of the mirror forms, whichever comes out.
import { NUMBERS } from './units.js';
import { fromGraph, neighbours, hydrogens, bondLength } from './editor.js';

// the usual numbers of bonds of the elements written without brackets
const NORMAL = { B: [3], C: [4], N: [3, 5], O: [2], P: [3, 5], S: [2, 4, 6], F: [1], Cl: [1], Br: [1], I: [1] };
const AROMATIC = { b: 'B', c: 'C', n: 'N', o: 'O', p: 'P', s: 'S', se: 'Se', as: 'As' };
const ORDER = { '-': 1, '=': 2, '#': 3, ':': 0, '/': 1, '\\': 1 }; // 0: aromatic, decided later

// hydrogens of an unbracketed atom with `used` bonds: up to the nearest usual number
const implicit = (el, used) => { const v = NORMAL[el].find((x) => x >= used); return v === undefined ? 0 : v - used; };
// the number of bonds an aromatic atom has room for, given its charge
function room(el, charge) {
  const z = NUMBERS[el];
  if (el === 'B') return 3 - charge;
  if (el === 'C') return 4 - Math.abs(charge);
  if (z === 7 || z === 15 || z === 33) return 3 + charge;
  return 2 + charge; // O, S, Se
}

class SmilesError extends Error {}
const fail = (text, pos, what) => { throw new SmilesError(`SMILES: ${what} ${pos < text.length ? `at position ${pos + 1} (“${text[pos]}”)` : 'at the end'}.`); };

// text -> { atoms: [{ el, h, charge, aromatic }], bonds: [[i, j, order]], charge }
export function parseSMILES(input) {
  const text = String(input).trim().split(/\s/)[0]; // what follows a space is a name
  if (!text) throw new SmilesError('SMILES: the text is empty.');
  const atoms = [], bonds = [];
  const stack = [], open = new Map(); // branch points; ring numbers waiting for their other end
  let prev = -1, pending = null, pos = 0;
  const element = (sym, at) => { if (!(sym in NUMBERS)) throw new SmilesError(`SMILES: element ${sym} at position ${at + 1} is not among the 36 the program works with.`); return sym; };
  const link = (i, j, symbol, at) => {
    if (i === j || bonds.some(([a, b]) => (a === i && b === j) || (a === j && b === i))) fail(text, at, 'two atoms are bonded twice');
    const both = atoms[i].aromatic && atoms[j].aromatic;
    bonds.push([i, j, symbol === null ? (both ? 0 : 1) : ORDER[symbol]]);
  };
  const add = (atom, at) => {
    atoms.push(atom);
    const k = atoms.length - 1;
    if (prev >= 0) link(prev, k, pending, at);
    else if (pending !== null) fail(text, at, 'a bond with no atom before it');
    prev = k; pending = null;
  };
  while (pos < text.length) {
    const ch = text[pos], start = pos;
    if (ch === '[') {
      const m = /^\[(\d*)(se|as|[A-Z][a-z]?|[bcnops])(@{0,2}(?:[A-Z]{2}\d+)?)(H\d*)?(\+\+|--|[+-]\d*)?(?::\d+)?\]/.exec(text.slice(pos));
      if (!m) fail(text, pos, 'an atom in brackets that cannot be read');
      const aromatic = m[2] in AROMATIC, sign = m[5] ?? '';
      const charge = !sign ? 0 : sign === '++' ? 2 : sign === '--' ? -2 : (sign[0] === '+' ? 1 : -1) * (sign.length > 1 ? +sign.slice(1) : 1);
      add({ el: element(aromatic ? AROMATIC[m[2]] : m[2], pos + 1 + m[1].length), h: m[4] ? (m[4].length > 1 ? +m[4].slice(1) : 1) : 0, charge, aromatic, bracket: true }, pos);
      pos += m[0].length;
    } else if (/[A-Za-z*]/.test(ch)) {
      const two = text.slice(pos, pos + 2);
      const sym = two === 'Cl' || two === 'Br' ? two : ch;
      if (!(sym in NORMAL) && !(sym in AROMATIC)) fail(text, pos, sym === 'H' ? 'a hydrogen outside brackets (write [H])' : 'an element that has to be written in brackets, or an unknown one');
      const aromatic = sym in AROMATIC;
      add({ el: element(aromatic ? AROMATIC[sym] : sym, pos), h: -1, charge: 0, aromatic, bracket: false }, pos);
      pos += sym.length;
    } else if (ch in ORDER) {
      if (pending !== null || prev < 0) fail(text, pos, 'a bond sign out of place');
      pending = ch; pos++;
    } else if (ch === '$') fail(text, pos, 'a quadruple bond, which the editor does not have');
    else if (/\d|%/.test(ch)) {
      const m = /^(\d|%\d\d|%\(\d+\))/.exec(text.slice(pos));
      if (!m || prev < 0) fail(text, pos, 'a ring number out of place');
      const n = +m[0].replace(/\D/g, '');
      if (open.has(n)) {
        const o = open.get(n);
        if (pending !== null && o.symbol !== null && ORDER[pending] !== ORDER[o.symbol]) fail(text, pos, `the two ends of ring bond ${n} disagree`);
        link(o.atom, prev, pending ?? o.symbol, pos);
        open.delete(n);
      } else open.set(n, { atom: prev, symbol: pending, pos });
      pending = null; pos += m[0].length;
    } else if (ch === '(') {
      if (prev < 0 || pending !== null) fail(text, pos, 'a branch with nothing to branch from');
      stack.push(prev); pos++;
    } else if (ch === ')') {
      if (!stack.length || pending !== null) fail(text, pos, 'a closing bracket out of place');
      prev = stack.pop(); pos++;
    } else if (ch === '.') {
      if (pending !== null || prev < 0) fail(text, pos, 'a dot out of place');
      prev = -1; pos++;
    } else fail(text, pos, 'a character that does not belong to SMILES');
    if (pos === start) fail(text, pos, 'something that cannot be read');
  }
  if (pending !== null) fail(text, text.length, 'a bond that leads nowhere');
  if (stack.length) fail(text, text.length, 'a branch that is not closed');
  if (open.size) { const [n, o] = [...open][0]; fail(text, o.pos, `ring bond ${n} that is never closed`); }
  if (!atoms.length) throw new SmilesError('SMILES: no atoms found.');

  kekulize(atoms, bonds);
  const used = atoms.map(() => 0);
  for (const [i, j, order] of bonds) { used[i] += order; used[j] += order; }
  atoms.forEach((a, k) => { if (a.h < 0) a.h = implicit(a.el, used[k]); });
  return { atoms: atoms.map(({ el, h, charge }) => ({ el, h, charge })), bonds, charge: atoms.reduce((s, a) => s + a.charge, 0) };
}

// A perfect matching: every atom that `wants` gets one of its edges ([other atom, bond]),
// no two atoms the same one. Returns the bond of each atom (-1 for the rest), or null.
function matching(wants, edges) {
  const n = wants.length, mate = new Array(n).fill(-1);
  const match = () => {
    // the unmatched atom with the fewest partners left goes first
    let k = -1, fewest = Infinity;
    for (let a = 0; a < n; a++) {
      if (!wants[a] || mate[a] >= 0) continue;
      const free = edges[a].filter(([o]) => mate[o] < 0).length;
      if (free < fewest) { fewest = free; k = a; }
    }
    if (k < 0) return true;
    for (const [o, c] of edges[k]) {
      if (mate[o] >= 0) continue;
      mate[k] = c; mate[o] = c;
      if (match()) return true;
      mate[k] = -1; mate[o] = -1;
    }
    return false;
  };
  return match() ? mate : null;
}

// Aromatic bonds (order 0) become single and double. An aromatic atom takes a double bond
// when it has room for one beside the bonds and hydrogens it already has; the double bonds
// are then a perfect matching of those atoms along the ring bonds.
function kekulize(atoms, bonds) {
  if (!bonds.some((b) => b[2] === 0)) return;
  const n = atoms.length, used = new Array(n).fill(0);
  for (const [i, j, order] of bonds) { used[i] += order || 1; used[j] += order || 1; }
  const wants = atoms.map((a, k) => a.aromatic && room(a.el, a.charge) - used[k] - Math.max(0, a.h) >= 1);
  // a bond written between two aromatic atoms without a sign is aromatic only inside a ring
  const inRing = (b) => {
    const [from, to] = bonds[b], seen = new Set([from]), queue = [from];
    while (queue.length) {
      const k = queue.pop();
      for (let c = 0; c < bonds.length; c++) {
        if (c === b) continue;
        const o = bonds[c][0] === k ? bonds[c][1] : bonds[c][1] === k ? bonds[c][0] : -1;
        if (o < 0 || seen.has(o)) continue;
        if (o === to) return true;
        seen.add(o); queue.push(o);
      }
    }
    return false;
  };
  const edges = atoms.map(() => []);
  bonds.forEach((b, c) => {
    if (b[2] !== 0) return;
    if (wants[b[0]] && wants[b[1]] && inRing(c)) { edges[b[0]].push([b[1], c]); edges[b[1]].push([b[0], c]); }
  });
  const mate = matching(wants, edges);
  if (!mate) throw new SmilesError('SMILES: the double bonds of an aromatic ring cannot be placed. A ring nitrogen that carries a hydrogen is written [nH].');
  const double = new Set(mate.filter((c) => c >= 0));
  bonds.forEach((b, c) => { if (b[2] === 0) b[2] = double.has(c) ? 2 : 1; });
}

// text -> an editor molecule with a starting geometry (to be relaxed by the engine)
export const fromSMILES = (text) => fromGraph(parseSMILES(text));

// Does the text read as SMILES? Used to tell a formula in SMILES from a name.
export function isSMILES(text) {
  try { parseSMILES(text); return true; } catch { return false; }
}

// ---- out

// Where the charge of the molecule sits, atom by atom. The editor keeps one charge for the
// whole molecule, so this is read from the bonds: an atom with one bond too many is positive
// (N with four, O with three), boron with four is negative, and what is left of the charge
// goes to the atoms with a bond too few, the most electronegative first, or to bare metals.
// An atom with one bond too many can also be a radical (NO2 read as N(=O)=O): when the
// charges do not add up the first way, they are tried without that rule (onium = false).
function formalCharges(m, heavy, used, onium = true) {
  const charge = new Map(heavy.map((k) => [k, 0]));
  let left = m.charge;
  const give = (k, q) => { charge.set(k, q); left -= q; };
  for (const k of onium ? heavy : []) {
    const el = m.elements[k], z = NUMBERS[el], v = used[k];
    if ((z === 7 || z === 15) && v === 4) give(k, 1);
    else if ((z === 8 || z === 16) && v === 3) give(k, 1);
    else if ((z === 5 || z === 13) && v === 4) give(k, -1);
  }
  const lacking = (els, bonds) => heavy.filter((k) => charge.get(k) === 0 && els.includes(m.elements[k]) && used[k] === bonds);
  if (left < 0) for (const k of [...lacking(['F', 'Cl', 'Br'], 0), ...lacking(['O', 'S', 'Se'], 1), ...lacking(['N', 'P'], 2), ...lacking(['C', 'Si'], 3), ...lacking(['O', 'S', 'Se'], 0), ...lacking(['H'], 0)]) {
    if (left >= 0) break;
    give(k, m.elements[k] !== 'H' && used[k] === 0 && ['O', 'S', 'Se'].includes(m.elements[k]) && left <= -2 ? -2 : -1);
  }
  if (left > 0) {
    for (const k of heavy) {
      if (left <= 0) break;
      const z = NUMBERS[m.elements[k]];
      if (used[k] !== 0 || charge.get(k) !== 0) continue;
      if ([3, 11, 19].includes(z) || z === 1) give(k, 1);
      else if ((z === 4 || z === 12 || z === 20 || (z > 20 && z < 31)) && left >= 2) give(k, 2);
      else if (z === 13 && left >= 3) give(k, 3);
    }
    for (const k of lacking(['C', 'Si'], 3)) { if (left <= 0) break; give(k, 1); }
  }
  // the ions of a salt, written without their charges: a bare metal and a bare halogen
  if (left === 0) {
    const bare = (zs) => heavy.filter((k) => used[k] === 0 && charge.get(k) === 0 && zs.includes(NUMBERS[m.elements[k]]));
    const metals = bare([3, 11, 19]), halogens = bare([9, 17, 35]);
    for (let c = 0; c < Math.min(metals.length, halogens.length); c++) { give(metals[c], 1); give(halogens[c], -1); }
  }
  if (left !== 0 && onium) return formalCharges(m, heavy, used, false);
  if (left !== 0) throw new Error(`SMILES needs a charge on particular atoms, and the bonds of this structure do not show where its charge of ${m.charge > 0 ? '+' : ''}${m.charge} sits.`);
  return charge;
}

// The bonds of a molecule whose geometry came without them are read from the distances, and
// the bonds of an aromatic ring, all of one length between single and double, are read as
// single. Here they get their double bonds back: where neighbouring atoms are each one bond
// short of the usual number and the bond between them is shorter than a single one, and every
// such atom can be paired off, the pairs are written as double bonds.
function delocalized(m) {
  const n = m.elements.length;
  // the length of each bond against a single one; between 0.885 and 0.95 it may be aromatic
  const ratio = m.bonds.map(([i, j]) => Math.hypot(m.positions[3 * i] - m.positions[3 * j], m.positions[3 * i + 1] - m.positions[3 * j + 1], m.positions[3 * i + 2] - m.positions[3 * j + 2]) / bondLength(m.elements[i], m.elements[j]));
  // a bond read as double that is long for one (C-N in pyridine) joins the bonds to be decided
  const bonds = m.bonds.map((b, c) => (b[2] === 2 && ratio[c] > 0.885 ? [b[0], b[1], 1] : b));
  const count = new Array(n).fill(0);
  for (const [i, j, order] of bonds) { count[i] += order; count[j] += order; }
  const short = m.elements.map((el, k) => el in NORMAL && NORMAL[el][0] - count[k] === 1);
  const edges = m.elements.map(() => []);
  bonds.forEach(([i, j, order], c) => {
    if (order === 1 && short[i] && short[j] && ratio[c] < 0.95) { edges[i].push([j, c]); edges[j].push([i, c]); }
  });
  const mate = matching(edges.map((e) => e.length > 0), edges);
  if (!mate) return m.bonds;
  const double = new Set(mate.filter((c) => c >= 0)), out = bonds.map((b, c) => (double.has(c) ? [b[0], b[1], 2] : b));
  // kept only if it leaves no more atoms short of bonds than there were (in nitrate, all of
  // whose bonds are of that middle length, it would)
  const lacking = (list) => {
    const sum = new Array(n).fill(0);
    for (const [i, j, order] of list) { sum[i] += order; sum[j] += order; }
    return m.elements.filter((el, k) => el in NORMAL && sum[k] < NORMAL[el][0]).length;
  };
  return lacking(out) <= lacking(m.bonds) ? out : m.bonds;
}

// an editor molecule -> SMILES with whole bond orders (rings are not written as aromatic)
export function toSMILES(m) {
  const n = m.elements.length;
  // hydrogens with one bond to another atom are written into that atom
  const folded = new Set(), count = new Array(n).fill(0);
  for (let k = 0; k < n; k++) for (const h of hydrogens(m, k)) { folded.add(h); count[k]++; }
  const heavy = m.elements.map((_, k) => k).filter((k) => !folded.has(k));
  m = { ...m, bonds: delocalized(m) };
  const used = new Array(n).fill(0);
  for (const [i, j, order] of m.bonds) { used[i] += order; used[j] += order; }
  // An anion like acetate or nitrate has its C-O or N-O bonds all alike, and read from the
  // distances they are all double: one bond too many for the central atom. For each unit of
  // negative charge one such bond to an end oxygen becomes single, and that oxygen takes
  // the charge.
  let spare = -m.charge;
  m.bonds = m.bonds.map((b) => {
    const [i, j, order] = b, most = (k) => NORMAL[m.elements[k]]?.at(-1) ?? Infinity, end = (k) => ['O', 'S'].includes(m.elements[k]) && used[k] === 2;
    if (spare <= 0 || order !== 2 || !((used[i] > most(i) && end(j)) || (used[j] > most(j) && end(i)))) return b;
    spare--; used[i]--; used[j]--;
    return [i, j, 1];
  });
  const charge = formalCharges(m, heavy, used);
  const nb = (k) => neighbours(m, k).filter((x) => !folded.has(x.atom));

  // one pass to find the tree and the bonds that close rings, a second one to write
  const seen = new Set(), children = new Map(), rings = new Map(heavy.map((k) => [k, []]));
  let ringCount = 0;
  const walk = (k, from) => {
    seen.add(k); children.set(k, []);
    for (const { atom, order } of nb(k)) {
      if (atom === from) continue;
      if (!seen.has(atom)) { children.get(k).push({ atom, order }); walk(atom, k); }
      else if (!rings.get(k).some((r) => r.other === atom)) {
        const number = ++ringCount;
        rings.get(atom).push({ number, order, other: k, opens: true });
        rings.get(k).push({ number, order, other: atom, opens: false });
      }
    }
  };
  const sign = (order) => (order === 2 ? '=' : order === 3 ? '#' : '');
  const atomText = (k) => {
    const el = m.elements[k], q = charge.get(k), h = count[k];
    if (q === 0 && el in NORMAL && implicit(el, used[k] - h) === h) return el;
    return `[${el}${h ? (h > 1 ? `H${h}` : 'H') : ''}${q ? (q > 0 ? '+' : '-') + (Math.abs(q) > 1 ? Math.abs(q) : '') : ''}]`;
  };
  const write = (k) => {
    let s = atomText(k);
    for (const r of rings.get(k).sort((a, b) => a.number - b.number)) s += (r.opens ? sign(r.order) : '') + (r.number > 9 ? `%${r.number}` : r.number);
    const kids = children.get(k);
    kids.forEach((c, i) => { const part = sign(c.order) + write(c.atom); s += i < kids.length - 1 ? `(${part})` : part; });
    return s;
  };
  const parts = [];
  for (const k of heavy) if (!seen.has(k)) { walk(k, -1); parts.push(write(k)); }
  return parts.join('.');
}
