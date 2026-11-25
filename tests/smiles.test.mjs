// SMILES in and out, the geometry built from it, and what comes from PubChem:  node tests/smiles.test.mjs
// (PubChem itself is not asked: its answers are a saved record and answers made up here.)
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { createXtbEngine } from '../app/src/engine/xtb.js';
import { minimize } from '../app/src/relax.js';
import * as E from '../app/src/editor.js';
import { parseSMILES, fromSMILES, toSMILES, isSMILES } from '../app/src/smiles.js';
import { fromRecord, findCompound } from '../app/src/pubchem.js';
import { composition } from '../app/src/scene.js';
import { BOHR_ANGSTROM, NUMBERS } from '../app/src/units.js';

const root = new URL('../', import.meta.url);
const xtb = await createXtbEngine(readFileSync(new URL('app/engine/xtb_wasm.wasm', root)));
const library = JSON.parse(readFileSync(new URL('app/src/library.json', root))).molecules;
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };

const formula = (m) => composition({ elements: m.elements, positions: Float64Array.from(m.positions) });
const dist = (m, i, j) => Math.hypot(m.positions[3 * i] - m.positions[3 * j], m.positions[3 * i + 1] - m.positions[3 * j + 1], m.positions[3 * i + 2] - m.positions[3 * j + 2]) * BOHR_ANGSTROM;
const graph = (m) => m.bonds.map(([i, j]) => `${Math.min(i, j)}-${Math.max(i, j)}`).sort().join(' ');
const orders = (m) => m.bonds.filter(([i, j]) => m.elements[i] !== 'H' && m.elements[j] !== 'H').map((b) => b[2]).sort().join('');
const closest = (m) => { let d = 9; for (let i = 0; i < m.elements.length; i++) for (let j = i + 1; j < m.elements.length; j++) d = Math.min(d, dist(m, i, j)); return d; };
const worstBond = (m) => Math.max(...m.bonds.map(([i, j, o]) => Math.abs(dist(m, i, j) / BOHR_ANGSTROM / E.bondLength(m.elements[i], m.elements[j], o) - 1)));
function relax(m) {
  const numbers = m.elements.map((el) => NUMBERS[el]), unpaired = (numbers.reduce((s, z) => s + z, 0) - m.charge) % 2;
  const walk = minimize((x) => {
    const r = xtb.compute({ numbers, positions: x, charge: m.charge, etemp: 300, forces: true, orbitals: false, unpaired });
    return r.ok ? { e: r.energy, g: r.forces.map((f) => -f) } : { e: Infinity };
  }, m.positions, { gtol: 5e-5, maxIter: 600 });
  for (;;) { const s = walk.next(); if (s.done) return { ...s.value, mol: { ...m, positions: Array.from(s.value.x) } }; }
}
const fromLibrary = (key) => { const m = library.find((x) => x.key === key); return E.withBonds(m.elements, m.xyz.map((v) => v / BOHR_ANGSTROM), m.charge ?? 0); };

console.log('reading SMILES');
for (const [text, want, bondOrders, charge = 0] of [
  ['C', 'CH4', ''], ['CCO', 'C2H6O', '11'], ['C=C', 'C2H4', '2'], ['C#N', 'CHN', '3'], ['CC(=O)O', 'C2H4O2', '112'],
  ['c1ccccc1', 'C6H6', '111222'], ['c1ccncc1', 'C5H5N', '111222'], ['c1cc[nH]c1', 'C4H5N', '11122'], ['c1ccoc1', 'C4H4O', '11122'],
  ['c1ccc2ccccc2c1', 'C10H8', '11111122222'], ['O=c1cccc[nH]1', 'C5H5NO', '1111222'], ['c1ccc(cc1)-c1ccccc1', 'C12H10', '1111111222222'],
  ['Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C8H10N4O2', '111111111112222'], ['CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'C8H10N4O2', '111111111112222'],
  ['[NH4+]', 'H4N', '', 1], ['[OH-]', 'HO', '', -1], ['C[N+](=O)[O-]', 'CH3NO2', '112'], ['[Na+].[Cl-]', 'ClNa', ''], ['[CH3]', 'CH3', ''],
  ['N[C@@H](C)C(=O)O', 'C3H7NO2', '11112'], ['F/C=C/F', 'C2H2F2', '112'], ['OS(=O)(=O)O', 'H2O4S', '1122'], ['[H][H]', 'H2', ''],
  ['C%10CC%10', 'C3H6', '111'], ['[13CH4]', 'CH4', ''], ['CCO ethanol', 'C2H6O', '11'], ['[Fe+2]', 'Fe', '', 2], ['ClC(Cl)(Cl)Cl', 'CCl4', '1111'],
]) {
  const m = fromSMILES(text);
  ok(formula(m) === want && orders(m) === bondOrders && m.charge === charge, `${text}: ${want}${charge ? `, charge ${charge}` : ''}${bondOrders ? `, bond orders ${bondOrders}` : ''}`);
}
for (const [text, reason] of [['CI', /not among the 36/], ['c1ccnc1', /\[nH\]/], ['C(C', /not closed/], ['C1CC', /never closed/], ['water', /position 1/], ['', /empty/],
  ['C=', /leads nowhere/], ['C$C', /quadruple/], ['[Xe]', /not among/], ['C11', /bonded twice/], ['H2O', /\[H\]/], ['C'.repeat(33), /up to 100 atoms; this one has 101/]]) {
  assert.throws(() => fromSMILES(text), reason);
  checks++;
}
ok(true, 'text that is not SMILES, elements outside the 36 and molecules that are too large are refused, each with its reason');
ok(isSMILES('CO') && isSMILES('c1ccccc1') && !isSMILES('aspirin') && !isSMILES('co') && !isSMILES('2244') && !isSMILES('Benzene'), 'a name is told from SMILES: aspirin, co, 2244 and Benzene are names');

console.log('the geometry built from it');
for (const [name, text] of [['ethanol', 'CCO'], ['isobutane', 'CC(C)C'], ['benzene', 'c1ccccc1'], ['naphthalene', 'c1ccc2ccccc2c1'], ['caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'],
  ['aspirin', 'CC(=O)Oc1ccccc1C(=O)O'], ['cubane', 'C12C3C4C1C5C2C3C45'], ['adamantane', 'C1C2CC3CC1CC(C2)C3'], ['norbornadiene', 'C1=CC2C=CC1C2'], ['spiropentane', 'C1CC11CC1']]) {
  const m = fromSMILES(text);
  ok(closest(m) > 0.85 && worstBond(m) < 0.2 && graph(E.withBonds(m.elements, m.positions)) === graph(m),
    `${name}: no bond more than ${(100 * worstBond(m)).toFixed(0)}% off its length, closest pair ${closest(m).toFixed(2)} A, the atoms near each other are exactly the bonded ones`);
}
for (const [name, text, atoms] of [['dotriacontane, a chain of 32 carbons', 'C'.repeat(32), 98], ['cholesterol', 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 74], ['a branched alkane', 'CC(C)(C)CC(C)(C)CC(C)(C)CC(C)(C)CC(C)(C)C', 65]]) {
  const m = fromSMILES(text);
  ok(m.elements.length === atoms && closest(m) > 0.85 && worstBond(m) < 0.2 && graph(E.withBonds(m.elements, m.positions)) === graph(m),
    `${name}, ${atoms} atoms: nothing runs into anything (closest pair ${closest(m).toFixed(2)} A), the atoms near each other are exactly the bonded ones`);
}
{
  const chain = fromSMILES('C'.repeat(32)), cs = chain.elements.map((el, k) => (el === 'C' ? k : -1)).filter((k) => k >= 0);
  const span = Math.max(...cs.flatMap((a) => cs.map((b) => dist(chain, a, b))));
  ok(span > 38, `the chain is laid out straight, ${span.toFixed(1)} A from end to end`);
}
{
  const salt = fromSMILES('[Na+].[Cl-]'), pair = fromSMILES('O.O');
  ok(Math.abs(dist(salt, 0, 1) - 2.5) < 0.01 && closest(pair) > 0.9 && pair.bonds.length === 4, `separate pieces stand side by side: Na and Cl ${dist(salt, 0, 1).toFixed(2)} A apart, two waters apart and whole`);
}

console.log('relaxing it with the engine');
for (const [name, text, key] of [['ethanol', 'CCO', 'c2h5oh'], ['benzene', 'c1ccccc1', 'c6h6'], ['pyrrole', 'c1cc[nH]c1', 'c4h5n'], ['pyridine', 'c1ccncc1', 'c5h5n'],
  ['furan', 'c1ccoc1', 'c4h4o'], ['ammonium', '[NH4+]', 'nh4+'], ['sulfuric acid', 'OS(=O)(=O)O', 'h2so4'], ['hydrogen cyanide', 'C#N', 'hcn']]) {
  const r = relax(fromSMILES(text)), de = r.e - library.find((x) => x.key === key).energy;
  ok(r.converged && Math.abs(de) < 2e-5, `${name} from ${text}: ${r.it} steps, energy ${de >= 0 ? '+' : ''}${de.toExponential(1)} Eh from the library's`);
}
for (const [name, text] of [['cubane', 'C12C3C4C1C5C2C3C45'], ['adamantane', 'C1C2CC3CC1CC(C2)C3'], ['norbornadiene', 'C1=CC2C=CC1C2'], ['caffeine', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C']]) {
  const m = fromSMILES(text), r = relax(m);
  ok(r.converged && graph(E.withBonds(m.elements, r.mol.positions)) === graph(m), `${name}: relaxed in ${r.it} steps with every bond in place and no new ones`);
}

console.log('writing SMILES');
for (const text of ['C', 'CCO', 'C=C', 'C#N', 'CC(=O)O', 'C1=CC=CC=C1', '[NH4+]', '[OH-]', 'C[N+](=O)[O-]', '[CH3]', 'OS(=O)(=O)O', '[H][H]', 'O=O', 'C1CC1', 'CC(C)(C)C', 'FC=CF', '[Na+].[Cl-]', 'CC(=O)[O-]'])
  assert.equal(toSMILES(fromSMILES(text)), text), checks++;
ok(true, 'simple molecules come back letter for letter (CCO, CC(=O)O, C1=CC=CC=C1, [NH4+], C[N+](=O)[O-], [CH3], ...)');
for (const text of ['c1ccc2ccccc2c1', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'C12C3C4C1C5C2C3C45', 'C1C2CC3CC1CC(C2)C3', 'CC(=O)Oc1ccccc1C(=O)O', 'c1ccc(cc1)-c1ccccc1', 'C1CC11CC1']) {
  const a = fromSMILES(text), b = parseSMILES(toSMILES(a)), count = (g) => g.atoms.reduce((s, x) => s + 1 + x.h, 0);
  ok(count(b) === a.elements.length && b.bonds.length === a.bonds.filter(([i, j]) => a.elements[i] !== 'H' && a.elements[j] !== 'H').length && b.bonds.map((x) => x[2]).sort().join('') === orders(a),
    `${text} -> ${toSMILES(a)}: read again, the same atoms, bonds and bond orders`);
}
for (const [key, want] of [['c6h6', 'C1=CC=CC=C1'], ['c5h5n', 'N1=CC=CC=C1'], ['c4h5n', 'N1C=CC=C1'], ['co', '[C-]#[O+]'], ['h3o+', '[OH3+]'], ['oh-', '[OH-]'], ['ch3', '[CH3]'], ['glycine', 'NCC(=O)O']])
  ok(toSMILES(fromLibrary(key)) === want, `library ${key}, whose bonds are read from the distances: ${want}`);
for (const [key, want] of [['ch3coo-', 'C([O-])(=O)C'], ['no3-', 'N([O-])(=O)=O'], ['no2', '[N](=O)=O'], ['hno3', 'N(=O)(=O)O']])
  ok(toSMILES(fromLibrary(key)) === want, `library ${key}, whose like bonds are all read as double: ${want}`);
{
  // every molecule of the library can be written, and what is written reads back as the same atoms and charge
  const wrong = library.filter((entry) => {
    const m = fromLibrary(entry.key), g = parseSMILES(toSMILES(m));
    return g.atoms.reduce((s, a) => s + 1 + a.h, 0) !== m.elements.length || g.charge !== m.charge;
  }).map((entry) => entry.key);
  ok(wrong.length === 0, `all ${library.length} molecules of the library are written, and each reads back with the same atoms and charge${wrong.length ? ` (not: ${wrong.join(' ')})` : ''}`);
}

console.log('PubChem');
{
  const record = JSON.parse(readFileSync(new URL('tests/cases/pubchem-ethanol-3d.json', root)));
  const m = fromRecord(record.PC_Compounds[0]);
  ok(formula(m) === 'C2H6O' && m.bonds.length === 8 && m.charge === 0 && worstBond(m) < 0.1 && toSMILES(m) === 'OCC', `a 3D record of PubChem (ethanol, CID 702) is taken with its bonds: ${toSMILES(m)}, bonds within ${(100 * worstBond(m)).toFixed(0)}% of the editor's lengths`);
  const r = relax(m), de = r.e - library.find((x) => x.key === 'c2h5oh').energy;
  // PubChem keeps the gauche conformer; the library's ethanol is the anti one, a little higher in GFN2-xTB
  ok(r.converged && de < 0 && de > -4e-3 && graph(E.withBonds(m.elements, r.mol.positions)) === graph(m), `relaxed from there in ${r.it} steps to its own conformer, ${de.toExponential(1)} Eh from the library's`);

  // PubChem's answers, as they were on the day this was written
  const asked = [];
  const answers = {
    'name/ethanol/property/Title,SMILES/JSON': [200, { PropertyTable: { Properties: [{ CID: 702, SMILES: 'CCO', Title: 'Ethanol' }] } }],
    'cid/702/property/Title,SMILES/JSON': [200, { PropertyTable: { Properties: [{ CID: 702, SMILES: 'CCO', Title: 'Ethanol' }] } }],
    'cid/702/record/JSON?record_type=3d': [200, record],
    'name/sodium%20chloride/property/Title,SMILES/JSON': [200, { PropertyTable: { Properties: [{ CID: 5234, SMILES: '[Na+].[Cl-]', Title: 'Sodium Chloride' }] } }],
    'cid/5234/record/JSON?record_type=3d': [404, { Fault: { Code: 'PUGREST.NotFound' } }],
    'name/qwertyzzz/property/Title,SMILES/JSON': [404, { Fault: { Code: 'PUGREST.NotFound' } }],
    'name/busy/property/Title,SMILES/JSON': [503, {}],
  };
  const fetch = async (url) => {
    const path = url.replace('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/', '');
    asked.push(path);
    if (path.startsWith('name/offline')) throw new TypeError('Failed to fetch');
    const [status, body] = answers[path] ?? [500, {}];
    return { status, ok: status === 200, json: async () => body };
  };
  let found = await findCompound(' ethanol ', { fetch });
  ok(found.cid === 702 && found.name === 'Ethanol' && formula(found.molecule) === 'C2H6O' && asked.length === 2, 'a compound by name: two requests, its title, its 3D structure');
  found = await findCompound('702', { fetch });
  ok(found.name === 'Ethanol' && found.molecule.elements.length === 9, 'the same by its number');
  found = await findCompound('sodium chloride', { fetch });
  ok(found.molecule === null && formula(fromSMILES(found.smiles)) === 'ClNa', 'a salt has no 3D structure there: its SMILES comes instead and is built here');
  await assert.rejects(findCompound('qwertyzzz', { fetch }), /knows no compound/);
  await assert.rejects(findCompound('busy', { fetch }), /busy/);
  await assert.rejects(findCompound('offline', { fetch }), /cannot be reached/);
  await assert.rejects(findCompound('  ', { fetch }), /Type the name/);
  ok(true, 'an unknown name, a busy server, no connection and an empty field each give their own message');
  const big = { atoms: { aid: Array.from({ length: 101 }, (_, k) => k + 1), element: new Array(101).fill(6) }, coords: [{ aid: [], conformers: [{ x: [], y: [], z: [] }] }] };
  assert.throws(() => fromRecord(big), /up to 100 atoms; this one has 101/);
  assert.throws(() => fromRecord({ ...record.PC_Compounds[0], atoms: { ...record.PC_Compounds[0].atoms, element: [53, 6, 6, 1, 1, 1, 1, 1, 1] } }), /number 53/);
  ok(true, 'a compound that is too large or has an element past krypton is refused');
}
console.log(`\n${checks} checks passed`);
