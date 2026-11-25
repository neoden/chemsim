// Molecules from PubChem (https://pubchem.ncbi.nlm.nih.gov), through its PUG REST interface.
// Nothing here runs unless the user asks for a compound: the page itself requests nothing
// from other sites.
//
// A compound is looked up by name or by its number (CID). PubChem keeps a computed 3D
// structure for most small molecules: it comes with the bonds and is taken as it is. Salts,
// metal compounds and the like have none; for those the SMILES is returned, and the editor
// builds a geometry from it.
import { SYMBOLS, BOHR_ANGSTROM } from './units.js';
import { MAX_ATOMS } from './editor.js';

const API = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound';

// one compound record of PubChem (JSON, record_type=3d) -> an editor molecule
export function fromRecord(record) {
  const numbers = record?.atoms?.element, ids = record?.atoms?.aid, conformer = record?.coords?.[0]?.conformers?.[0];
  if (!numbers?.length || !conformer?.z) throw new Error('PubChem sent a record without a 3D structure.');
  if (numbers.length > MAX_ATOMS) throw new Error(`The editor takes molecules of up to ${MAX_ATOMS} atoms; this one has ${numbers.length}.`);
  const row = new Map(ids.map((id, k) => [id, k])), place = record.coords[0].aid.map((id) => row.get(id));
  const elements = numbers.map((z) => {
    if (!(z in SYMBOLS)) throw new Error(`This compound has an element (number ${z}) that is not among the 86 (H to Rn) the program works with.`);
    return SYMBOLS[z];
  });
  const positions = new Array(3 * elements.length).fill(0);
  place.forEach((k, i) => { positions[3 * k] = conformer.x[i] / BOHR_ANGSTROM; positions[3 * k + 1] = conformer.y[i] / BOHR_ANGSTROM; positions[3 * k + 2] = conformer.z[i] / BOHR_ANGSTROM; });
  const b = record.bonds ?? {};
  // orders above 3 are PubChem's marks for dative and complex bonds
  const bonds = (b.aid1 ?? []).map((id, c) => [row.get(id), row.get(b.aid2[c]), b.order[c] >= 1 && b.order[c] <= 3 ? b.order[c] : 1]);
  return { elements, positions, bonds, charge: record.charge ?? 0 };
}

// Looks a compound up. `query` is a name ("aspirin") or a CID ("2244").
// Resolves with { cid, name, smiles, molecule }, molecule being null when PubChem has no 3D
// structure for the compound; rejects with an Error whose message is meant for the user.
export async function findCompound(query, { fetch = globalThis.fetch, signal } = {}) {
  const text = String(query).trim();
  if (!text) throw new Error('Type the name of a compound.');
  const get = async (path) => {
    let r;
    try { r = await fetch(`${API}/${path}`, { signal }); } catch (e) {
      if (e.name === 'AbortError') throw e;
      throw new Error('PubChem cannot be reached. Check the connection, or build the molecule from SMILES instead.');
    }
    if (r.status === 404) return null;
    if (r.status === 503 || r.status === 429) throw new Error('PubChem is busy right now. Try again in a moment.');
    if (!r.ok) throw new Error(`PubChem answered with an error (${r.status}).`);
    return r.json();
  };
  const by = /^\d+$/.test(text) ? `cid/${text}` : `name/${encodeURIComponent(text)}`;
  const found = await get(`${by}/property/Title,SMILES/JSON`);
  const info = found?.PropertyTable?.Properties?.[0];
  if (!info) throw new Error(`PubChem knows no compound called “${text}”.`);
  const record = await get(`cid/${info.CID}/record/JSON?record_type=3d`);
  const compound = record?.PC_Compounds?.[0];
  return { cid: info.CID, name: info.Title ?? text, smiles: info.SMILES ?? '', molecule: compound ? fromRecord(compound) : null };
}
