import { createEngineClient } from './engine/client.js';
import { FogRenderer, project, cameraBasis, turnCamera } from './render/fog.js';
import * as P from './plane.js';
import { drawSlice } from './slice2d.js';
import { drawEnergy } from './chart.js';
import { buildScene, translated, rotated, pieces, composition } from './scene.js';
import { createEditMode, savedMolecules } from './editmode.js';
import { EXPERIMENTS } from './experiments.js';
import { fixSigns, frontier, firstShown, overlapMatrix, followOrbitals } from './orbitals.js';
import { chargeSites, potentialScale, ESP_FIXED_SCALE } from './esp.js';
import { BOHR_ANGSTROM, HARTREE_EV, NUMBERS, COVALENT, VALENCE } from './units.js';

const $ = (id) => document.getElementById(id);
const params = new URLSearchParams(location.search);
const canvas = $('view'), labels = $('labels'), statusEl = $('status');
const setStatus = (text, error = false) => { statusEl.textContent = text; statusEl.classList.toggle('error', error); };

const reduceMotion = window.matchMedia?.('(prefers-reduced-motion: reduce)').matches;

const state = {
  mol: params.get('mol') || 'h2o',
  second: params.get('second') || '', // key of the second molecule in the scene, or ''
  unpaired: 0, // number of unpaired electrons in the scene (set from the scene, changeable in the settings)
  experiment: null, // the chosen entry of EXPERIMENTS, until the molecules are changed by hand
  body: params.get('body') !== '0',
  // what is drawn: the selected orbital; 'esp', the electrostatic potential on the density
  // surface; 'equi', the surfaces where the potential is plus and minus equiVolts
  view: ['esp', 'equi'].includes(params.get('view')) ? params.get('view') : 'orbital',
  equiVolts: Math.max(0.1, Math.min(10, +params.get('ev') || 1)),
  espRange: params.get('scale') === 'own' ? 'own' : 'fixed', // the colours of the potential: on the scale that is the same for every molecule, or stretched to the scene's own
  espScale: 0.08, // Eh per unit charge at which the colours are full when stretched to the scene; set for every new scene
  wobble: params.has('wobble') ? params.get('wobble') === '1' : !reduceMotion,
  selected: -1,
  result: null,
  molecule: null, // the scene: see scene.js
  camera: { yaw: +(params.get('yaw') ?? 0.6), pitch: +(params.get('pitch') ?? 0.3), distance: 10 },
  // cutting plane: geometry (see plane.js) or null, and how it is shown
  plane: null,
  planeMode: params.get('pmode') === 'sheet' ? 'sheet' : 'knife',
  planeField: params.get('pfield') === 'rho' ? 'rho' : 'psi',
  mouse: 'camera', // what dragging does: 'camera', 'shift' or 'tilt' the plane
  picking: null, // { kind: 'three' | 'along' | 'across', atoms: [] } while atoms are being chosen
  // dynamics: settings, and the running session (null while the molecule is static)
  settings: {
    temperature: +(params.get('temp') ?? 300), // K, target of the thermostat
    thermostat: params.get('thermostat') !== '0',
    speed: params.has('speed') ? +params.get('speed') : 30, // fs of simulation per second; 0 = as fast as possible
    dt: +(params.get('dt') ?? 0.5), // fs
    tau: 100, // fs, coupling time of the thermostat
    etemp: 300, // K, electronic temperature
    wall: true, // the soft wall around the scene during dynamics
    wallRadius: null, // its radius set by hand, angstrom; null: by the size of the scene, see wallFor
  },
  // wallRadius: of the wall, bohr, also while it is switched off; extent: the radius the nuclei are drawn within
  md: null, // { session, running, wallRadius, extent, planeRadius, time, steps, epot, ekin, etot, temperature, rate }
  frames: 0, // dynamics frames shown so far
};
let shownCamera = state.camera; // the camera of the last drawn frame, wobble included

let freshScene = false; // set when an experiment is chosen: its spin setting applies even if the molecules stay the same
let renderer, library, dirty = true, overlay = false, drag = null, last = 0, slow = 0, lastMolKey = null;
// The cloud is drawn from the orbital grid while the picture changes. Once it has stood still
// for STILL_MS it is drawn once more, exactly: `refined` says that this has been done.
const STILL_MS = 200;
let drawnAt = 0, refined = false;
let scale = Math.min(window.devicePixelRatio || 1, 1.25);
let resolveReady;
const ready = new Promise((r) => (resolveReady = r));

// test hook: ?canvas=496x420 pins the canvas size
if (params.has('canvas')) {
  const [w, h] = params.get('canvas').split('x');
  canvas.style.cssText = `width:${w}px;height:${h}px;flex:none`;
}

try {
  // test hooks: ?bodygrid=0 evaluates the body density exactly at every ray step,
  // ?orbgrid=0 the selected orbital as well (every frame is then an exact one)
  renderer = new FogRenderer(canvas, { bodyGrid: params.get('bodygrid') !== '0', orbGrid: params.get('orbgrid') !== '0' });
} catch (e) {
  setStatus(e.message, true);
  throw e;
}
canvas.addEventListener('webglcontextlost', (e) => {
  e.preventDefault();
  setStatus('The browser lost its graphics context (WebGL). Reload the page; if that does not help, restart the browser.', true);
});
const engine = createEngineClient();

const withBohr = (m) => ({ ...m, bohr: Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM) });
// a molecule in both lists, under the heading of its group; one that is there already is renamed
function listMolecule(m) {
  for (const select of [$('mol'), ...(m.optimized !== false ? [$('mol2')] : [])]) {
    const old = [...select.options].find((o) => o.value === m.key);
    if (old) { old.text = m.name; continue; }
    const name = m.group || 'Other';
    let group = [...select.querySelectorAll('optgroup')].find((g) => g.label === name);
    if (!group) { group = document.createElement('optgroup'); group.label = name; select.append(group); }
    group.append(new Option(m.name, m.key));
  }
}

async function loadLibrary() {
  const lib = await (await fetch(new URL('./library.json', import.meta.url))).json();
  // the library, then the molecules the user has made in the editor
  library = [...lib.molecules, ...savedMolecules()].map(withBohr);
  $('mol2').add(new Option('— none —', ''));
  library.forEach(listMolecule);
  if (!library.some((m) => m.key === state.mol)) state.mol = library[0].key;
  $('mol').value = state.mol;
  if (!library.some((m) => m.key === state.second)) state.second = '';
  $('mol2').value = state.second;
  // formulas as people write them, for the "now in the scene" line: taken from the library
  // names (their composition in the program's own order -> the text before the comma)
  for (const m of library) {
    if (m.custom) continue;
    const key = composition({ elements: m.elements, positions: m.bohr });
    if (!key.includes(' ')) FORMULAS[key] ??= m.name.split(',')[0].replace(/[⁺⁻]/g, '');
  }
  const list = $('experiment');
  list.add(new Option('— custom —', ''));
  const groups = new Map();
  for (const e of EXPERIMENTS) {
    if (!library.some((m) => m.key === e.first) || (e.second && !library.some((m) => m.key === e.second))) continue;
    if (!groups.has(e.group)) { const g = document.createElement('optgroup'); g.label = e.group; groups.set(e.group, g); list.append(g); }
    groups.get(e.group).append(new Option(e.title, e.key));
  }
}

// products that are not in the library, by the program's own composition string
const FORMULAS = {
  HO: 'OH', H: 'H', F: 'F', O: 'O', FHLi: 'LiF·H', H2LiO: 'Li·H₂O', CN: 'CN', O2: 'O₂', H5O2: 'H₅O₂', C2H5: 'C₂H₅', CH2N: 'HCNH', HO2: 'HO₂', H2O3: 'H₂O₃',
  FH2NaO: 'NaF·H₂O', FH2LiO: 'LiF·H₂O', CHNaO3: 'NaHCO₃', F2Li2: 'Li₂F₂', F2Na2: 'Na₂F₂', H2Li2: 'Li₂H₂', BFH3: 'BH₃F', CH3BO: 'BH₃CO', BF3H3N: 'BF₃NH₃',
  H3LiO: 'LiH·H₂O', H3NaO: 'NaH·H₂O', H4NaO2: 'NaOH·H₂O', CHLiO2: 'LiH·CO₂', Be2H4: 'Be₂H₄', C2H3O2: 'CH₃COO', CHO2: 'HCOO',
};
const SUB = '₀₁₂₃₄₅₆₇₈₉';
// 'H2O + H4N' -> 'H₂O + NH₄'
const pretty = (text) => text.split(' + ').map((part) => {
  const [, count, f] = part.match(/^(\d+ )?(.+)$/);
  return (count ? count.trim() + ' ' : '') + (FORMULAS[f] ?? f.replace(/\d/g, (d) => SUB[d]));
}).join(' + ');

const TWO_HINT = 'Two molecules: drag any nucleus to move its molecule (with Shift, to rotate it). After a breakup each fragment moves separately. While the dynamics is running, the molecule follows the cursor by the nucleus you hold; point at a nucleus of the other molecule to lead it there.';
const sceneNote = () => {
  const e = state.experiment;
  if (e) return `${e.pull ?? e.todo} What should happen: ${e.expect}`;
  if (state.molecule.nFragments < 2) return state.molecule.note || '';
  return TWO_HINT;
};
const failure = (res) => (res.converged === false ? `The electronic problem did not converge: ${res.error}` : `Calculation error: ${res.error}`);

// builds the scene from the library and shows it: new molecule(s), view reset
async function recompute() {
  const scene = buildScene(library.find((m) => m.key === state.mol), library.find((m) => m.key === state.second));
  // a new scene starts with its own number of unpaired electrons; an experiment may ask for another
  if (freshScene || scene.key !== state.molecule?.key) state.unpaired = state.experiment?.unpaired ?? scene.unpaired;
  freshScene = false;
  setStatus('Calculating…');
  const res = await engine.compute({ numbers: scene.numbers, positions: scene.positions, charge: scene.charge, etemp: state.settings.etemp, unpaired: state.unpaired, forces: false });
  if (!res) return; // superseded by a newer request
  if (!res.ok) { setStatus(failure(res), true); return; }
  fixSigns(res);
  state.result = res;
  state.molecule = scene;
  try {
    renderer.setScene(Float32Array.from(scene.positions), res, { fragments: scene.fragments });
    refreshPotential(true);
    state.camera.distance = framing();
  } catch (e) {
    setStatus(e.message, true);
    return;
  }
  const { lumo } = frontier(res);
  const want = params.get('orb');
  state.selected = want === 'lumo' && lumo >= 0 ? lumo : want && !/^(homo|lumo)$/.test(want) ? Math.min(res.nmo - 1, Math.max(0, +want)) : firstShown(res);
  renderer.select(state.selected);
  if (scene.key !== lastMolKey) { lastMolKey = scene.key; setPlane(null); cancelPicking(); state.settings.wallRadius = null; } // a wall set by hand was for the other scene
  for (const el of document.querySelectorAll('.two-only')) el.hidden = scene.nFragments < 2;
  buildLabels();
  buildLevels();
  buildInfo();
  refreshView();
  setStatus(sceneNote());
  dirty = true;
  refreshPlane();
  refreshMd();
  if (!planeFromUrlDone) { planeFromUrlDone = true; planeFromUrl(); }
  resolveReady();
}

// the same scene with the nuclei somewhere else (dynamics, or moved by hand): orbitals keep
// their identity and sign, the camera stays where it is
function showMoved(positions, res) {
  state.selected = followOrbitals(state.result, res, overlapMatrix(res.shells, positions), state.selected);
  state.result = res;
  state.molecule = { ...state.molecule, positions };
  const md = state.md;
  renderer.setScene(Float32Array.from(positions), res, {
    fragments: state.molecule.fragments, planeRadius: md ? md.planeRadius : renderer.scene.planeRadius, extent: md?.extent,
  });
  renderer.select(state.selected);
  refreshPotential();
  updateLevels();
  refreshPlane();
  buildInfo();
  dirty = true;
  window.chemsim?.onShown?.(); // lets the checks look at every frame
}

// nuclei placed by hand while nothing is running
async function recomputeMoved(positions) {
  state.molecule = { ...state.molecule, positions };
  dirty = true;
  const m = state.molecule;
  const res = await engine.compute({ numbers: m.numbers, positions, charge: m.charge, etemp: state.settings.etemp, unpaired: state.unpaired, forces: false });
  if (!res) return;
  if (!res.ok) { setStatus(failure(res), true); return; }
  if (statusEl.classList.contains('error')) setStatus(sceneNote());
  showMoved(positions, res);
}

const fmt = (v, d) => v.toFixed(d).replace('-', '−');

function buildInfo() {
  const r = state.result, md = state.md;
  const { homo, lumo } = r ? frontier(r) : { homo: -1, lumo: -1 };
  const rows = !r ? [['Method', 'GFN2-xTB'], ['Total energy', '—']] : md ? [
    ['Time', `${md.time.toFixed(md.time < 100 ? 1 : 0)} fs`],
    ['Temperature', `${md.temperature.toFixed(0)} K`],
    ['E potential', `${fmt(md.epot, 5)} Eh`],
    ['E kinetic', `${fmt(md.ekin, 5)} Eh`],
    ['E total', `${fmt(md.etot, 5)} Eh`],
    ['Steps per second', md.running && md.rate ? md.rate.toFixed(0) : '—'],
  ] : [
    ['Method', 'GFN2-xTB'],
    ['Total energy', `${fmt(r.energy, 5)} Eh`],
  ];
  if (homo >= 0 && lumo >= 0) rows.push(['HOMO–LUMO gap', `${fmt((r.eps[lumo] - r.eps[homo]) * HARTREE_EV, 2)} eV`]);
  if (!md && r) {
    rows.push(['Basis functions', String(r.nao)]);
    rows.push(['Calculation time', `${r.ms < 10 ? r.ms.toFixed(1) : r.ms.toFixed(0)} ms`]);
  }
  rows.push(['Now in the scene', pretty(composition(state.molecule))]);
  if (state.unpaired) rows.push(['Unpaired e⁻', String(state.unpaired)]);
  const info = $('info');
  info.textContent = '';
  for (const [k, v] of rows) {
    const dt = document.createElement('dt'), dd = document.createElement('dd');
    dt.textContent = k; dd.textContent = v;
    info.append(dt, dd);
  }
}

// alpha: occupation by spin-up electrons alone, when the engine reports it
function occupationMark(occ, alpha) {
  if (occ > 1.999) return '●●';
  if (occ < 0.001) return '○○';
  if (Math.abs(occ - 1) < 0.005 && alpha > 0.995) return '↑○'; // one electron, unpaired
  return occ.toFixed(2); // fractional occupations are shown as they are
}

function buildLevels() {
  const r = state.result, { homo, lumo } = frontier(r), box = $('levels');
  box.textContent = '';
  for (let m = r.nmo - 1; m >= 0; m--) {
    const b = document.createElement('button');
    b.type = 'button';
    b.dataset.orbital = m;
    const e = document.createElement('span'), t = document.createElement('span');
    e.textContent = `${fmt(r.eps[m] * HARTREE_EV, 2)} eV`;
    t.className = 'tag';
    t.textContent = (m === homo ? 'HOMO ' : m === lumo ? 'LUMO ' : '') + occupationMark(r.occ[m], r.occAlpha?.[m]);
    b.append(e, t);
    b.addEventListener('click', () => selectOrbital(m));
    box.appendChild(b);
  }
  markLevels();
  box.querySelector('[aria-pressed="true"]')?.scrollIntoView({ block: 'center' });
}

// the same list with new numbers, without rebuilding it (every frame of the dynamics)
function updateLevels() {
  const r = state.result, box = $('levels');
  if (box.children.length !== r.nmo) return buildLevels();
  const { homo, lumo } = frontier(r);
  for (const b of box.children) {
    const m = +b.dataset.orbital;
    b.firstChild.textContent = `${fmt(r.eps[m] * HARTREE_EV, 2)} eV`;
    b.lastChild.textContent = (m === homo ? 'HOMO ' : m === lumo ? 'LUMO ' : '') + occupationMark(r.occ[m], r.occAlpha?.[m]);
  }
  markLevels();
}

function markLevels() {
  for (const b of $('levels').children) b.setAttribute('aria-pressed', String(+b.dataset.orbital === state.selected));
}

function selectOrbital(m) {
  state.selected = m;
  renderer.select(m);
  markLevels();
  refreshPlane();
  dirty = true;
}

// the potential at which the colours are full now
const espFull = () => (state.espRange === 'own' ? state.espScale : ESP_FIXED_SCALE);
// the camera distance that frames what is drawn
const framing = () => renderer.scene[{ esp: 'espDistance', equi: 'equiDistance' }[state.view] ?? 'distance'];

// The charges the potential is drawn from, for what is in the scene now: worked out when the
// potential is on the screen, or (rescale) for a new scene, which also sets the scale of the
// colours. The scale then stays while the nuclei of the scene move.
function refreshPotential(rescale = false) {
  if (!state.result || (!rescale && state.view === 'orbital')) return;
  const m = state.molecule, sites = chargeSites(m.numbers, m.positions, state.result);
  if (rescale) state.espScale = potentialScale(m.numbers, m.positions, sites);
  renderer.setPotential(sites);
  dirty = true;
}

// pushes what is drawn (orbital or potential) to its switch and to the legend
function refreshView() {
  const esp = state.view === 'esp', equi = state.view === 'equi';
  for (const b of $('view-mode').children) b.setAttribute('aria-pressed', String(b.dataset.v === state.view));
  $('body').disabled = esp || equi; // the surface of the cloud is drawn in place of the body
  $('legend-orbital').hidden = esp || equi;
  $('legend-esp').hidden = !esp;
  $('legend-equi').hidden = !equi;
  $('equi-field').hidden = !equi;
  // the slider runs from 0.1 to 10 V, evenly in the logarithm
  const ev = state.equiVolts, shown = ev.toFixed(ev < 1 ? 2 : 1);
  $('equi-value').value = Math.round(50 * Math.log10(ev / 0.1));
  $('equi-out').textContent = `±${shown} V`;
  $('equi-neg').textContent = `−${shown} V`;
  $('equi-pos').textContent = `+${shown} V`;
  $('esp-range-field').hidden = !esp;
  for (const b of $('esp-range').children) b.setAttribute('aria-pressed', String(b.dataset.v === state.espRange));
  const volts = espFull() * HARTREE_EV;
  $('esp-min').textContent = `−${volts.toFixed(volts < 1 ? 2 : 1)} V`;
  $('esp-max').textContent = `+${volts.toFixed(volts < 1 ? 2 : 1)} V`;
  dirty = true;
}

// ---- cutting plane
let planeFromUrlDone = false;
const atomPosition = (k) => Array.from(state.molecule.positions.subarray(3 * k, 3 * k + 3));

// pushes the plane state to the renderer, the controls and the flat plot
function refreshPlane() {
  const pl = state.plane;
  renderer.setPlane(pl && { normal: pl.normal, d: P.planeOffset(pl), mode: state.planeMode, field: state.planeField });
  $('plane-on').checked = !!pl;
  $('plane-controls').hidden = !pl;
  $('slice-box').hidden = !pl;
  for (const [id, value] of [['plane-mode', state.planeMode], ['plane-field', state.planeField], ['plane-mouse', state.mouse], ['plane-snap', state.picking?.kind]])
    for (const b of $(id).children) b.setAttribute('aria-pressed', String(b.dataset.v === value));
  canvas.classList.toggle('picking', !!state.picking);
  if (pl && state.result) drawSlice($('slice'), $('slice-caption'), renderer, pl, state.planeField, state.molecule);
  dirty = true;
}

function setPlane(pl) {
  state.plane = pl;
  if (!pl) state.mouse = 'camera';
  refreshPlane();
}

// keep the plane inside the drawn disc
function movePlane(pl) {
  const limit = renderer.scene.planeRadius - 0.5, off = P.planeOffset(pl);
  setPlane(Math.abs(off) > limit ? P.shifted(pl, Math.sign(off) * limit - off) : pl);
}

const SNAP = {
  three: { need: 3, hint: 'Pick three atoms', make: (a, cam) => P.planeThroughThree(a[0], a[1], a[2], cam) },
  along: { need: 2, hint: 'Pick two atoms: the plane will run along the bond between them', make: (a, cam) => P.planeAlongBond(a[0], a[1], cam) },
  across: { need: 2, hint: 'Pick two atoms: the plane will stand across the bond between them', make: (a, cam) => P.planeAcrossBond(a[0], a[1], cam) },
};

function markPicked() {
  for (let k = 0; k < labels.children.length; k++) labels.children[k].classList.toggle('picked', !!state.picking?.atoms.includes(k));
}

function startPicking(kind) {
  if (state.molecule.elements.length < SNAP[kind].need) { setStatus('This molecule does not have enough atoms for this snap.', true); return; }
  state.picking = { kind, atoms: [] };
  setStatus(`${SNAP[kind].hint}. Click the nuclei (0 of ${SNAP[kind].need}).`);
  markPicked();
  refreshPlane();
}

function cancelPicking() {
  if (!state.picking) return;
  state.picking = null;
  markPicked();
}

function snapTo(kind, atoms) {
  const pl = SNAP[kind].make(atoms.map(atomPosition), cameraBasis(shownCamera));
  if (!pl) { setStatus('These three atoms lie on one line, so no plane through them is defined.', true); return false; }
  setPlane(pl);
  return true;
}

function pickAtom(k) {
  const pk = state.picking;
  if (pk.atoms.includes(k)) return;
  pk.atoms.push(k);
  markPicked();
  if (pk.atoms.length < SNAP[pk.kind].need) {
    setStatus(`${SNAP[pk.kind].hint}. Click the nuclei (${pk.atoms.length} of ${SNAP[pk.kind].need}).`);
    return;
  }
  state.picking = null;
  markPicked();
  if (snapTo(pk.kind, pk.atoms)) setStatus(sceneNote());
  else refreshPlane();
}

// nearest nucleus to a point of the canvas (CSS pixels), within `radius` px
function atomAt(x, y, radius = 28) {
  const w = canvas.clientWidth, h = canvas.clientHeight;
  let best = -1, bestD = radius;
  for (let k = 0; k < state.molecule.elements.length; k++) {
    const [px, py] = project(atomPosition(k), shownCamera, w, h), d = Math.hypot(px - x, py - y);
    if (d < bestD) { bestD = d; best = k; }
  }
  return best;
}

// test hook: ?plane=three:0,1,2 | along:0,1 | across:0,1 | view  and  &poffset=<angstrom>
function planeFromUrl() {
  const spec = params.get('plane');
  if (!spec) return;
  const [kind, list] = spec.split(':');
  if (kind === 'view') setPlane(P.planeFacingCamera(cameraBasis(shownCamera)));
  else if (SNAP[kind]) snapTo(kind, list.split(',').map(Number));
  if (state.plane && params.has('poffset')) movePlane(P.shifted(state.plane, +params.get('poffset') / BOHR_ANGSTROM));
}

function buildLabels() {
  labels.textContent = '';
  for (const el of state.molecule.elements) {
    const s = document.createElement('span');
    s.textContent = el;
    labels.appendChild(s);
  }
}

function placeLabels(camera) {
  const w = canvas.clientWidth, h = canvas.clientHeight, p = state.molecule.positions;
  const ox = canvas.offsetLeft, oy = canvas.offsetTop;
  for (let k = 0; k < labels.children.length; k++) {
    const [x, y] = project(p.subarray(3 * k, 3 * k + 3), camera, w, h);
    labels.children[k].style.transform = `translate(${(ox + x + 7).toFixed(1)}px,${(oy + y - 19).toFixed(1)}px)`;
  }
}

function resize() {
  const w = Math.max(1, Math.round(canvas.clientWidth * scale)), h = Math.max(1, Math.round(canvas.clientHeight * scale));
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; dirty = true; }
}

// exact: the orbital evaluated at every ray step rather than read from its grid (see fog.js)
function render(ts, exact = false) {
  const wob = state.wobble && !drag && !edit.active;
  const camera = wob ? turnCamera(state.camera, 0.22 * Math.sin(ts * 0.0008), 0) : { ...state.camera };
  if (state.result && !(edit.active && !edit.cloud)) renderer.draw(camera, { body: state.body, exact, view: state.view, espScale: state.view === 'equi' ? state.equiVolts / HARTREE_EV : espFull() }); else renderer.clear();
  placeLabels(camera);
  shownCamera = camera;
  if (edit.active) edit.draw(camera);
  // the wall: a sphere about the centre, which is a circle on the screen. Shown while there
  // is a run for it to act in, or while its setting is in sight.
  {
    const R = state.settings.wall && state.result && !edit.active && (state.md || $('md-settings').open) ? (state.md ? state.md.wallRadius : wallFor(state.molecule.positions)) : 0;
    const D = camera.distance, h = canvas.clientHeight, c = $('wall-ring').firstElementChild;
    c.setAttribute('r', R && D > 1.02 * R ? 2.3 * h / 2 * R / Math.sqrt(D * D - R * R) : 0);
    c.setAttribute('cx', canvas.offsetLeft + canvas.clientWidth / 2);
    c.setAttribute('cy', canvas.offsetTop + h / 2);
  }
  // the spring between the held nucleus and the cursor
  const line = $('grab-line'), hold = drag?.mode === 'hold' && drag.cursor;
  line.classList.toggle('on', !!hold);
  if (hold) {
    const w = canvas.clientWidth, h = canvas.clientHeight, ox = canvas.offsetLeft, oy = canvas.offsetTop;
    const [ax, ay] = project(atomPosition(drag.atom), camera, w, h);
    // the line ends where the nucleus is being led: at the nucleus it is aimed at (ringed), or at the cursor
    const [tx, ty] = drag.to >= 0 ? project(atomPosition(drag.to), camera, w, h) : drag.cursor;
    const l = line.firstElementChild, c = line.lastElementChild;
    l.setAttribute('x1', ox + ax); l.setAttribute('y1', oy + ay); l.setAttribute('x2', ox + tx); l.setAttribute('y2', oy + ty);
    c.setAttribute('cx', ox + tx); c.setAttribute('cy', oy + ty); c.setAttribute('r', drag.to >= 0 ? 11 : 4);
  }
}

function frame(ts) {
  resize();
  const wob = state.wobble && !drag && !edit.active;
  // frames that follow one another: the wobble, the view being turned by hand, the editor at work
  const live = wob || (dirty && (drag || edit.active));
  if ((state.result || edit.active) && (dirty || wob)) {
    // if continuous rendering is too slow (under some 28 frames a second), lower the resolution
    if (live && last) { slow = ts - last > 36 ? slow + 1 : 0; if (slow > 6 && scale > 0.55) { scale *= 0.8; slow = 0; resize(); } }
    last = live ? ts : 0;
    render(ts);
    dirty = false;
    overlay = false;
    drawnAt = ts;
    refined = !renderer.orbGrid || state.view !== 'orbital'; // without the grid the frame is exact already; the potential has no exact frame
  } else if (!refined && !drag && !state.md?.running && ts - drawnAt > STILL_MS) {
    // nothing has changed for a while: the same frame once more, exactly
    last = 0;
    render(ts, true);
    refined = true;
    overlay = false;
  } else {
    last = 0;
    // the editor's skeleton alone (the cursor moved over an atom): the cloud is not redrawn for that
    if (overlay && edit.active) edit.draw(shownCamera);
    overlay = false;
  }
  requestAnimationFrame(frame);
}

// ---- dynamics
const speedFromSlider = (v) => (v >= 100 ? 0 : Math.round(2 * 10 ** (3 * v / 100))); // 2 … 2000 fs/s, then "max"
const sliderFromSpeed = (fs) => (fs <= 0 ? 100 : Math.round(100 * Math.log10(fs / 2) / 3));
const mdParams = () => ({ dt: state.settings.dt, thermostat: state.settings.thermostat, temperature: state.settings.temperature, tau: state.settings.tau });
let sessions = 0;

// pushes the dynamics state to its controls
function refreshMd() {
  const st = state.settings, md = state.md;
  $('md-run').textContent = md?.running ? '⏸ Pause' : '▶ Start';
  $('md-run').setAttribute('aria-pressed', String(!!md?.running));
  $('md-reset').disabled = !md;
  $('md-chart').hidden = !md;
  $('md-temp').value = st.temperature;
  $('md-temp-out').textContent = `${st.temperature} K`;
  $('md-speed').value = sliderFromSpeed(st.speed);
  $('md-speed-out').textContent = st.speed > 0 ? `${st.speed} fs/s` : 'max';
  $('md-thermo').checked = st.thermostat;
  $('md-dt').value = st.dt; $('md-tau').value = st.tau; $('md-etemp').value = st.etemp;
  $('md-wall').checked = st.wall;
  $('md-wall-r').disabled = !st.wall;
  if (state.molecule) $('md-wall-r').value = ((md ? md.wallRadius : wallFor(state.molecule.positions)) * BOHR_ANGSTROM).toFixed(1);
  // unpaired electrons: only the values of the right parity make sense
  const pick = $('unpaired'), n = state.molecule?.electrons ?? 0;
  pick.textContent = '';
  for (let u = n % 2; u <= Math.min(n, 6); u += 2) pick.add(new Option(u === 0 ? '0 — all paired' : String(u), u));
  pick.value = state.unpaired;
}

// how far from the centre the farthest nucleus is, bohr
function farthest(p) {
  let mr = 0;
  for (let k = 0; k < p.length; k += 3) mr = Math.max(mr, Math.hypot(p[k], p[k + 1], p[k + 2]));
  return mr;
}
// The radius of the soft wall for nuclei at p, bohr: the one set by hand, but not inside the
// nuclei; otherwise just inside what the camera frames for the molecule at rest.
function wallFor(p) {
  const set = state.settings.wallRadius;
  return set ? Math.max(set / BOHR_ANGSTROM, farthest(p) + 0.5) : farthest(p) + 2.5;
}
// what the dynamics is told: no wall is a wall infinitely far
const wallSent = () => (state.settings.wall ? state.md.wallRadius : Infinity);
// the wall of the running session has changed: tell the dynamics, keep the drawn region around it
function wallChanged() {
  const md = state.md;
  if (md) {
    if (state.settings.wall) md.extent = md.wallRadius + 1;
    engine.md('set', { session: md.session, params: { wallRadius: wallSent() } });
  }
  refreshMd();
  dirty = true;
}
// Without the wall a piece that comes off flies away for good. The scene glides after what is
// left (the worker keeps it about the origin), the camera is not touched, and a piece that is
// DROP bohr from everything else is taken out of the run.
const DROP = 10 / BOHR_ANGSTROM;
function followPieces(f) {
  const md = state.md;
  md.extent = Math.max(md.extent, 1.3 * (farthest(f.positions) + 1)); // only the region that is drawn; the zoom stays
  if (md.dropping || md.keepAll || !f.result?.coeff) return;
  const m = state.molecule, n = m.numbers.length, p = f.positions, parts = pieces({ ...m, positions: p });
  if (parts.count < 2) return;
  // pieces closer than DROP to one another stay together; the largest such group is what is kept
  const group = Array.from({ length: parts.count }, (_, k) => k);
  const find = (k) => (group[k] === k ? k : (group[k] = find(group[k])));
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) {
    if (parts.piece[i] !== parts.piece[j] && Math.hypot(p[3 * i] - p[3 * j], p[3 * i + 1] - p[3 * j + 1], p[3 * i + 2] - p[3 * j + 2]) < DROP) group[find(parts.piece[i])] = find(parts.piece[j]);
  }
  const size = new Map();
  for (let i = 0; i < n; i++) size.set(find(parts.piece[i]), (size.get(find(parts.piece[i])) ?? 0) + 1);
  if (size.size < 2) return;
  const best = [...size].reduce((a, b) => (b[1] > a[1] ? b : a))[0]; // of equal ones, the one with the first atom
  const keep = [], gone = [];
  for (let i = 0; i < n; i++) (find(parts.piece[i]) === best ? keep : gone).push(i);
  // The charge that leaves: the electrons the engine counts on those atoms, to a whole number.
  // The unpaired electrons that stay: the spin the orbitals put on the atoms kept, to the
  // nearest number of the right parity.
  const r = f.result, charge = m.charge - Math.round(gone.reduce((q, i) => q + r.charges[i], 0));
  const electrons = keep.reduce((e, i) => e + VALENCE[m.numbers[i]], 0) - charge;
  let spin = 0;
  if (state.unpaired) {
    const S = overlapMatrix(r.shells, p), stays = new Uint8Array(n), atomOf = [];
    for (const i of keep) stays[i] = 1;
    for (const sh of r.shells) for (let k = 0; k < 2 * sh.l + 1; k++) atomOf.push(sh.atom);
    for (let o = 0; o < r.nmo; o++) {
      const w = 2 * r.occAlpha[o] - r.occ[o];
      if (Math.abs(w) < 1e-9) continue;
      for (let i = 0; i < r.nao; i++) {
        if (!stays[atomOf[i]]) continue;
        let t = 0;
        for (let j = 0; j < r.nao; j++) t += S[i * r.nao + j] * r.coeff[o * r.nao + j];
        spin += w * r.coeff[o * r.nao + i] * t;
      }
    }
  }
  let unpaired = electrons % 2;
  while (unpaired + 2 <= electrons && Math.abs(unpaired + 2 - Math.abs(spin)) < Math.abs(unpaired - Math.abs(spin))) unpaired += 2;
  md.dropping = pretty(composition({ elements: gone.map((i) => m.elements[i]), positions: Float64Array.from(gone.flatMap((i) => [p[3 * i], p[3 * i + 1], p[3 * i + 2]])) }));
  engine.md('drop', { session: md.session, keep, charge, unpaired });
}

// the worker has left out what flew away: the scene is now the atoms f.kept
function applyDrop(f) {
  const m = state.molecule, md = state.md, keep = f.kept, fragments = Uint8Array.from(keep, (i) => m.fragments[i]);
  const two = fragments.includes(0) && fragments.includes(1), numbers = keep.map((i) => m.numbers[i]);
  state.molecule = {
    ...m, elements: keep.map((i) => m.elements[i]), numbers, positions: f.positions, charge: f.charge,
    fragments: two ? fragments : new Uint8Array(keep.length), nFragments: two ? 2 : 1,
    electrons: numbers.reduce((e, z) => e + VALENCE[z], 0) - f.charge,
  };
  state.unpaired = f.unpaired;
  state.result = null; // the orbitals of the smaller scene are new ones
  for (const k of ['t', 'etot', 'epot']) md.history[k] = []; // its energy is another energy
  cancelPicking();
  drag = null;
  buildLabels();
  for (const el of document.querySelectorAll('.two-only')) el.hidden = !two;
  setStatus(`${md.dropping} has flown 10 Å away from the rest and is taken out of the scene. “Reset” brings everything back.`);
  md.dropping = null;
}

function mdStart() {
  const p = state.molecule.positions, wallRadius = wallFor(p);
  state.md = {
    session: ++sessions, running: false, wallRadius, extent: wallRadius + 1, planeRadius: renderer.scene.planeRadius,
    time: 0, steps: 0, epot: state.result.energy, ekin: 0, etot: state.result.energy, temperature: 0, rate: 0,
    history: { t: [], etot: [], epot: [] },
  };
  engine.md('init', {
    session: state.md.session, numbers: state.molecule.numbers, positions: p, charge: state.molecule.charge,
    etemp: state.settings.etemp, unpaired: state.unpaired, params: mdParams(), speed: state.settings.speed, wallRadius: wallSent(),
    seed: params.has('seed') ? +params.get('seed') : (Math.random() * 2 ** 31) | 0, thermalize: true,
  });
}

function mdCommand(cmd, data = {}) {
  if (!state.md) mdStart();
  engine.md(cmd, { session: state.md.session, ...data });
  if (cmd === 'run' || cmd === 'pause' || cmd === 'step') { state.md.running = cmd === 'run'; setStatus(sceneNote()); refreshMd(); }
}

function mdStop() {
  if (!state.md) return;
  engine.md('stop', { session: state.md.session });
  state.md = null;
}

engine.onFrame = (f) => {
  const md = state.md;
  if (f.fatal) { setStatus(f.error, true); if (md) md.running = false; refreshMd(); return; }
  if (!md || f.session !== md.session) return;
  // the frame sent on 'init' predates the 'run' command that follows it: it does not stop anything
  const running = f.steps === 0 && !f.error ? md.running : f.running, changed = running !== md.running;
  md.running = running;
  if (f.error) {
    if (md.dropping) { md.dropping = null; md.keepAll = true; } // what is left cannot be solved alone: nothing is taken out in this run
    setStatus(f.converged === false
      ? `The electronic problem did not converge (${f.error}). The dynamics has been stopped: it cannot go on with unconverged forces. This usually means that atoms came too close together or that an occupied and an empty level have almost the same energy. Raise the electronic temperature under “Settings” (to 1500–3000 K, say) and press “Start” again, or press “Reset”.`
      : `Calculation error: ${f.error}. The dynamics has been stopped.`, true);
    refreshMd();
    return;
  }
  if (f.kept) applyDrop(f);
  Object.assign(md, { time: f.time, steps: f.steps, epot: f.epot, ekin: f.ekin, etot: f.etot, temperature: f.temperature, rate: f.rate });
  const hist = md.history;
  if (!hist.t.length || f.time > hist.t[hist.t.length - 1]) {
    hist.t.push(f.time); hist.etot.push(f.etot); hist.epot.push(f.epot);
    if (hist.t.length > 1500) for (const k of ['t', 'etot', 'epot']) hist[k] = hist[k].filter((_, i) => i % 2 === 0); // keep it light
  }
  drawEnergy($('md-chart'), hist);
  if (!state.settings.wall) followPieces(f);
  if (f.result?.coeff) { f.result.ms = 0; showMoved(f.positions, f.result); }
  if (f.kept) { selectOrbital(firstShown(state.result)); refreshMd(); }
  buildInfo();
  if (changed) refreshMd();
  state.frames++;
  dirty = true;
};

// ---- controls
const leaveExperiment = () => {
  // settings an experiment asked for do not outlive it
  if (state.experiment) Object.assign(state.settings, { etemp: 300, temperature: 300 });
  state.experiment = null;
  $('experiment').value = '';
};
$('mol').addEventListener('change', (e) => { state.mol = e.target.value; params.delete('orb'); leaveExperiment(); mdStop(); recompute(); });
$('mol2').addEventListener('change', (e) => { state.second = e.target.value; params.delete('orb'); leaveExperiment(); mdStop(); recompute(); });
// an experiment sets the molecules and everything the run needs
function startExperiment(key) {
  const e = EXPERIMENTS.find((x) => x.key === key);
  state.experiment = e ?? null;
  $('experiment').value = e ? key : '';
  if (!e) { setStatus(sceneNote()); return; }
  state.mol = e.first; state.second = e.second || '';
  $('mol').value = state.mol; $('mol2').value = state.second;
  Object.assign(state.settings, { temperature: e.temperature ?? 300, etemp: e.etemp ?? 300, speed: e.speed ?? 60, thermostat: true });
  freshScene = true; // the scene's (or the experiment's) number of unpaired electrons is taken afresh
  params.delete('orb');
  mdStop();
  return recompute();
}
$('experiment').addEventListener('change', (e) => startExperiment(e.target.value));
$('body').checked = state.body;
$('wobble').checked = state.wobble;
$('body').addEventListener('change', (e) => { state.body = e.target.checked; dirty = true; });
$('wobble').addEventListener('change', (e) => { state.wobble = e.target.checked; dirty = true; });
for (const b of $('view-mode').children) b.addEventListener('click', () => {
  if (b.dataset.v === state.view) return;
  state.view = b.dataset.v; // the camera stays where it is: the views are compared at one zoom
  refreshPotential();
  refreshView();
});
$('equi-value').addEventListener('input', (e) => { state.equiVolts = 0.1 * 10 ** (+e.target.value / 50); refreshView(); });
for (const b of $('esp-range').children) b.addEventListener('click', () => { state.espRange = b.dataset.v; refreshView(); });
$('about-open').addEventListener('click', () => $('about').showModal());

$('md-run').addEventListener('click', () => mdCommand(state.md?.running ? 'pause' : 'run'));
$('md-step').addEventListener('click', () => mdCommand('step'));
$('md-reset').addEventListener('click', () => { mdStop(); recompute(); });
const setting = (id, key, read, apply) => $(id).addEventListener('input', (e) => {
  const v = read(e.target);
  if (typeof v === 'number' && !Number.isFinite(v)) return;
  state.settings[key] = v;
  refreshMd();
  apply?.();
});
const sendParams = () => { if (state.md) engine.md('set', { session: state.md.session, params: mdParams() }); };
setting('md-temp', 'temperature', (t) => +t.value, sendParams);
setting('md-thermo', 'thermostat', (t) => t.checked, sendParams);
setting('md-speed', 'speed', (t) => speedFromSlider(+t.value), () => { if (state.md) engine.md('set', { session: state.md.session, speed: state.settings.speed }); });
for (const [id, key, lo, hi] of [['md-dt', 'dt', 0.05, 2], ['md-tau', 'tau', 5, 5000], ['md-etemp', 'etemp', 0, 30000]]) {
  $(id).addEventListener('change', (e) => {
    const v = Math.max(lo, Math.min(hi, +e.target.value || state.settings[key]));
    state.settings[key] = v;
    refreshMd();
    if (key !== 'etemp') return sendParams();
    // the electronic temperature also changes the static picture
    if (state.md) engine.md('set', { session: state.md.session, etemp: v }); else recompute();
  });
}

$('md-wall').addEventListener('change', (e) => {
  state.settings.wall = e.target.checked;
  // put back around pieces that have gone far, the wall must not slam them inwards
  if (state.md && state.settings.wall) state.md.wallRadius = Math.max(state.md.wallRadius, farthest(state.molecule.positions) + 2.5);
  wallChanged();
});
$('md-wall-r').addEventListener('change', (e) => {
  const v = +e.target.value;
  if (Number.isFinite(v) && v > 0) {
    state.settings.wallRadius = Math.min(30, v);
    if (state.md) state.md.wallRadius = wallFor(state.molecule.positions);
    // a wall asked for is a wall to be seen: the camera backs away until its circle is in the picture
    const R = state.md ? state.md.wallRadius : wallFor(state.molecule.positions);
    state.camera.distance = Math.max(state.camera.distance, 2.6 * R);
  }
  wallChanged();
});
$('md-settings').addEventListener('toggle', () => { dirty = true; }); // the wall is shown while the settings are open

$('unpaired').addEventListener('change', (e) => {
  state.unpaired = +e.target.value;
  if (state.md) {
    // the running session goes on with the new spin state, from the same nuclei
    engine.md('set', { session: state.md.session, unpaired: state.unpaired });
    engine.md('move', { session: state.md.session, positions: state.molecule.positions, wallRadius: wallSent() });
  } else recomputeMoved(state.molecule.positions);
  buildInfo();
});

$('plane-on').addEventListener('change', (e) => {
  cancelPicking();
  setPlane(e.target.checked ? P.planeFacingCamera(cameraBasis(shownCamera)) : null);
  setStatus(sceneNote());
});
const segmented = (id, onPick) => { for (const b of $(id).children) b.addEventListener('click', () => onPick(b.dataset.v)); };
segmented('plane-mode', (v) => { state.planeMode = v; refreshPlane(); });
segmented('plane-field', (v) => { state.planeField = v; refreshPlane(); });
segmented('plane-mouse', (v) => { state.mouse = v; refreshPlane(); });
segmented('plane-snap', (v) => {
  if (state.picking?.kind === v) { cancelPicking(); setStatus(sceneNote()); refreshPlane(); } else startPicking(v);
});
$('plane-flip').addEventListener('click', () => setPlane(P.flipped(state.plane)));
$('plane-reset').addEventListener('click', () => { cancelPicking(); setPlane(P.planeFacingCamera(cameraBasis(shownCamera))); });
window.addEventListener('keydown', (e) => {
  if (e.key === 'Escape' && state.picking) { cancelPicking(); setStatus(sceneNote()); refreshPlane(); }
});
new ResizeObserver(() => { if (state.plane && state.result) refreshPlane(); }).observe($('slice'));

// ---- the molecule editor (editmode.js): what it needs from the page
let shows = 0;
// The molecule under the editor's hands, with the cloud the engine gives for it (or the one
// that came with a relaxation frame). The camera stays where it is and only backs off when
// the molecule outgrows the view. A structure the engine cannot solve is shown without a
// cloud: the skeleton alone. Resolves with true when the cloud is up to date.
async function showEdited(scene, { result, follow } = {}) {
  const seq = ++shows, n = scene.electrons;
  if ((state.unpaired - n) % 2 !== 0 || state.unpaired > n) state.unpaired = n % 2;
  state.molecule = scene;
  buildLabels();
  dirty = true;
  let res = result;
  if (!res) {
    res = await engine.compute({ numbers: scene.numbers, positions: scene.positions, charge: scene.charge, etemp: state.settings.etemp, unpaired: state.unpaired, forces: false });
    if (!res || seq !== shows) return false; // a newer state is on its way
    if (!res.ok) {
      state.result = null;
      $('levels').textContent = '';
      buildInfo();
      refreshMd();
      setStatus(`${failure(res)}. The cloud is hidden until the structure can be calculated; the structure itself can still be edited.`, true);
      dirty = true;
      return false;
    }
  }
  const prev = state.result, same = follow && prev && prev.nao === res.nao && prev.nmo === res.nmo;
  if (same) state.selected = followOrbitals(prev, res, overlapMatrix(res.shells, scene.positions), state.selected);
  else { fixSigns(res); state.selected = firstShown(res); }
  state.result = res;
  try {
    renderer.setScene(Float32Array.from(scene.positions), res);
    // a molecule being built gets its scale anew whenever its atoms change, not while they only move
    refreshPotential(!same);
    if (!same) refreshView();
    if (framing() > state.camera.distance) state.camera.distance = framing();
  } catch (e) {
    setStatus(e.message, true);
    return false;
  }
  renderer.select(state.selected);
  if (same) updateLevels(); else buildLevels();
  buildInfo();
  refreshMd();
  dirty = true;
  return true;
}

const EDIT_LOCKS = '.bar .field, #edit-open, #edit-new, #find-open, #md-panel, #plane-panel';
const edit = createEditMode({
  state, canvas, engine, setStatus, turnView,
  show: showEdited,
  camera: () => shownCamera,
  // everything, or (cloud: false) only the skeleton over a cloud that has not changed
  redraw: ({ cloud = true } = {}) => { if (cloud) dirty = true; else overlay = true; },
  unproject: (x, y, depth) => unproject(x, y, depth),
  library: () => library,
  name: (scene) => pretty(composition(scene)),
  barHeight: (h) => { statusEl.style.top = h ? `${h + 20}px` : ''; },
  // everything that works on a finished scene stands still while the molecule is being made
  enter() {
    mdStop();
    leaveExperiment();
    cancelPicking();
    setPlane(null);
    for (const el of document.querySelectorAll(EDIT_LOCKS)) el.inert = true;
    for (const el of document.querySelectorAll('.two-only')) el.hidden = true;
    $('stage').classList.add('editing');
  },
  // back to an ordinary scene: the molecule with this key, or the one that was there before
  leave(key) {
    $('stage').classList.remove('editing');
    shows++;
    if (key) { state.mol = key; state.second = ''; }
    if (!library.some((m) => m.key === state.mol)) state.mol = library[0].key;
    if (!library.some((m) => m.key === state.second)) state.second = '';
    $('mol').value = state.mol; $('mol2').value = state.second;
    params.delete('orb');
    // the controls come back once the scene they work on is there
    const unlock = () => { for (const el of document.querySelectorAll(EDIT_LOCKS)) el.inert = false; };
    return recompute().then(unlock, unlock);
  },
  addMolecule(entry) {
    const m = withBohr(entry), k = library.findIndex((x) => x.key === m.key);
    if (k >= 0) library[k] = m; else library.push(m);
    listMolecule(m);
  },
  removeMolecule(key) {
    library = library.filter((m) => m.key !== key);
    for (const select of [$('mol'), $('mol2')]) {
      const o = [...select.options].find((x) => x.value === key), group = o?.parentElement;
      o?.remove();
      if (group?.tagName === 'OPTGROUP' && !group.children.length) group.remove();
    }
  },
});

// ---- mouse on the canvas: rotate the view, move the cutting plane, pick nuclei, and
// take hold of a nucleus (a spring while the dynamics runs, the whole molecule otherwise)
const local = (e) => { const r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };

// dragging by (dx, dy) pixels turns the view about the axes of the screen, without end
function turnView(dx, dy) {
  Object.assign(state.camera, turnCamera(state.camera, -dx * 0.008, -dy * 0.008));
  dirty = true;
}

// the point under the cursor at distance `depth` in front of the camera
function unproject(x, y, depth) {
  const cam = cameraBasis(shownCamera), h = canvas.clientHeight, u = (x - canvas.clientWidth / 2) / (h / 2), v = -(y - h / 2) / (h / 2);
  return [0, 1, 2].map((c) => cam.ro[c] + (cam.fw[c] * 2.3 + cam.rt[c] * u + cam.up[c] * v) * depth / 2.3);
}

// While the dynamics runs, the held nucleus is led by the cursor. Over a nucleus of another
// molecule the target is that nucleus, wherever it is in depth: the point at bonding distance
// from it. Anywhere else the target is the point under the cursor in the plane that faces
// the screen and passed through the held nucleus when it was taken.
function holdAt(x, y) {
  drag.cursor = [x, y];
  const w = canvas.clientWidth, h = canvas.clientHeight, a = drag.atom, numbers = state.molecule.numbers;
  let to = -1, near = 22;
  for (let k = 0; k < numbers.length; k++) {
    if (drag.piece.includes(k)) continue;
    const [px, py] = project(atomPosition(k), shownCamera, w, h), d = Math.hypot(px - x, py - y);
    if (d < near) { near = d; to = k; }
  }
  // a single molecule has nowhere to be carried: only the spring on the nucleus acts
  const grab = { atom: a, piece: drag.piece.length < numbers.length ? drag.piece : undefined };
  if (to >= 0) { grab.toAtom = to; grab.rest = (COVALENT[numbers[a]] + COVALENT[numbers[to]]) / BOHR_ANGSTROM; }
  else grab.target = unproject(x, y, drag.depth);
  engine.md('set', { session: state.md.session, grab });
  if (to !== drag.to) {
    const el = state.molecule.elements;
    setStatus(to >= 0 ? `Leading ${el[a]} to ${el[to]}: it is pulled straight towards that nucleus, to bonding distance.`
      : `Leading ${el[a]} to the cursor, in the plane of the screen. To bring it to an atom of the other molecule, point at that atom's nucleus.`);
  }
  drag.to = to;
  dirty = true;
}

// a molecule put somewhere else by hand
function placeNuclei(positions) {
  for (let k = 0; k < positions.length; k += 3) if (Math.hypot(positions[k], positions[k + 1], positions[k + 2]) > 16) return; // not out of sight
  if (!state.md) return recomputeMoved(positions);
  state.molecule = { ...state.molecule, positions };
  state.md.wallRadius = Math.max(state.md.wallRadius, wallFor(positions));
  state.md.extent = Math.max(state.md.extent, state.md.wallRadius + 1);
  engine.md('move', { session: state.md.session, positions, wallRadius: wallSent() });
  dirty = true;
}

canvas.addEventListener('pointerdown', (e) => {
  if (edit.active) return edit.pointerdown(e);
  const [x, y] = local(e);
  drag = { x: e.clientX, y: e.clientY, moved: 0, mode: 'view' };
  try { canvas.setPointerCapture(e.pointerId); } catch {} // synthetic events in the checks have no real pointer
  const k = state.picking || !state.result ? -1 : atomAt(x, y, 18);
  // what is moved by hand is the piece the nucleus is bonded into right now, not the
  // molecule it once belonged to; with a single piece there is nothing to rearrange
  const parts = k >= 0 ? pieces(state.molecule) : null;
  if (k >= 0 && (state.md?.running || parts.count > 1)) {
    drag.atom = k;
    drag.depth = project(atomPosition(k), shownCamera, canvas.clientWidth, canvas.clientHeight)[2];
    if (state.md?.running) { drag.mode = 'hold'; drag.piece = Array.from(parts.piece.keys()).filter((i) => parts.piece[i] === parts.piece[k]); holdAt(x, y); }
    else { drag.mode = e.shiftKey ? 'turn' : 'move'; drag.member = (a) => parts.piece[a] === parts.piece[k]; }
  }
  if (!state.picking) canvas.style.cursor = 'grabbing';
});
canvas.addEventListener('pointermove', (e) => {
  if (edit.active) return edit.pointermove(e);
  const [x, y] = local(e);
  if (!drag) {
    // show that a nucleus can be taken hold of
    const can = state.result && !state.picking && atomAt(x, y, 18) >= 0 && (state.md?.running || pieces(state.molecule).count > 1);
    canvas.style.cursor = can ? 'pointer' : '';
    return;
  }
  const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
  drag.x = e.clientX; drag.y = e.clientY; drag.moved += Math.abs(dx) + Math.abs(dy);
  const cam = cameraBasis(shownCamera);
  if (drag.mode === 'hold') {
    if (state.md?.running) holdAt(x, y);
    return;
  }
  if (drag.mode === 'move') {
    const perPixel = drag.depth / (2.3 * canvas.clientHeight / 2);
    placeNuclei(translated(state.molecule, drag.member, [0, 1, 2].map((c) => (cam.rt[c] * dx - cam.up[c] * dy) * perPixel)));
    return;
  }
  if (drag.mode === 'turn') {
    placeNuclei(rotated(state.molecule, drag.member, cam.up, dx * 0.01, cam.rt, dy * 0.01));
    return;
  }
  const action = !state.plane ? 'camera' : e.shiftKey ? 'shift' : e.altKey ? 'tilt' : state.mouse;
  if (action === 'camera') {
    turnView(dx, dy);
    dirty = true;
    return;
  }
  const pl = state.plane;
  if (action === 'shift') {
    // drag along the direction in which the normal points on the screen; if the plane
    // faces the camera, dragging up pulls it towards the viewer
    const nx = P.dot(pl.normal, cam.rt), ny = P.dot(pl.normal, cam.up), len = Math.hypot(nx, ny);
    const perPixel = shownCamera.distance / (2.3 * canvas.clientHeight / 2);
    const pixels = len > 0.35 ? (dx * nx - dy * ny) / len : -dy * Math.sign(P.dot(pl.normal, cam.ro) || 1);
    movePlane(P.shifted(pl, pixels * perPixel));
  } else {
    movePlane(P.tilted(P.tilted(pl, cam.up, dx * 0.008), cam.rt, dy * 0.008));
  }
});
const endDrag = (e) => {
  if (edit.active) return edit.pointerup(e);
  if (drag?.mode === 'hold') { if (state.md) engine.md('set', { session: state.md.session, grab: null }); setStatus(sceneNote()); }
  if (drag && drag.moved < 4 && state.picking && e.type === 'pointerup') {
    const [x, y] = local(e), k = atomAt(x, y);
    if (k >= 0) pickAtom(k);
  }
  drag = null;
  canvas.style.cursor = '';
  $('grab-line').classList.remove('on');
  dirty = true;
};
canvas.addEventListener('pointerup', endDrag);
canvas.addEventListener('pointercancel', endDrag);
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  // with Ctrl, among the surfaces of equal potential: the wheel changes their value (up: higher)
  if (e.ctrlKey && state.view === 'equi') {
    state.equiVolts = Math.max(0.1, Math.min(10, state.equiVolts * Math.exp(-e.deltaY * 0.001)));
    refreshView();
    return;
  }
  const base = renderer.scene ? framing() : 10;
  // as far out as it takes to see the whole region of a run: the wall, or the pieces flying apart without it
  const out = Math.max(2.5 * base, 2.6 * (state.md?.extent ?? 0));
  state.camera.distance = Math.max(0.5 * base, Math.min(out, state.camera.distance * Math.exp(e.deltaY * 0.001)));
  dirty = true;
}, { passive: false });
window.addEventListener('keydown', (e) => {
  if (!state.result || (e.key !== 'ArrowUp' && e.key !== 'ArrowDown') || !$('levels').contains(document.activeElement)) return;
  e.preventDefault();
  const m = Math.max(0, Math.min(state.result.nmo - 1, state.selected + (e.key === 'ArrowUp' ? 1 : -1)));
  selectOrbital(m);
  $('levels').querySelector(`[data-orbital="${m}"]`)?.focus();
});

// ---- hooks for automated checks
window.chemsim = {
  ready,
  state,
  renderer,
  selectOrbital,
  frame: (f) => engine.onFrame(f), // feeds a dynamics frame by hand (checks of the error path)
  startExperiment,
  edit,
  composition: () => pretty(composition(state.molecule)),
  md: { run: () => mdCommand('run'), pause: () => mdCommand('pause'), step: () => mdCommand('step'), reset: () => { mdStop(); return recompute(); } },
  // field values on the cut (size x size, bottom row first) with the geometry needed to
  // recompute them elsewhere; null if the browser cannot read float pixels. field: what the
  // plane shows, or 'esp' for the electrostatic potential as the shader sums it
  sliceField(size = 96, field = state.planeField) {
    const pl = state.plane, centre = P.planeCentre(pl), e2 = P.planeUp(pl), half = renderer.scene.planeRadius;
    const values = renderer.renderSlice(size, centre, pl.e1, e2, half, field, true);
    return values && { size, centre, e1: pl.e1, e2, half, field, values };
  },
  // the flat plot of the cut as a PNG data URL
  sliceSnapshot: () => $('slice').toDataURL('image/png'),
  // renders one frame without wobble and returns it as a PNG data URL: the exact frame, or
  // with exact = false the one drawn from the orbital grid
  snapshot(exact = true) { resize(); render(0, exact); dirty = true; return canvas.toDataURL('image/png'); },
  // Median time of n full frames in ms, GPU work included: ms for a frame drawn from the
  // orbital grid (what is shown while the picture changes), exactMs for an exact one, gridMs
  // and orbGridMs for filling the body and the orbital grid.
  benchmark(n = 20) {
    const gl = renderer.gl, px = new Uint8Array(4);
    resize();
    const median = (exact) => {
      const times = [];
      render(0, exact); // the first frame of a kind builds its shader and fills its grid
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
      for (let i = 0; i < n; i++) {
        const t = performance.now();
        render(0, exact);
        gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); // waits for the GPU
        times.push(performance.now() - t);
      }
      return times.sort((a, b) => a - b)[n >> 1];
    };
    const ms = median(false), exactMs = median(true);
    dirty = true;
    let orbGridMs = 0;
    if (renderer.orbGrid) {
      const t = performance.now();
      for (let i = 0; i < 5; i++) { renderer._fillOrbital(); gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); }
      orbGridMs = (performance.now() - t) / 5;
    }
    let gridMs = 0;
    if (renderer.bodyGrid) {
      const t = performance.now();
      for (let i = 0; i < 5; i++) { renderer._fillGrid(); gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); }
      gridMs = (performance.now() - t) / 5;
    }
    return { gridMs, orbGridMs, ms, exactMs, orbGrid: renderer.orbGrid, width: canvas.width, height: canvas.height, gpu: renderer.rendererName(), bodyGrid: renderer.bodyGrid, gridSize: renderer.grid.n };
  },
};

await loadLibrary();
if (params.has('experiment')) await startExperiment(params.get('experiment')); // test hook
else await recompute();
if (params.get('md') === 'run') mdCommand('run'); // test hook
if (params.has('edit')) edit.enter(params.get('edit') === 'new'); // test hook: ?edit=1 opens the editor on the scene, ?edit=new on a new molecule
requestAnimationFrame(frame);
