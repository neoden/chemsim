import { createEngineClient } from './engine/client.js';
import { FogRenderer, project, cameraBasis } from './render/fog.js';
import * as P from './plane.js';
import { drawSlice } from './slice2d.js';
import { drawEnergy } from './chart.js';
import { buildScene, translated, rotated, pieces, composition } from './scene.js';
import { EXPERIMENTS } from './experiments.js';
import { fixSigns, frontier, overlapMatrix, followOrbitals } from './orbitals.js';
import { BOHR_ANGSTROM, HARTREE_EV, NUMBERS } from './units.js';

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
  },
  md: null, // { session, running, wallRadius, extent, planeRadius, time, steps, epot, ekin, etot, temperature, rate }
  frames: 0, // dynamics frames shown so far
};
let shownCamera = state.camera; // the camera of the last drawn frame, wobble included

let freshScene = false; // set when an experiment is chosen: its spin setting applies even if the molecules stay the same
let renderer, library, dirty = true, drag = null, last = 0, slow = 0, lastMolKey = null;
let scale = Math.min(window.devicePixelRatio || 1, 1.25);
let resolveReady;
const ready = new Promise((r) => (resolveReady = r));

// test hook: ?canvas=496x420 pins the canvas size
if (params.has('canvas')) {
  const [w, h] = params.get('canvas').split('x');
  canvas.style.cssText = `width:${w}px;height:${h}px;flex:none`;
}

try {
  // test hook: ?bodygrid=0 evaluates the body density exactly at every ray step
  renderer = new FogRenderer(canvas, { bodyGrid: params.get('bodygrid') !== '0' });
} catch (e) {
  setStatus(e.message, true);
  throw e;
}
canvas.addEventListener('webglcontextlost', (e) => {
  e.preventDefault();
  setStatus('The browser lost its graphics context (WebGL). Reload the page; if that does not help, restart the browser.', true);
});
const engine = createEngineClient();

async function loadLibrary() {
  const lib = await (await fetch(new URL('./library.json', import.meta.url))).json();
  library = lib.molecules.map((m) => ({ ...m, bohr: Float64Array.from(m.xyz, (v) => v / BOHR_ANGSTROM) }));
  // both lists, grouped under the headings of the library
  const fill = (select, entries) => {
    const groups = new Map();
    for (const m of entries) {
      const name = m.group || 'Other';
      if (!groups.has(name)) { const g = document.createElement('optgroup'); g.label = name; groups.set(name, g); select.append(g); }
      groups.get(name).append(new Option(m.name, m.key));
    }
  };
  fill($('mol'), library);
  if (!library.some((m) => m.key === state.mol)) state.mol = library[0].key;
  $('mol').value = state.mol;
  $('mol2').add(new Option('— none —', ''));
  fill($('mol2'), library.filter((m) => m.optimized !== false));
  if (!library.some((m) => m.key === state.second)) state.second = '';
  $('mol2').value = state.second;
  // formulas as people write them, for the "now in the scene" line: taken from the library
  // names (their composition in the program's own order -> the text before the comma)
  for (const m of library) {
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

const TWO_HINT = 'Two molecules: drag any nucleus to move its molecule (with Shift, to rotate it). After a breakup each fragment moves separately. While the dynamics is running, the nucleus is pulled to the cursor by a spring.';
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
    const view = renderer.setScene(Float32Array.from(scene.positions), res, { fragments: scene.fragments });
    state.camera.distance = view.distance;
  } catch (e) {
    setStatus(e.message, true);
    return;
  }
  const { homo, lumo } = frontier(res);
  const want = params.get('orb');
  state.selected = want === 'lumo' ? lumo : want && want !== 'homo' ? Math.min(res.nmo - 1, Math.max(0, +want)) : homo;
  renderer.select(state.selected);
  if (scene.key !== lastMolKey) { lastMolKey = scene.key; setPlane(null); cancelPicking(); }
  for (const el of document.querySelectorAll('.two-only')) el.hidden = scene.nFragments < 2;
  buildLabels();
  buildLevels();
  buildInfo();
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
  const r = state.result, { homo, lumo } = frontier(r), md = state.md;
  const rows = md ? [
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
  if (!md) {
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

function render(ts) {
  const wob = state.wobble && !drag;
  const camera = { ...state.camera, yaw: state.camera.yaw + (wob ? 0.22 * Math.sin(ts * 0.0008) : 0) };
  renderer.draw(camera, { body: state.body });
  placeLabels(camera);
  shownCamera = camera;
  // the spring between the held nucleus and the cursor
  const line = $('grab-line'), hold = drag?.mode === 'hold' && drag.cursor;
  line.classList.toggle('on', !!hold);
  if (hold) {
    const [ax, ay] = project(atomPosition(drag.atom), camera, canvas.clientWidth, canvas.clientHeight), ox = canvas.offsetLeft, oy = canvas.offsetTop;
    const l = line.firstElementChild, c = line.lastElementChild;
    l.setAttribute('x1', ox + ax); l.setAttribute('y1', oy + ay); l.setAttribute('x2', ox + drag.cursor[0]); l.setAttribute('y2', oy + drag.cursor[1]);
    c.setAttribute('cx', ox + drag.cursor[0]); c.setAttribute('cy', oy + drag.cursor[1]);
  }
}

function frame(ts) {
  resize();
  const wob = state.wobble && !drag;
  if (state.result && (dirty || wob)) {
    // if continuous rendering is too slow, lower the resolution
    if (wob && last) { slow = ts - last > 45 ? slow + 1 : 0; if (slow > 6 && scale > 0.55) { scale *= 0.8; slow = 0; resize(); } }
    last = wob ? ts : 0;
    render(ts);
    dirty = false;
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
  // unpaired electrons: only the values of the right parity make sense
  const pick = $('unpaired'), n = state.molecule?.electrons ?? 0;
  pick.textContent = '';
  for (let u = n % 2; u <= Math.min(n, 6); u += 2) pick.add(new Option(u === 0 ? '0 — all paired' : String(u), u));
  pick.value = state.unpaired;
}

// the soft wall sits just inside what the camera frames for the molecule at rest
function wallFor(p) {
  let mr = 0;
  for (let k = 0; k < p.length; k += 3) mr = Math.max(mr, Math.hypot(p[k], p[k + 1], p[k + 2]));
  return mr + 2.5;
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
    etemp: state.settings.etemp, unpaired: state.unpaired, params: mdParams(), speed: state.settings.speed, wallRadius,
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
    setStatus(f.converged === false
      ? `The electronic problem did not converge (${f.error}). The dynamics has been stopped: it cannot go on with unconverged forces. This usually means that atoms came too close together or that an occupied and an empty level have almost the same energy. Raise the electronic temperature under “Settings” (to 1500–3000 K, say) and press “Start” again, or press “Reset”.`
      : `Calculation error: ${f.error}. The dynamics has been stopped.`, true);
    refreshMd();
    return;
  }
  Object.assign(md, { time: f.time, steps: f.steps, epot: f.epot, ekin: f.ekin, etot: f.etot, temperature: f.temperature, rate: f.rate });
  const hist = md.history;
  if (!hist.t.length || f.time > hist.t[hist.t.length - 1]) {
    hist.t.push(f.time); hist.etot.push(f.etot); hist.epot.push(f.epot);
    if (hist.t.length > 1500) for (const k of ['t', 'etot', 'epot']) hist[k] = hist[k].filter((_, i) => i % 2 === 0); // keep it light
  }
  drawEnergy($('md-chart'), hist);
  if (f.result?.coeff) { f.result.ms = 0; showMoved(f.positions, f.result); }
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

$('unpaired').addEventListener('change', (e) => {
  state.unpaired = +e.target.value;
  if (state.md) {
    // the running session goes on with the new spin state, from the same nuclei
    engine.md('set', { session: state.md.session, unpaired: state.unpaired });
    engine.md('move', { session: state.md.session, positions: state.molecule.positions, wallRadius: state.md.wallRadius });
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

// ---- mouse on the canvas: rotate the view, move the cutting plane, pick nuclei, and
// take hold of a nucleus (a spring while the dynamics runs, the whole molecule otherwise)
const local = (e) => { const r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };

// the point under the cursor at distance `depth` in front of the camera
function unproject(x, y, depth) {
  const cam = cameraBasis(shownCamera), h = canvas.clientHeight, u = (x - canvas.clientWidth / 2) / (h / 2), v = -(y - h / 2) / (h / 2);
  return [0, 1, 2].map((c) => cam.ro[c] + (cam.fw[c] * 2.3 + cam.rt[c] * u + cam.up[c] * v) * depth / 2.3);
}

function holdAt(x, y) {
  drag.cursor = [x, y];
  engine.md('set', { session: state.md.session, grab: { atom: drag.atom, target: unproject(x, y, drag.depth) } });
  dirty = true;
}

// a molecule put somewhere else by hand
function placeNuclei(positions) {
  for (let k = 0; k < positions.length; k += 3) if (Math.hypot(positions[k], positions[k + 1], positions[k + 2]) > 16) return; // not out of sight
  if (!state.md) return recomputeMoved(positions);
  state.molecule = { ...state.molecule, positions };
  state.md.wallRadius = Math.max(state.md.wallRadius, wallFor(positions));
  state.md.extent = state.md.wallRadius + 1;
  engine.md('move', { session: state.md.session, positions, wallRadius: state.md.wallRadius });
  dirty = true;
}

canvas.addEventListener('pointerdown', (e) => {
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
    if (state.md?.running) { drag.mode = 'hold'; holdAt(x, y); }
    else { drag.mode = e.shiftKey ? 'turn' : 'move'; drag.member = (a) => parts.piece[a] === parts.piece[k]; }
  }
  if (!state.picking) canvas.style.cursor = 'grabbing';
});
canvas.addEventListener('pointermove', (e) => {
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
    state.camera.yaw -= dx * 0.008;
    state.camera.pitch = Math.max(-1.4, Math.min(1.4, state.camera.pitch + dy * 0.008));
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
  if (drag?.mode === 'hold' && state.md) engine.md('set', { session: state.md.session, grab: null });
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
  const base = renderer.scene?.distance ?? 10;
  state.camera.distance = Math.max(0.5 * base, Math.min(2.5 * base, state.camera.distance * Math.exp(e.deltaY * 0.001)));
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
  composition: () => pretty(composition(state.molecule)),
  md: { run: () => mdCommand('run'), pause: () => mdCommand('pause'), step: () => mdCommand('step'), reset: () => { mdStop(); return recompute(); } },
  // field values on the cut (size x size, bottom row first) with the geometry needed to
  // recompute them elsewhere; null if the browser cannot read float pixels
  sliceField(size = 96) {
    const pl = state.plane, centre = P.planeCentre(pl), e2 = P.planeUp(pl), half = renderer.scene.planeRadius;
    const values = renderer.renderSlice(size, centre, pl.e1, e2, half, state.planeField, true);
    return values && { size, centre, e1: pl.e1, e2, half, field: state.planeField, values };
  },
  // the flat plot of the cut as a PNG data URL
  sliceSnapshot: () => $('slice').toDataURL('image/png'),
  // renders one frame without wobble and returns it as a PNG data URL
  snapshot() { resize(); render(0); dirty = true; return canvas.toDataURL('image/png'); },
  // median time of n full frames in ms, GPU work included
  benchmark(n = 20) {
    const gl = renderer.gl, px = new Uint8Array(4), times = [];
    resize();
    for (let i = 0; i < n; i++) {
      const t = performance.now();
      render(0);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); // waits for the GPU
      times.push(performance.now() - t);
    }
    times.sort((a, b) => a - b);
    let gridMs = 0;
    if (renderer.bodyGrid) {
      const t = performance.now();
      for (let i = 0; i < 5; i++) { renderer._fillGrid(); gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px); }
      gridMs = (performance.now() - t) / 5;
    }
    return { gridMs, ms: times[n >> 1], width: canvas.width, height: canvas.height, gpu: renderer.rendererName(), bodyGrid: renderer.bodyGrid, gridSize: renderer.grid.n };
  },
};

await loadLibrary();
if (params.has('experiment')) await startExperiment(params.get('experiment')); // test hook
else await recompute();
if (params.get('md') === 'run') mdCommand('run'); // test hook
requestAnimationFrame(frame);
