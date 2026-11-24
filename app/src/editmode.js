// The molecule editor on the page: the toolbar, the ball-and-stick skeleton over the cloud,
// what the mouse does while editing, undo, relaxation, and the user's own molecules kept in
// the browser. The model and its operations are in editor.js; the page around it (the scene,
// the renderer, the lists of molecules) is reached through `app`, see main.js.
import * as E from './editor.js';
import { project, cameraBasis } from './render/fog.js';
import { looseScene } from './scene.js';
import { BOHR_ANGSTROM } from './units.js';

const PALETTE = ['H', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl'];
const COLOURS = { H: '#e9ebf0', C: '#8f96a3', N: '#5f84ff', O: '#ff5f52', F: '#86e36f', P: '#ff9d42', S: '#f2da46', Cl: '#4fd98a', B: '#f0b6a0', Si: '#d9c29a', Br: '#c0603c' };
const STORE = 'chemsim.molecules';
export const MY_GROUP = 'My molecules';

// the user's molecules: library entries (see tools/build-library.mjs) kept in localStorage
export function savedMolecules() {
  try {
    const list = JSON.parse(localStorage.getItem(STORE)) ?? [];
    return list.filter((m) => m && typeof m.key === 'string' && Array.isArray(m.elements) && Array.isArray(m.xyz) && m.xyz.length === 3 * m.elements.length);
  } catch { return []; }
}
function store(list) {
  try { localStorage.setItem(STORE, JSON.stringify(list)); return true; } catch { return false; }
}

const HINTS = {
  draw: (el) => (el === 'H'
    ? 'Hydrogen: click an atom to give it one more hydrogen. Drag from atom to atom for a bond; click a bond to make it double or triple.'
    : `Click a hydrogen to put ${el} in its place, or another atom to change it to ${el}. Drag from atom to atom for a bond; click a bond to make it double or triple. Click empty space for a separate atom.`),
  erase: () => 'Click an atom to remove it: its neighbours get hydrogens instead. Click a hydrogen to remove just that hydrogen (a radical is left). Click a bond to break it.',
  move: () => 'Drag an atom: it moves with its hydrogens in the plane of the screen. Drag empty space to turn the view.',
};

export function createEditMode(app) {
  const $ = (id) => document.getElementById(id);
  const { state, canvas } = app;
  const bar = $('edit-bar'), svg = $('skeleton');
  const edit = { active: false, cloud: true }; // cloud: whether the electron cloud is drawn while editing
  let mol = null, history = [], future = [], tool = 'draw', element = 'C';
  let origin = null; // the key of the user's molecule being edited, or null for a new one
  let drag = null, hover = { atom: -1, bond: -1 }, busy = null, relaxed = false, sessions = 0, armed = false;

  const at = (k) => [mol.positions[3 * k], mol.positions[3 * k + 1], mol.positions[3 * k + 2]];
  const sceneNow = () => looseScene(mol.elements, mol.positions, mol.charge);
  const hint = () => app.setStatus(HINTS[tool](element));

  function refreshBar() {
    for (const b of $('edit-tool').children) b.setAttribute('aria-pressed', String(b.dataset.v === tool));
    for (const b of $('edit-element').children) b.setAttribute('aria-pressed', String(b.dataset.v === element));
    $('edit-undo').disabled = !history.length || !!busy;
    $('edit-redo').disabled = !future.length || !!busy;
    $('edit-relax').textContent = busy ? 'Stop' : 'Relax';
    $('edit-relax').setAttribute('aria-pressed', String(!!busy));
    for (const id of ['edit-done', 'edit-xyz', 'edit-charge']) $(id).disabled = !!busy;
    $('edit-charge').value = mol.charge;
    $('edit-delete').hidden = !origin;
    $('edit-delete').textContent = armed ? 'Really delete?' : 'Delete';
    canvas.classList.toggle('editing-erase', tool === 'erase');
  }

  // shows the molecule as it is now; the hint returns unless the calculation failed
  async function show(follow = false) {
    refreshBar();
    app.redraw();
    if (await app.show(sceneNow(), { follow })) hint();
  }

  // one editing step: `change` makes the new molecule from the current one
  function apply(change) {
    if (busy) return false;
    let next;
    try { next = change(mol); } catch (e) { app.setStatus(e.message, true); return false; }
    if (next === mol) return false;
    history.push(mol);
    future = [];
    // The molecule is kept centred: the renderer marches rays through a sphere about the
    // origin that has to hold all of it, so one that grows away from the origin costs more
    // and more to draw. The price is a small shift of the picture after a change.
    mol = E.centred(next);
    relaxed = false;
    hover = { atom: -1, bond: -1 };
    show();
    return true;
  }
  function undo() {
    if (busy || !history.length) return;
    future.push(mol); mol = history.pop(); relaxed = false; hover = { atom: -1, bond: -1 }; show();
  }
  function redo() {
    if (busy || !future.length) return;
    history.push(mol); mol = future.pop(); relaxed = false; hover = { atom: -1, bond: -1 }; show();
  }

  // ---- relaxation: the engine takes the nuclei downhill, frame by frame
  function endRelax(job, changed) {
    busy = null;
    if (changed) { history.push(job.before); future = []; }
    refreshBar();
  }
  function relax() {
    if (busy) { // the button reads "Stop" now
      const job = busy;
      app.engine.relax('stop');
      endRelax(job, job.moved);
      hint();
      job.resolve(false);
      return Promise.resolve(false);
    }
    return new Promise((resolve) => {
      const scene = sceneNow();
      busy = { session: ++sessions, before: mol, moved: false, resolve };
      refreshBar();
      app.setStatus('Relaxing: the nuclei go downhill in energy…');
      app.engine.relax('start', { session: busy.session, numbers: scene.numbers, positions: scene.positions, charge: mol.charge, etemp: state.settings.etemp, unpaired: state.unpaired });
    });
  }
  app.engine.onRelax = (f) => {
    const job = busy;
    if (!job || (f.session !== job.session && !f.fatal)) return;
    if (f.positions && f.result) {
      mol = { ...mol, positions: Array.from(f.positions) };
      job.moved = true;
      f.result.ms = 0;
      app.show(sceneNow(), { result: f.result, follow: true });
    }
    if (!f.done) { app.setStatus(`Relaxing… step ${f.steps}, largest force ${f.gmax.toExponential(1)} Eh/bohr.`); return; }
    endRelax(job, job.moved);
    if (f.error) app.setStatus(`The relaxation stopped: the electronic problem cannot be solved for this structure (${f.error}). Undo the last change or move the atoms apart.`, true);
    else if (!f.converged) app.setStatus(`Stopped after ${f.steps} steps, not fully relaxed yet. Press “Relax” again to go on.`);
    else { relaxed = true; app.setStatus(`Relaxed in ${f.steps} steps.`); }
    job.resolve(!f.error);
  };

  // ---- into the editor and out of it
  function enter(fresh) {
    if (edit.active) return;
    const scene = state.molecule;
    if (!fresh && scene.elements.length > E.MAX_ATOMS) { app.setStatus(`The editor takes molecules of up to ${E.MAX_ATOMS} atoms.`, true); return; }
    const own = !fresh && scene.nFragments === 1 ? app.library().find((m) => m.key === scene.key && m.custom) : null;
    origin = own?.key ?? null;
    mol = fresh ? E.placeAtom(E.emptyMolecule(), 'C', [0, 0, 0]) : E.withBonds(scene.elements, scene.positions, scene.charge);
    history = []; future = []; tool = 'draw'; element = 'C'; relaxed = false; armed = false; drag = null; busy = null;
    hover = { atom: -1, bond: -1 };
    $('edit-name').value = own?.name ?? '';
    $('edit-cloud').checked = edit.cloud;
    edit.active = true;
    app.enter();
    bar.hidden = false;
    svg.style.display = 'block';
    show();
  }
  function leave(key) {
    if (busy) { app.engine.relax('stop'); busy = null; }
    edit.active = false;
    bar.hidden = true;
    svg.style.display = '';
    svg.textContent = '';
    canvas.classList.remove('editing-erase');
    return app.leave(key);
  }

  async function done() {
    if (busy) return;
    const renamed = origin && $('edit-name').value.trim() !== app.library().find((m) => m.key === origin).name;
    if (!history.length && !renamed && state.result) return leave(null); // nothing was changed
    if (!relaxed && !(await relax())) return;
    if (!state.result) { app.setStatus('This structure cannot be calculated, so there is nothing to save yet.', true); return; }
    mol = E.centred(mol);
    const scene = sceneNow();
    const entry = {
      key: origin ?? `my-${Date.now().toString(36)}`, name: $('edit-name').value.trim() || app.name(scene), group: MY_GROUP, custom: true,
      charge: mol.charge, unpaired: state.unpaired !== scene.unpaired ? state.unpaired : undefined,
      elements: mol.elements, xyz: mol.positions.map((v) => +(v * BOHR_ANGSTROM).toFixed(6)), optimized: true,
    };
    const list = savedMolecules().filter((m) => m.key !== entry.key);
    list.push(entry);
    const kept = store(list);
    app.addMolecule(entry);
    await leave(entry.key);
    if (!kept) app.setStatus('The molecule is in the list, but the browser did not let it be stored: it will be gone after a reload.', true);
  }

  function remove() {
    if (!origin || busy) return;
    if (!armed) { armed = true; refreshBar(); return; }
    store(savedMolecules().filter((m) => m.key !== origin));
    app.removeMolecule(origin);
    leave(null);
  }

  // ---- the skeleton: bonds as lines, atoms as discs, far ones first
  function draw(camera) {
    const w = canvas.clientWidth, h = canvas.clientHeight, ox = canvas.offsetLeft, oy = canvas.offsetTop;
    const pts = mol.elements.map((_, k) => project(at(k), camera, w, h)), parts = [];
    mol.bonds.forEach(([i, j, order], b) => {
      const [x1, y1] = pts[i], [x2, y2] = pts[j], l = Math.hypot(x2 - x1, y2 - y1) || 1, nx = -(y2 - y1) / l, ny = (x2 - x1) / l;
      for (const off of order === 1 ? [0] : order === 2 ? [-2.5, 2.5] : [-4.5, 0, 4.5])
        parts.push(`<line class="${hover.bond === b ? 'hot' : ''}" x1="${(ox + x1 + off * nx).toFixed(1)}" y1="${(oy + y1 + off * ny).toFixed(1)}" x2="${(ox + x2 + off * nx).toFixed(1)}" y2="${(oy + y2 + off * ny).toFixed(1)}"/>`);
    });
    if (drag?.mode === 'bond' && drag.cursor) {
      const [x1, y1] = pts[drag.atom], [x2, y2] = drag.to >= 0 ? pts[drag.to] : drag.cursor;
      parts.push(`<line class="rubber" x1="${ox + x1}" y1="${oy + y1}" x2="${ox + x2}" y2="${oy + y2}"/>`);
    }
    const order = pts.map((_, k) => k).sort((a, b) => pts[b][2] - pts[a][2]);
    for (const k of order) {
      const el = mol.elements[k], r = (el === 'H' ? 4.5 : 7.5) * camera.distance / pts[k][2];
      const hot = hover.atom === k || drag?.to === k || (drag?.mode && drag.mode !== 'view' && drag.atom === k);
      parts.push(`<circle class="${hot ? 'hot' : ''}" cx="${(ox + pts[k][0]).toFixed(1)}" cy="${(oy + pts[k][1]).toFixed(1)}" r="${r.toFixed(1)}" fill="${COLOURS[el] ?? '#c99be0'}"/>`);
    }
    svg.innerHTML = parts.join('');
  }

  // ---- the mouse
  const local = (e) => { const r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
  const projected = () => mol.elements.map((_, k) => project(at(k), app.camera(), canvas.clientWidth, canvas.clientHeight));
  // the atom under the cursor: within 13 px, the nearer to the viewer of those almost equally close
  function hitAtom(x, y, except = -1) {
    let best = -1, bestD = 13, bestZ = Infinity;
    projected().forEach(([px, py, z], k) => {
      const d = Math.hypot(px - x, py - y);
      if (k === except || d >= 13) return;
      if (d < bestD - 4 || (d < bestD + 4 && z < bestZ)) { best = k; bestD = d; bestZ = z; }
    });
    return best;
  }
  // the bond under the cursor: within 6 px of its middle part
  function hitBond(x, y) {
    const pts = projected();
    let best = -1, bestD = 6;
    mol.bonds.forEach(([i, j], b) => {
      const [x1, y1] = pts[i], [x2, y2] = pts[j], l2 = (x2 - x1) ** 2 + (y2 - y1) ** 2 || 1;
      const t = ((x - x1) * (x2 - x1) + (y - y1) * (y2 - y1)) / l2;
      if (t < 0.2 || t > 0.8) return;
      const d = Math.hypot(x - x1 - t * (x2 - x1), y - y1 - t * (y2 - y1));
      if (d < bestD) { bestD = d; best = b; }
    });
    return best;
  }
  const target = (x, y) => { const atom = hitAtom(x, y); return { atom, bond: atom < 0 ? hitBond(x, y) : -1 }; };

  function clickAtom(k) {
    const was = mol.elements[k];
    if (tool === 'erase') return apply((m) => E.erase(m, k));
    if (tool !== 'draw') return false;
    if (was === 'H') {
      if (element === 'H') return false;
      return apply((m) => (E.neighbours(m, k).length === 1 ? E.replaceHydrogen(m, k, element) : E.setElement(m, k, element)));
    }
    if (element === 'H') return apply((m) => E.addHydrogen(m, k));
    if (was === element) { app.setStatus(`This atom is ${element} already. Click one of its hydrogens to grow the molecule from it.`); return false; }
    return apply((m) => E.setElement(m, k, element));
  }
  function clickBond(b) {
    if (tool === 'draw') return apply((m) => E.cycleBond(m, b));
    if (tool === 'erase') return apply((m) => E.removeBond(m, b));
    return false;
  }

  edit.pointerdown = (e) => {
    const [x, y] = local(e), hit = busy ? { atom: -1, bond: -1 } : target(x, y);
    drag = { x: e.clientX, y: e.clientY, moved: 0, mode: 'view', ...hit, to: -1 };
    try { canvas.setPointerCapture(e.pointerId); } catch {} // synthetic events in the checks have no real pointer
    if (hit.atom >= 0 && tool === 'move') {
      drag.mode = 'move';
      drag.before = mol;
      drag.depth = project(at(hit.atom), app.camera(), canvas.clientWidth, canvas.clientHeight)[2];
    } else if (hit.atom >= 0 && tool === 'draw') drag.mode = 'bond';
    app.redraw({ cloud: false });
  };
  edit.pointermove = (e) => {
    const [x, y] = local(e);
    if (!drag) {
      const now = busy ? { atom: -1, bond: -1 } : target(x, y);
      if (now.atom !== hover.atom || now.bond !== hover.bond) { hover = now; app.redraw({ cloud: false }); }
      canvas.style.cursor = now.atom >= 0 || now.bond >= 0 ? (tool === 'move' ? (now.atom >= 0 ? 'move' : '') : 'pointer') : '';
      return;
    }
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.x = e.clientX; drag.y = e.clientY; drag.moved += Math.abs(dx) + Math.abs(dy);
    if (drag.mode === 'move') {
      const cam = cameraBasis(app.camera()), perPixel = drag.depth / (2.3 * canvas.clientHeight / 2);
      mol = E.moveAtom(mol, drag.atom, [0, 1, 2].map((c) => (cam.rt[c] * dx - cam.up[c] * dy) * perPixel));
      relaxed = false;
      app.show(sceneNow(), { follow: true });
    } else if (drag.mode === 'bond') {
      drag.cursor = [x + canvas.offsetLeft, y + canvas.offsetTop];
      drag.to = drag.moved > 4 ? hitAtom(x, y, drag.atom) : -1;
    } else {
      app.turnView(dx, dy);
    }
    app.redraw({ cloud: drag.mode === 'view' });
  };
  edit.pointerup = (e) => {
    const d = drag;
    drag = null;
    canvas.style.cursor = '';
    app.redraw({ cloud: false });
    if (!d || e.type !== 'pointerup') { if (d?.mode === 'move' && d.before !== mol) { history.push(d.before); future = []; refreshBar(); } return; }
    if (d.mode === 'move') {
      if (d.before !== mol) { history.push(d.before); future = []; mol = E.centred(mol); show(true); }
    } else if (d.mode === 'bond' && d.to >= 0) apply((m) => E.connect(m, d.atom, d.to));
    else if (d.moved < 4) {
      if (d.atom >= 0) clickAtom(d.atom);
      else if (d.bond >= 0) clickBond(d.bond);
      else if (tool === 'draw' && !busy) {
        const [x, y] = local(e);
        apply((m) => E.placeAtom(m, element, app.unproject(x, y, app.camera().distance)));
      }
    }
  };
  edit.draw = draw;

  // ---- the toolbar
  for (const el of PALETTE) {
    const b = document.createElement('button');
    b.type = 'button'; b.dataset.v = el; b.textContent = el;
    b.style.setProperty('--dot', COLOURS[el]);
    $('edit-element').append(b);
  }
  for (const b of $('edit-tool').children) b.addEventListener('click', () => { tool = b.dataset.v; refreshBar(); hint(); });
  for (const b of $('edit-element').children) b.addEventListener('click', () => { element = b.dataset.v; tool = 'draw'; refreshBar(); hint(); });
  $('edit-cloud').addEventListener('change', (e) => { edit.cloud = e.target.checked; app.redraw(); });
  $('edit-undo').addEventListener('click', undo);
  $('edit-redo').addEventListener('click', redo);
  $('edit-relax').addEventListener('click', relax);
  $('edit-done').addEventListener('click', done);
  $('edit-cancel').addEventListener('click', () => leave(null));
  $('edit-delete').addEventListener('click', remove);
  $('edit-open').addEventListener('click', () => enter(false));
  $('edit-new').addEventListener('click', () => enter(true));
  $('edit-charge').addEventListener('change', (e) => {
    const charge = Math.max(-4, Math.min(4, Math.round(+e.target.value) || 0));
    if (charge === mol.charge || !apply((m) => ({ ...m, charge }))) refreshBar();
  });
  window.addEventListener('keydown', (e) => {
    if (!edit.active || !(e.ctrlKey || e.metaKey) || e.target.matches?.('input, textarea, select')) return;
    const key = e.key.toLowerCase();
    if (key === 'z' && !e.shiftKey) { e.preventDefault(); undo(); }
    else if (key === 'y' || (key === 'z' && e.shiftKey)) { e.preventDefault(); redo(); }
  });
  // the status line sits under the toolbar, however many rows the toolbar takes
  new ResizeObserver(() => app.barHeight(bar.hidden ? 0 : bar.offsetHeight)).observe(bar);

  // XYZ in and out through one window: the text of the molecule, which can be replaced
  const dialog = $('xyz'), text = $('xyz-text'), problem = $('xyz-problem');
  $('edit-xyz').addEventListener('click', () => {
    text.value = E.toXYZ(mol, $('edit-name').value.trim() || app.name(sceneNow()));
    problem.textContent = '';
    dialog.showModal();
  });
  $('xyz-copy').addEventListener('click', async () => {
    try { await navigator.clipboard.writeText(text.value); problem.textContent = 'Copied.'; } catch { text.select(); problem.textContent = 'Select the text and copy it by hand.'; }
  });
  $('xyz-close').addEventListener('click', () => dialog.close());
  $('xyz-apply').addEventListener('click', () => {
    try {
      const next = E.fromXYZ(text.value, mol.charge);
      dialog.close();
      apply(() => next);
    } catch (e) { problem.textContent = e.message; }
  });

  // for the checks
  Object.assign(edit, {
    enter, leave, done, relax, undo, redo,
    molecule: () => mol,
    tool: (t, el) => { tool = t; if (el) element = el; refreshBar(); },
    busy: () => !!busy,
  });
  return edit;
}
