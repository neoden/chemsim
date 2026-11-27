// Runs the quantum engine (GFN2-xTB) and the molecular dynamics off the main thread.
//
// Engine interface. Request: { numbers, positions (bohr, flat 3n),
// charge, etemp (K), forces, unpaired (number of unpaired electrons) }. Result: { ok, converged, error?, energy (Eh), forces (Eh/bohr)
// or null, eps (Eh), occ, coeff (nmo x nao, orbital by orbital), nao, nmo, shells, iterations }.
// A shell is { atom, l, alpha[], coeff[] }: the radial part is sum_k coeff_k exp(-alpha_k r^2),
// multiplied by 1 for l = 0, by x, y, z (in this AO order) for l = 1, and by the five
// quadratic polynomials listed in orbitals.js for l = 2.
//
// Messages in:
//   { type: 'compute', id, ...request }   one calculation; if several pile up, only the newest runs
//   { type: 'md', cmd: 'init', numbers, positions, charge, etemp, unpaired, params, speed, wallRadius, seed, thermalize }
//   { type: 'md', cmd: 'run' | 'pause' | 'step' | 'stop' }
//   { type: 'md', cmd: 'set', params?, etemp?, speed?, unpaired?, grab? }   grab: { atom, target | toAtom + rest, piece? } or null, see md.js
//   { type: 'md', cmd: 'move', positions, wallRadius }           nuclei placed by hand
//   { type: 'relax', cmd: 'start', session, numbers, positions, charge, etemp, unpaired }   geometry optimization
//   { type: 'relax', cmd: 'stop' }
// Messages out: compute results (with the request's id), { type: 'md-frame', ... } and
// { type: 'relax-frame', session, done, converged, error?, positions, energy, gmax, steps, result }.
import { createXtbEngine } from './xtb.js';
import { Dynamics } from '../md.js';
import { minimize, modelGuess } from '../relax.js';

let xtb = null, xtbLoading = null;
async function engine() {
  if (!xtb) {
    xtbLoading ??= fetch(new URL('../../engine/xtb_wasm.wasm', import.meta.url))
      .then((r) => { if (!r.ok) throw new Error(`xtb_wasm.wasm: HTTP ${r.status}`); return r.arrayBuffer(); })
      .then(createXtbEngine);
    xtb = await xtbLoading;
  }
  return xtb;
}
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

// ---- single calculations
let pending = null, running = false;
async function drain() {
  if (running) return;
  running = true;
  while (pending) {
    const req = pending;
    pending = null;
    let res;
    const t0 = performance.now();
    try {
      res = (await engine()).compute(req);
    } catch (e) {
      res = { ok: false, converged: true, error: String(e?.message ?? e) };
    }
    res.id = req.id;
    res.ms = performance.now() - t0;
    postMessage(res);
    await sleep(0); // let queued messages replace `pending`
  }
  running = false;
}

// ---- molecular dynamics
// SCF accuracy during dynamics: 10 times tighter than axtb's default. Measured: the energy
// drift without thermostat falls from ~4e-6 to ~2e-7 Eh/atom/ps for ~15 % more time per step.
const MD_ACCURACY = 0.1;
const TICK = 16; // ms between frames sent to the page
let md = null; // { dyn, speed (fs of simulation per second, 0 = as fast as possible), running, etemp, session }
let looping = false;

function frame(extra = {}) {
  const d = md.dyn;
  postMessage({
    type: 'md-frame', session: md.session, running: md.running, time: d.time, steps: d.steps,
    epot: d.epot, ekin: d.ekin, etot: d.etot, temperature: d.temperature,
    positions: d.x.slice(), result: d.result, rate: md.rate, ...extra,
  });
}

function fail() {
  md.running = false;
  const e = md.dyn.error;
  frame({ error: e.error, converged: e.converged !== false });
}

async function mdInit(msg) {
  const xtb = await engine();
  const session = msg.session;
  const state = { speed: msg.speed, running: false, etemp: msg.etemp, unpaired: msg.unpaired, session, rate: 0 };
  state.dyn = new Dynamics({
    numbers: msg.numbers, positions: msg.positions, seed: msg.seed, wallRadius: msg.wallRadius,
    force: (x, orbitals) => xtb.compute({ numbers: msg.numbers, positions: x, charge: msg.charge, etemp: state.etemp, unpaired: state.unpaired, forces: true, orbitals, accuracy: MD_ACCURACY }),
  });
  Object.assign(state.dyn.params, msg.params);
  md = state;
  if (!state.dyn.init(true)) return fail();
  if (msg.thermalize) state.dyn.thermalize(state.dyn.params.temperature);
  frame();
}

async function mdLoop() {
  if (looping) return;
  looping = true;
  let carry = 0, last = performance.now(), count = 0, since = last;
  while (md?.running) {
    const start = performance.now();
    // how many steps this tick: the simulation advances `speed` fs per second of wall time
    carry = md.speed > 0 ? Math.min(carry + md.speed / md.dyn.params.dt * Math.min(0.1, (start - last) / 1000), 1e4) : Infinity;
    last = start;
    const want = Math.floor(carry);
    let done = 0, failed = false;
    while (done < want) {
      // the last step of the tick also returns the orbitals, for drawing
      const lastOne = done === want - 1 || performance.now() - start > TICK - 4;
      if (!md.dyn.step(lastOne)) { failed = true; break; }
      done++;
      if (lastOne) break;
    }
    carry = md.speed > 0 ? carry - done : 0;
    count += done;
    if (start - since > 500) { md.rate = count / (performance.now() - since) * 1000; count = 0; since = performance.now(); }
    if (failed) { fail(); break; }
    if (done) frame();
    await sleep(Math.max(0, TICK - (performance.now() - start)));
  }
  looping = false;
}

let initializing = Promise.resolve();
async function onMd(msg) {
  if (msg.cmd === 'init') { initializing = mdInit(msg); return initializing; }
  await initializing; // commands sent right after 'init' wait for it
  if (!md || !md.dyn.ready || msg.session !== md.session) return;
  if (msg.cmd === 'run') { md.running = true; mdLoop(); }
  else if (msg.cmd === 'pause') { md.running = false; frame(); }
  else if (msg.cmd === 'step') { md.running = false; if (md.dyn.step(true)) frame(); else fail(); }
  else if (msg.cmd === 'stop') { md.running = false; md = null; }
  else if (msg.cmd === 'set') {
    if (msg.params) Object.assign(md.dyn.params, msg.params);
    if (msg.etemp !== undefined) md.etemp = msg.etemp;
    if (msg.speed !== undefined) md.speed = msg.speed;
    if (msg.grab !== undefined) md.dyn.grab = msg.grab;
    if (msg.unpaired !== undefined) md.unpaired = msg.unpaired;
  } else if (msg.cmd === 'move') {
    md.dyn.params.wallRadius = msg.wallRadius;
    if (md.dyn.moveTo(msg.positions, true)) frame(); else fail();
  }
}

// ---- geometry optimization: the nuclei go downhill to the nearest minimum
// The largest force left, Eh/bohr: a little looser than the library build (2e-5).
const RELAX_GTOL = 5e-5, RELAX_STEPS = 600;
let relaxing = null;
async function onRelax(msg) {
  if (msg.cmd !== 'start') { relaxing = null; return; }
  const job = relaxing = { session: msg.session };
  const xtb = await engine();
  let error = null;
  const walk = minimize((x) => {
    const r = xtb.compute({ numbers: msg.numbers, positions: x, charge: msg.charge, etemp: msg.etemp, unpaired: msg.unpaired, forces: true });
    if (!r.ok) { error = r.error; return { e: Infinity }; }
    return { e: r.energy, g: r.forces.map((f) => -f), r };
  }, msg.positions, { gtol: RELAX_GTOL, maxIter: RELAX_STEPS, guess: (x) => modelGuess(msg.numbers, x) });
  let sent = performance.now();
  for (;;) {
    const step = walk.next(), v = step.value, now = performance.now();
    if (step.done || now - sent > TICK) {
      postMessage({ type: 'relax-frame', session: job.session, done: step.done, converged: !!v.converged, error: v.failed ? error : undefined,
        positions: v.x, energy: v.e, gmax: v.gmax, steps: v.it, result: v.point.r });
      if (step.done) break;
      await sleep(0); // lets a 'stop' through
      if (relaxing !== job) return;
      sent = performance.now();
    }
  }
  if (relaxing === job) relaxing = null;
}

onmessage = (e) => {
  if (e.data.type === 'relax') onRelax(e.data).catch((err) => postMessage({ type: 'relax-frame', session: e.data.session, done: true, converged: false, error: String(err?.message ?? err) }));
  else if (e.data.type === 'md') onMd(e.data).catch((err) => postMessage({ type: 'md-frame', session: e.data.session, running: false, error: String(err?.message ?? err), converged: true }));
  else { pending = e.data; drain(); }
};
