// Stage 3 acceptance checks of the dynamics, in node:  node tests/md.test.mjs
import assert from 'node:assert/strict';
import { setup, slope } from './md-explore.mjs';
import { Dynamics } from '../app/src/md.js';
import { BOHR_ANGSTROM } from '../app/src/units.js';

const dist = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]) * BOHR_ANGSTROM;
const ACC = 0.1; // SCF accuracy used by the app during dynamics
let checks = 0;
const ok = (cond, msg) => { assert.ok(cond, msg); checks++; console.log('  ok  ' + msg); };

console.log('energy conservation without thermostat (dt 0.5 fs, 2 ps, start at 300 K)');
for (const key of ['h2o', 'bh3nh3', 'c6h6']) {
  const { md } = setup(key, { accuracy: ACC });
  md.params.thermostat = false;
  md.init(false); md.thermalize(300);
  const t = [], e = [], t0 = performance.now();
  while (md.time < 2000) { assert.ok(md.step(), 'SCF failed'); t.push(md.time / 1000); e.push(md.etot); }
  const perStep = (performance.now() - t0) / md.steps, drift = slope(t, e) / md.n;
  ok(Math.abs(drift) < 1e-5, `${key}: drift ${drift.toExponential(1)} Eh/atom/ps (limit 1e-5); ${perStep.toFixed(3)} ms/step = ${Math.round(1000 / perStep)} steps/s in node`);
}

console.log('water at 300 K for 10 ps (thermostat on)');
{
  const { md } = setup('h2o', { accuracy: ACC, seed: 3 });
  md.init(false); md.thermalize(300);
  const r0 = [dist(md.x, 0, 1), dist(md.x, 0, 2)];
  let lo = 9, hi = 0, sum = 0, n = 0, tsum = 0;
  while (md.time < 10000) {
    assert.ok(md.step(), 'SCF failed');
    for (const h of [1, 2]) { const r = dist(md.x, 0, h); lo = Math.min(lo, r); hi = Math.max(hi, r); sum += r; n++; }
    tsum += md.temperature;
  }
  ok(hi < 1.25 && lo > 0.8, `O-H stays between ${lo.toFixed(3)} and ${hi.toFixed(3)} A, mean ${(sum / n).toFixed(4)} A (equilibrium ${r0[0].toFixed(4)} A): the molecule is intact`);
  ok(Math.abs(sum / n - r0[0]) < 0.02, 'the bond length oscillates about the equilibrium value');
  const tmean = tsum / md.steps;
  ok(tmean > 240 && tmean < 360, `mean temperature ${tmean.toFixed(0)} K for a target of 300 K (3 atoms: large fluctuations are expected)`);
}

console.log('BH3NH3: intact at 300 K, breaks at the B-N bond when heated');
{
  const cold = setup('bh3nh3', { accuracy: ACC, seed: 1, wallRadius: 6.5 }).md;
  cold.init(false); cold.thermalize(300);
  let hi = 0;
  while (cold.time < 5000) { assert.ok(cold.step(), 'SCF failed'); hi = Math.max(hi, dist(cold.x, 0, 4)); }
  ok(hi < 2.2, `300 K, 5 ps: B-N never exceeds ${hi.toFixed(2)} A`);
  const hot = setup('bh3nh3', { accuracy: ACC, seed: 1, wallRadius: 6.5 }).md;
  hot.params.temperature = 1500;
  hot.init(false); hot.thermalize(1500);
  let when = null;
  while (hot.time < 10000 && when === null) { assert.ok(hot.step(), 'SCF failed'); if (dist(hot.x, 0, 4) > 3.0) when = hot.time; }
  const bh = Math.max(...[1, 2, 3].map((h) => dist(hot.x, 0, h))), nh = Math.max(...[5, 6, 7].map((h) => dist(hot.x, 4, h)));
  ok(when !== null && bh < 1.8 && nh < 1.6, `1500 K: B-N passes 3 A after ${(when / 1000).toFixed(2)} ps while B-H (max ${bh.toFixed(2)} A) and N-H (max ${nh.toFixed(2)} A) hold`);
}

console.log('the soft wall keeps fragments inside');
{
  const { md } = setup('bh3nh3', { accuracy: ACC, seed: 2, wallRadius: 6.5 });
  md.params.temperature = 3000;
  md.init(false); md.thermalize(3000);
  let far = 0;
  while (md.time < 3000) { assert.ok(md.step(), 'SCF failed'); for (let i = 0; i < md.n; i++) far = Math.max(far, Math.hypot(md.x[3 * i], md.x[3 * i + 1], md.x[3 * i + 2])); }
  ok(far < 6.5 + 2.5, `3000 K, 3 ps: no atom gets farther than ${far.toFixed(2)} bohr from the centre (wall at 6.5)`);
}

console.log('a failed electronic step stops the dynamics and leaves the state untouched');
{
  let calls = 0;
  const md = new Dynamics({
    numbers: [1, 1], positions: [0, 0, -0.7, 0, 0, 0.7],
    force: (x) => (++calls > 3 ? { ok: false, converged: false, error: 'SCF not converged' } : { ok: true, energy: 0, forces: Float64Array.from([0, 0, 0.01, 0, 0, -0.01]) }),
  });
  md.params.thermostat = false;
  md.init();
  md.step(); md.step();
  const x = md.x.slice(), v = md.v.slice(), time = md.time;
  ok(md.step() === false && md.error.converged === false, 'step() reports the failure');
  ok(md.x.every((q, i) => q === x[i]) && md.v.every((q, i) => q === v[i]) && md.time === time, 'positions, velocities and time are as before the failed step');
}
console.log(`${checks} checks passed`);
