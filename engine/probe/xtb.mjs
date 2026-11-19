// Helpers for the stage-0 probes; the binding itself lives in app/src/engine/xtb-binding.js.
export const BOHR = 0.529177210903;
const Z = { H: 1, B: 5, C: 6, N: 7, O: 8 };

export function parseXyz(text) {
  const lines = text.split('\n');
  const nat = parseInt(lines[0]);
  const numbers = [], pos = [];
  for (const line of lines.slice(2, 2 + nat)) {
    const [el, x, y, z] = line.trim().split(/\s+/);
    numbers.push(Z[el]);
    pos.push(+x / BOHR, +y / BOHR, +z / BOHR);
  }
  return { numbers, pos: Float64Array.from(pos) };
}

export { loadXtb, unpack } from '../../app/src/engine/xtb-binding.js';

export function jitter(pos, step) {
  const q = Float64Array.from(pos);
  for (let i = 0; i < q.length / 3; i++)
    for (let k = 0; k < 3; k++) q[3 * i + k] += 0.01 * Math.sin(1.7 * step + 2.3 * i + 0.9 * k);
  return q;
}

// cold run + `reps` warm runs on slightly displaced geometries (same as the native binary)
export function bench(xtb, { numbers, pos }, reps = 20) {
  let t = performance.now();
  const eng = xtb.create(numbers, pos);
  const first = eng.singlepoint(pos);
  const cold = performance.now() - t;
  if (first.status) return { error: first.error };
  const warm = [];
  let last = first;
  const tAll = performance.now();
  for (let step = 1; step <= reps; step++) {
    const p = jitter(pos, step);
    t = performance.now();
    last = eng.singlepoint(p);
    warm.push(performance.now() - t);
    if (last.status) return { error: last.error };
  }
  // browsers round performance.now() to 0.1 ms, so the mean over all runs is the precise figure
  const warm_mean_ms = (performance.now() - tAll) / reps;
  warm.sort((a, b) => a - b);
  eng.destroy();
  return { cold_ms: cold, warm_median_ms: warm[reps >> 1], warm_mean_ms, out: Array.from(first.raw), last: Array.from(last.raw) };
}
