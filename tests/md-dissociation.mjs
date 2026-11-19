// How hot must BH3NH3 be to fall apart at the B-N bond within a watchable time?
//   node tests/md-dissociation.mjs
// "Broken" = B-N beyond 3.0 A (equilibrium 1.66 A) and no B-H or N-H bond broken first.
import { setup } from './md-explore.mjs';
import { BOHR_ANGSTROM } from '../app/src/units.js';
const d = (x, i, j) => Math.hypot(x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]) * BOHR_ANGSTROM;
const LIMIT = 20; // ps
for (const T of [1000, 1500, 2000, 2500, 3000, 4000]) {
  const times = [], other = [];
  let failed = 0;
  for (let seed = 1; seed <= 8; seed++) {
    const { md } = setup('bh3nh3', { accuracy: 0.1, seed, wallRadius: 6.5 });
    md.params.temperature = T;
    md.init(false); md.thermalize(T);
    let when = null;
    while (md.time < LIMIT * 1000) {
      if (!md.step()) { failed++; break; }
      if (d(md.x, 0, 4) > 3.0) { when = md.time / 1000; break; }
      // B is atom 0 with H 1-3, N is atom 4 with H 5-7
      if ([1, 2, 3].some((h) => d(md.x, 0, h) > 2.2) || [5, 6, 7].some((h) => d(md.x, 4, h) > 2.0)) { other.push(md.time / 1000); break; }
    }
    if (when !== null) times.push(when);
  }
  times.sort((a, b) => a - b);
  console.log(`${String(T).padStart(5)} K: B-N broke in ${times.length}/8 runs within ${LIMIT} ps` +
    (times.length ? `, times ${times.map((t) => t.toFixed(2)).join(' ')} ps, median ${times[times.length >> 1].toFixed(2)} ps` : '') +
    (other.length ? `; another bond broke first in ${other.length} (at ${other.map((t) => t.toFixed(2)).join(' ')} ps)` : '') + (failed ? `; SCF failed in ${failed}` : ''));
}
