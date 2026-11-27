// Geometry optimization: BFGS in Cartesian coordinates with a step cap and backtracking.
// No DOM and no engine of its own, so it serves the Worker and the library build alike.
//
// evaluate(x) returns { e, g } for positions x (bohr, flat 3n): the energy and its gradient.
// It may return anything else next to them; a point that cannot be calculated is { e: Infinity }.
//
// A generator: it yields { x, e, gmax, it, point } before every iteration (point is what
// evaluate returned for x) and returns the same with `converged`, and with `failed` when the
// walk ran into points that cannot be calculated.
//
// guess(x): the inverse Hessian to start from at x, as rows (see modelGuess); without it every
// coordinate is taken to be as stiff as a bond, which makes a floppy molecule creep. The guess
// is taken afresh every `refresh` steps: over a long way what was learned about the stiffness
// at the start no longer holds, and the model at the new place is the better beginning.
// Measured from the editor's geometries, without a guess and with it: C32H66 59 steps and 9,
// cholesterol 348 and about 90, oxytocin (135 atoms) not there after 600 and about 300.
const SOFTEST = 2e-3; // Eh/bohr^2, see below
export function* minimize(evaluate, x0, { gtol = 2e-5, maxIter = 800, guess = null, refresh = 40 } = {}) {
  const n = x0.length;
  let x = Float64Array.from(x0), point = evaluate(x), { e, g } = point;
  if (!Number.isFinite(e)) return { x, e, gmax: Infinity, it: 0, point, converged: false, failed: true };
  let H = guess ? guess(x) : Array.from({ length: n }, (_, i) => { const row = new Float64Array(n); row[i] = 2; return row; }); // inverse Hessian
  for (let it = 0; ; it++) {
    if (guess && it && it % refresh === 0) H = guess(x);
    const gmax = g.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    if (gmax < gtol) return { x, e, gmax, it, point, converged: true };
    if (it >= maxIter) return { x, e, gmax, it, point, converged: false };
    yield { x, e, gmax, it, point };
    let d = new Float64Array(n);
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) d[i] -= H[i][j] * g[j];
    if (d.reduce((s, v, i) => s + v * g[i], 0) >= 0) d = g.map((v) => -v); // not a descent direction
    const dmax = d.reduce((m, v) => Math.max(m, Math.abs(v)), 0);
    let scale = Math.min(1, 0.3 / dmax), next, xn;
    for (;;) {
      xn = x.map((v, i) => v + scale * d[i]);
      next = evaluate(xn);
      if (next.e <= e + 1e-10 || scale < 1e-4) break;
      scale *= 0.4;
    }
    if (!Number.isFinite(next.e)) return { x, e, gmax, it, point, converged: false, failed: true };
    const s = xn.map((v, i) => v - x[i]), y = next.g.map((v, i) => v - g[i]);
    let sy = s.reduce((t, v, i) => t + v * y[i], 0);
    // Along a flat or downward-curving stretch (a group turning over its barrier) the step
    // teaches nothing about the stiffness, or says it is negative. Left alone, the walk keeps
    // whatever stiffness it assumed there and creeps. Such a step is taken to have met the
    // least stiffness worth believing in, SOFTEST: the next one along it is that much longer.
    const ss = s.reduce((t, v) => t + v * v, 0);
    if (sy < SOFTEST * ss) {
      const more = (SOFTEST * ss - sy) / ss;
      for (let i = 0; i < n; i++) y[i] += more * s[i];
      sy = SOFTEST * ss;
    }
    if (sy > 1e-12) {
      const Hy = new Float64Array(n);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) Hy[i] += H[i][j] * y[j];
      const yHy = y.reduce((t, v, i) => t + v * Hy[i], 0);
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++)
        H[i][j] += ((sy + yHy) * s[i] * s[j]) / (sy * sy) - (Hy[i] * s[j] + s[i] * Hy[j]) / sy;
    }
    x = xn; e = next.e; g = next.g; point = next;
  }
}

// ---- a guess at the Hessian from the geometry alone
// The model Hessian of Lindh, Bernhardsson, Karlström and Malmqvist (Chem. Phys. Lett. 241
// (1995) 423): every pair of atoms is a spring along the line between them, every three a
// spring on the angle, every four a spring on the torsion, with stiffness that falls off
// with distance as rho = exp(alpha (r0^2 - r^2)), so that only neighbours count and no list
// of bonds is needed. alpha and r0 go by the periods of the two atoms (1, 2, 3 and beyond).
const ALPHA = [[1.0, 0.3949, 0.3949], [0.3949, 0.28, 0.28], [0.3949, 0.28, 0.28]];
const R0 = [[1.35, 2.10, 2.53], [2.10, 2.87, 3.40], [2.53, 3.40, 3.40]]; // bohr
const K_STRETCH = 0.45, K_BEND = 0.15, K_TORSION = 0.005; // Eh/bohr^2, Eh/rad^2
const RHO_CUT = 1e-3; // pairs, triples and fours weaker than this are left out
const period = (z) => (z <= 2 ? 0 : z <= 10 ? 1 : 2);

// numbers: atomic numbers; x: positions (bohr, flat 3n). Returns the 3n x 3n matrix, flat.
export function modelHessian(numbers, x) {
  const nat = numbers.length, n = 3 * nat, H = new Float64Array(n * n);
  const rho = new Float64Array(nat * nat), near = Array.from({ length: nat }, () => []);
  const vec = (i, j) => [x[3 * i] - x[3 * j], x[3 * i + 1] - x[3 * j + 1], x[3 * i + 2] - x[3 * j + 2]];
  const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
  for (let i = 0; i < nat; i++) for (let j = 0; j < i; j++) {
    const d = vec(i, j), a = period(numbers[i]), b = period(numbers[j]);
    const w = Math.exp(ALPHA[a][b] * (R0[a][b] * R0[a][b] - dot(d, d)));
    if (w < RHO_CUT) continue;
    rho[i * nat + j] = rho[j * nat + i] = w;
    near[i].push(j); near[j].push(i);
  }
  // a spring of stiffness k on a coordinate q adds k (dq/dx)(dq/dx)^T; grads: dq/dx of each atom
  const add = (k, atoms, grads) => {
    for (let a = 0; a < atoms.length; a++) for (let b = 0; b < atoms.length; b++)
      for (let p = 0; p < 3; p++) {
        const row = (3 * atoms[a] + p) * n + 3 * atoms[b], g = k * grads[a][p];
        H[row] += g * grads[b][0]; H[row + 1] += g * grads[b][1]; H[row + 2] += g * grads[b][2];
      }
  };
  for (let j = 0; j < nat; j++) for (const k of near[j]) {
    const rjk = rho[j * nat + k];
    if (k < j) { // the distance j-k
      const d = vec(j, k), l = Math.sqrt(dot(d, d)), u = d.map((v) => v / l);
      add(K_STRETCH * rjk, [j, k], [u, u.map((v) => -v)]);
    }
    for (const i of near[j]) { // the angle i-j-k at j
      if (i >= k) continue;
      const w = rho[i * nat + j] * rjk;
      if (w < RHO_CUT) continue;
      const a = vec(i, j), b = vec(k, j), la = Math.sqrt(dot(a, a)), lb = Math.sqrt(dot(b, b));
      const u = a.map((v) => v / la), v = b.map((t) => t / lb), c = dot(u, v), s = Math.sqrt(Math.max(0, 1 - c * c));
      if (s < 0.15) continue; // nearly straight: the angle has no direction to open in
      const gi = u.map((t, p) => (c * t - v[p]) / (la * s)), gk = v.map((t, p) => (c * t - u[p]) / (lb * s));
      add(K_BEND * w, [i, j, k], [gi, gi.map((t, p) => -t - gk[p]), gk]);
    }
    if (k > j) continue;
    for (const i of near[j]) { // the torsion i-j-k-l about j-k
      if (i === k) continue;
      const wi = rho[i * nat + j] * rjk;
      if (wi < RHO_CUT) continue;
      for (const l of near[k]) {
        if (l === j || l === i) continue;
        const w = wi * rho[k * nat + l];
        if (w < RHO_CUT) continue;
        const b1 = vec(j, i), b2 = vec(k, j), b3 = vec(l, k), n1 = cross(b1, b2), n2 = cross(b2, b3);
        const l2 = dot(b2, b2), m1 = dot(n1, n1), m2 = dot(n2, n2);
        // either end nearly in line with the axis: no torsion to speak of
        if (m1 < 0.0225 * dot(b1, b1) * l2 || m2 < 0.0225 * dot(b3, b3) * l2) continue;
        const lb = Math.sqrt(l2), gi = n1.map((t) => -lb / m1 * t), gl = n2.map((t) => lb / m2 * t);
        const f = dot(b1, b2) / l2, h = dot(b3, b2) / l2;
        add(K_TORSION * w, [i, j, k, l], [gi, gi.map((t, p) => h * gl[p] - (1 + f) * t), gl.map((t, p) => f * gi[p] - (1 + h) * t), gl]);
      }
    }
  }
  return H;
}

// The inverse of the model Hessian, as the rows minimize() starts from. floor (Eh/bohr^2) is
// added to every direction: the model gives none to the molecule moving or turning as a whole,
// and next to none to what is held by weak contacts alone.
export function modelGuess(numbers, x, { floor = 0.01 } = {}) {
  const n = x.length, A = modelHessian(numbers, x);
  for (let i = 0; i < n; i++) A[i * n + i] += floor;
  // Cholesky, A = L L^T, L in the lower triangle
  for (let j = 0; j < n; j++) {
    let d = A[j * n + j];
    for (let k = 0; k < j; k++) d -= A[j * n + k] * A[j * n + k];
    d = Math.sqrt(d);
    A[j * n + j] = d;
    for (let i = j + 1; i < n; i++) {
      let v = A[i * n + j];
      for (let k = 0; k < j; k++) v -= A[i * n + k] * A[j * n + k];
      A[i * n + j] = v / d;
    }
  }
  // the inverse, column by column: L L^T z = e_c
  const rows = Array.from({ length: n }, () => new Float64Array(n)), z = new Float64Array(n);
  for (let c = 0; c < n; c++) {
    z.fill(0, 0, c);
    for (let i = c; i < n; i++) {
      let v = i === c ? 1 : 0;
      for (let k = c; k < i; k++) v -= A[i * n + k] * z[k];
      z[i] = v / A[i * n + i];
    }
    for (let i = n - 1; i >= 0; i--) {
      let v = z[i];
      for (let k = i + 1; k < n; k++) v -= A[k * n + i] * z[k];
      z[i] = v / A[i * n + i];
    }
    for (let i = 0; i < n; i++) rows[i][c] = z[i];
  }
  return rows;
}
