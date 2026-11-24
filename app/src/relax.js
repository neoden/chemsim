// Geometry optimization: BFGS in Cartesian coordinates with a step cap and backtracking.
// No DOM and no engine of its own, so it serves the Worker and the library build alike.
//
// evaluate(x) returns { e, g } for positions x (bohr, flat 3n): the energy and its gradient.
// It may return anything else next to them; a point that cannot be calculated is { e: Infinity }.
//
// A generator: it yields { x, e, gmax, it, point } before every iteration (point is what
// evaluate returned for x) and returns the same with `converged`, and with `failed` when the
// walk ran into points that cannot be calculated.
export function* minimize(evaluate, x0, { gtol = 2e-5, maxIter = 800 } = {}) {
  const n = x0.length;
  let x = Float64Array.from(x0), point = evaluate(x), { e, g } = point;
  if (!Number.isFinite(e)) return { x, e, gmax: Infinity, it: 0, point, converged: false, failed: true };
  const H = Array.from({ length: n }, (_, i) => { const row = new Float64Array(n); row[i] = 2; return row; }); // inverse Hessian guess
  for (let it = 0; ; it++) {
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
    const sy = s.reduce((t, v, i) => t + v * y[i], 0);
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
