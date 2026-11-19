// GFN2-xTB engine (axtb compiled to WebAssembly) behind the common engine interface.
import { loadXtb, unpack } from './xtb-binding.js';

export async function createXtbEngine(wasmBytes) {
  const xtb = await loadXtb(wasmBytes);
  let calc = null, key = '';

  return {
    name: 'xtb',
    hasForces: true,
    // req: { numbers, positions (bohr, 3n), charge, etemp (K), forces (bool), orbitals (bool),
    // accuracy (SCF accuracy multiplier, default 1), unpaired (number of unpaired electrons;
    // omitted: the lowest the electron count allows) }. With orbitals: false the result has only
    // the energy, the forces and the iteration count.
    compute(req) {
      const k = `${req.numbers.join(',')}|${req.charge}|${req.unpaired ?? -1}`;
      if (k !== key) {
        calc?.destroy();
        calc = null; key = '';
        calc = xtb.create(req.numbers, req.positions, req.charge, req.etemp, req.unpaired ?? -1);
        key = k;
      } else {
        calc.setEtemp(req.etemp);
      }
      calc.setAccuracy(req.accuracy ?? 1);
      let r = calc.singlepoint(req.positions, req.forces !== false, req.orbitals !== false);
      if (r.status === 1 || r.status === 2) {
        // The iterations started from the previous geometry's electrons and did not settle
        // (or blew up numerically). One more attempt from a fresh start; if that converges,
        // the result is as good as any other.
        calc.resetGuess();
        r = calc.singlepoint(req.positions, req.forces !== false, req.orbitals !== false);
      }
      if (r.status) {
        if (r.status === 3) { calc = null; key = ''; } // the module trapped: do not reuse the instance state
        // 'numerical failure' is the iterations diverging: the same thing as not converging
        return { ok: false, converged: r.status !== 1 && !/numerical failure|non-finite/.test(r.error), error: r.error };
      }
      const u = unpack(r.raw);
      if (req.orbitals === false)
        return { ok: true, converged: true, engine: 'xtb', energy: u.energy, forces: u.hasGrad ? Float64Array.from(u.gradient, (g) => -g) : null, iterations: u.iterations };
      // axtb orders p functions y, z, x; the interface uses x, y, z. The d functions keep
      // axtb's order (m = -2..2: xy, yz, z2, xz, x2-y2; see orbitals.js). GFN2-xTB has no f.
      const perm = [];
      for (const sh of u.shells) {
        const b = perm.length;
        if (sh.l === 0) perm.push(b);
        else if (sh.l === 1) perm.push(b + 2, b, b + 1);
        else if (sh.l === 2) perm.push(b, b + 1, b + 2, b + 3, b + 4);
        else return { ok: false, converged: true, error: `shell with l = ${sh.l} is not supported` };
      }
      const coeff = new Float64Array(u.nmo * u.nao);
      for (let m = 0; m < u.nmo; m++)
        for (let i = 0; i < u.nao; i++) coeff[m * u.nao + i] = u.coeff[m * u.nao + perm[i]];
      return {
        ok: true,
        converged: true,
        engine: 'xtb',
        energy: u.energy,
        forces: u.hasGrad ? Float64Array.from(u.gradient, (g) => -g) : null,
        eps: Float64Array.from(u.eps),
        occ: Float64Array.from(u.occ), // both spins together, 0..2
        occAlpha: Float64Array.from(u.occAlpha), // spin up alone, 0..1
        coeff,
        nao: u.nao,
        nmo: u.nmo,
        shells: u.shells.map((s) => ({ atom: s.atom, l: s.l, alpha: Float64Array.from(s.alpha), coeff: Float64Array.from(s.coeff) })),
        iterations: u.iterations,
      };
    },
  };
}
