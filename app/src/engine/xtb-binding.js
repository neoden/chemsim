// Low-level binding for xtb_wasm.wasm (flat C ABI, see engine/xtb-wasm/src/lib.rs).
// Works in node, a page and a worker.

// The module is built for wasm32-wasip1 so that std::time works (axtb times its SCF phases).
// It needs six WASI calls; this is all of them.
function wasiShim(mem) {
  const view = () => new DataView(mem().buffer);
  return {
    clock_time_get(_id, _precision, ptr) {
      view().setBigUint64(ptr, BigInt(Math.round(performance.now() * 1e6)), true);
      return 0;
    },
    environ_sizes_get(countPtr, sizePtr) {
      view().setUint32(countPtr, 0, true);
      view().setUint32(sizePtr, 0, true);
      return 0;
    },
    environ_get: () => 0,
    fd_write(_fd, iovs, n, nwrittenPtr) {
      const v = view();
      let total = 0, text = '';
      for (let i = 0; i < n; i++) {
        const ptr = v.getUint32(iovs + 8 * i, true), len = v.getUint32(iovs + 8 * i + 4, true);
        text += new TextDecoder().decode(new Uint8Array(mem().buffer, ptr, len).slice());
        total += len;
      }
      if (text.trim()) console.warn('[xtb]', text.trimEnd());
      v.setUint32(nwrittenPtr, total, true);
      return 0;
    },
    proc_exit(code) { throw new Error(`xtb module exited with code ${code}`); },
    sched_yield: () => 0,
  };
}

export async function loadXtb(bytes) {
  let x;
  const { instance } = await WebAssembly.instantiate(bytes, { wasi_snapshot_preview1: wasiShim(() => x.memory) });
  x = instance.exports;
  const f64 = (ptr, len) => new Float64Array(x.memory.buffer, ptr, len);
  const lastError = () =>
    new TextDecoder().decode(new Uint8Array(x.memory.buffer, x.xtb_err_ptr(), x.xtb_err_len()).slice());

  // numbers: atomic numbers; pos: 3*nat in bohr; etemp: electronic temperature in kelvin;
  // unpaired: number of unpaired electrons, -1 for the lowest the electron count allows
  function create(numbers, pos, charge = 0, etemp = 300, unpaired = -1) {
    const nat = numbers.length;
    const pn = x.xtb_alloc(nat), pp = x.xtb_alloc(3 * nat);
    f64(pn, nat).set(numbers);
    f64(pp, 3 * nat).set(pos);
    const h = x.xtb_create(pn, nat, pp, charge, etemp, unpaired);
    x.xtb_free(pn, nat);
    if (!h) { x.xtb_free(pp, 3 * nat); throw new Error(lastError()); }
    return {
      // status: 0 ok, 1 SCF not converged, 2 other error, 3 module trapped.
      // raw: Float64Array copy of the result buffer (see unpack).
      singlepoint(pos, gradient = true, orbitals = true) {
        f64(pp, 3 * nat).set(pos);
        let status;
        try {
          status = x.xtb_singlepoint(h, pp, gradient ? 1 : 0, orbitals ? 1 : 0);
        } catch (e) {
          return { status: 3, error: `${e.message}; ${lastError()}` };
        }
        if (status) return { status, error: lastError() };
        return { status, raw: f64(x.xtb_out_ptr(h), x.xtb_out_len(h)).slice() };
      },
      setEtemp: (t) => x.xtb_set_etemp(h, t),
      // SCF accuracy multiplier: 1 is axtb's default, smaller is tighter
      setAccuracy: (a) => x.xtb_set_accuracy(h, a),
      resetGuess: () => x.xtb_reset_guess(h),
      destroy() { x.xtb_destroy(h); x.xtb_free(pp, 3 * nat); },
    };
  }
  return { create };
}

// Unpack the flat result buffer. AO order here is axtb's: p functions as y, z, x.
export function unpack(raw) {
  const [, energy, nat, nao, nmo, iterations, nshell, hasGrad] = raw;
  let o = 8;
  const take = (n) => raw.subarray(o, (o += n));
  const gradient = take(3 * nat), charges = take(nat), dipoles = take(3 * nat), quadrupoles = take(6 * nat);
  const eps = take(nmo), occ = take(nmo), occAlpha = take(nmo);
  const coeff = take(nmo * nao);
  const shells = [];
  for (let s = 0; s < nshell; s++) {
    const [atom, l, nprim] = take(3);
    shells.push({ atom, l, alpha: take(nprim), coeff: take(nprim) });
  }
  return { energy, nat, nao, nmo, iterations, hasGrad: !!hasGrad, gradient, charges, dipoles, quadrupoles, eps, occ, occAlpha, coeff, shells };
}
