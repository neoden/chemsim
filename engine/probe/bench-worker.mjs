import { loadXtb, parseXyz, bench } from './xtb.mjs';
try {
  const t0 = performance.now();
  const xtb = await loadXtb(await (await fetch('/xtb_wasm.wasm')).arrayBuffer());
  const load_ms = performance.now() - t0;
  const results = {};
  for (const name of ['h2o', 'bh3_nh3', 'hexane']) {
    const mol = parseXyz(await (await fetch(`/cases/${name}.xyz`)).text());
    bench(xtb, mol, 30); // warm up the JIT tiers
    const r = bench(xtb, mol, 200);
    results[name] = { cold_ms: r.cold_ms, warm_mean_ms: r.warm_mean_ms, energy: r.out?.[1], error: r.error };
  }
  postMessage({ load_ms, results });
} catch (e) {
  postMessage({ error: String(e) });
}
