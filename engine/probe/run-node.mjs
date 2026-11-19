import { readFileSync } from 'node:fs';
import { loadXtb, parseXyz, bench } from './xtb.mjs';
const wasm = readFileSync(new URL('../xtb-wasm/target/wasm32-wasip1/release/xtb_wasm.wasm', import.meta.url));
const xtb = await loadXtb(wasm);
const name = process.argv[2];
const mol = parseXyz(readFileSync(new URL(`../cases/${name}.xyz`, import.meta.url), 'utf8'));
console.log(JSON.stringify(bench(xtb, mol, +(process.argv[3] ?? 20))));
