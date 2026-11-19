#!/bin/sh
# Builds the GFN2-xTB WebAssembly module and puts it where the app loads it from.
set -e
cd "$(dirname "$0")/../engine/xtb-wasm"
cargo build --release --lib --target wasm32-wasip1
cargo build --release --bin native
cp target/wasm32-wasip1/release/xtb_wasm.wasm ../../app/engine/xtb_wasm.wasm
ls -l ../../app/engine/xtb_wasm.wasm
