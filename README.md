# chemsim

A browser-based chemistry simulator at the level of atoms.

Status: stage 0 (engine probe), stage 1 (static viewer on the new engine), stage 2 (cutting plane), stage 3 (dynamics), stage 4 (hands-on interaction) and stage 5 (build for demonstration) are complete. The first version is ready.

Elements: all that GFN2-xTB has parameters for, the first 86, H to Rn. It describes them with s, p and d functions (d for most elements from Ne on), and the renderer draws all three kinds. For the main-group elements the d functions only add flexibility; for the d metals they are the valence orbitals themselves. The method has no f functions: the 4f electrons of the lanthanides are counted with the core, so each lanthanide is a metal with three valence electrons. The tables of masses, radii and valence electrons are in `app/src/units.js`; the valence counts are the engine's own (an atom with the wrong parity of unpaired electrons is refused by it).

Unpaired electrons: the engine is given their number (two for O₂, one for radicals). The scene sets it, and it can be changed under "Settings" in the dynamics panel. Spin only decides how electrons are distributed over the orbitals; GFN2-xTB knows no energy difference between "up" and "down".

Your own molecules: "Edit" opens what is in the scene in the molecule editor, "New" starts from methane. A molecule is grown by its hydrogens: every atom arrives with its hydrogens, and a click on a hydrogen puts the chosen element in its place. Dragging from atom to atom draws a bond, a click on a bond makes it double or triple, the eraser removes atoms, bonds or single hydrogens (which leaves a radical). "Relax" lets the engine take the nuclei to the nearest energy minimum; "Done" relaxes and keeps the molecule under "My molecules" in the browser (`localStorage`), from where it works as any molecule of the library. "XYZ…" shows the coordinates as text and takes a pasted geometry. The editor keeps H, C, N, O, F, P, S and Cl on its toolbar; "…" next to them opens the periodic table with all 86 elements, and the last three taken from it stay on the toolbar. Elements with a usual number of bonds (B, C, N, O, F, Si, P, S, Cl, Br, I) arrive with their hydrogens; the others, the metals and the noble gases among them, arrive bare, and their bonds are drawn by dragging from atom to atom.

Molecules from elsewhere: "Find…" takes a compound from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) by its name or number (CID), or builds one from a formula in SMILES; either way it opens in the editor, where "Done" relaxes and keeps it. From PubChem comes the computed 3D structure with its bonds; for compounds that have none there (salts, metal compounds) the geometry is built from the SMILES PubChem gives. "SMILES…" in the editor shows the molecule being edited as SMILES and replaces it the same way. Typing Enter builds what reads as SMILES and looks up anything else. A search sends the typed text to PubChem and is the only request the page ever makes outside its own folder; nothing is sent until the button is pressed. SMILES is read in full except for what the program has no use for: the handedness marks (`@`, `/`, `\`) and isotopes are dropped, so of two mirror forms one comes out, whichever. It is written with whole bond orders (benzene as `C1=CC=CC=C1`); the charge of an ion is put on the atom whose bonds point to it.

## Running

```sh
node tools/serve.mjs 8080      # then open http://127.0.0.1:8080/app/
```

The `app/` folder is the site: static files, ES modules, no bundler and no dependencies. It can be put on any static hosting as is. The calculation runs in a single thread; `SharedArrayBuffer` and the COOP/COEP headers are not needed.

## Deployment

```sh
node tools/build-site.mjs      # builds dist/
```

The contents of `dist/` are the finished site: plain static files with relative links. It can be placed on any static hosting, at the root or in a subfolder. No server side, bundler or special headers are needed. The script also adds what the LGPL requires: `licenses/` (the GPL and LGPL texts, the `NOTICE` from `axtb`, and the notices of every Rust package that went into the module) and `source/` (the wrapper source and instructions for rebuilding the module with your own version of `axtb`). Both pages are linked from the "About the model" window.

## Layout

- `app/` — the application.
  - `src/engine/` — the engine `xtb.js` (GFN2-xTB in WebAssembly) behind the interface described in `worker.js`. It runs in a Web Worker; stale requests are dropped.
  - `src/render/fog.js` — the renderer: ray marching in WebGL2, with the basis and coefficients in data textures; the plane cut lives here too. While the picture changes, the selected orbital is read from a 3D texture filled once per geometry and orbital; a picture that stands still is drawn once more with the orbital evaluated exactly at every ray step.
  - `src/md.js` — molecular dynamics: integrator, thermostat, wall. No DOM, so it runs both in a Worker and in Node.
  - `src/orbitals.js` — the basis functions of a shell (s, p, d) and their order, signs, the overlap matrix, matching orbitals between frames.
  - `src/experiments.js` — the list of ready-made experiments: molecules, settings, a hint and the expected products.
  - `src/scene.js` — a scene of one or two molecules, moving and rotating a fragment.
  - `src/editor.js` — the molecule editor's model: atoms, bonds and the operations on them (hydrogens counted by difference, new atoms placed by repulsion on a sphere). No DOM.
  - `src/editmode.js` — the editor on the page: toolbar, skeleton over the cloud, mouse, undo, the user's molecules in `localStorage`, the "Find…" window.
  - `src/smiles.js` — SMILES: reading it into a bond graph (aromatic rings get their double bonds), writing a molecule back. The geometry for a graph is built by `fromGraph` in `editor.js`. No DOM.
  - `src/pubchem.js` — looking a compound up on PubChem (PUG REST) and taking its 3D record as a molecule.
  - `src/relax.js` — geometry optimization (BFGS), shared by the Worker and the library build.
  - `src/chart.js` — the plot of energy over time.
  - `src/plane.js`, `src/slice2d.js` — the geometry of the cutting plane and the flat 2D plot of the cut.
  - `src/library.json` — the molecule library, generated by `tools/build-library.mjs`.
  - `engine/xtb_wasm.wasm` — the built module, generated by `tools/build-engine.sh`.
- `engine/xtb-wasm/` — a Rust wrapper around the `axtb` crate with a flat C interface. It builds into a WebAssembly module and into a native program for cross-checking.
- `engine/cases/`, `engine/probe/` — geometries and checks from stage 0.
- `third_party/axtb/` — https://github.com/atomsforhumanity/axtb at tag `v0.1.2`, unmodified, as a git submodule (clone with `--recurse-submodules`).
- `tools/` — the engine build, the library build, the local server.
- `tests/` — the checks; everything they save goes to `tests/out/`.

## Building the engine and the library

Requires Rust 1.89+ with the `wasm32-wasip1` target (`rustup target add wasm32-wasip1`) and Node.

```sh
tools/build-engine.sh          # the WebAssembly module in app/engine/ and the native reference build
node tools/build-library.mjs   # optimizes geometries with the xTB engine, writes app/src/library.json
                               # new molecule: add a line to the `molecules` list in this file (coordinates or a z-matrix), then run it
```

The module is built for `wasm32-wasip1` rather than `wasm32-unknown-unknown`: `axtb` times the phases of a calculation with `std::time::Instant`, and under `unknown-unknown` that call aborts the module. Under WASI the clock works; the six WASI calls the module needs are implemented in `app/src/engine/xtb-binding.js`. The `axtb` sources are left unmodified.

## Checks

```sh
node tests/engines.test.mjs                 # the engine: the basis with its d functions, the character of xTB orbitals
node tests/md.test.mjs                      # dynamics: energy conservation, water for 10 ps at 300 K, BH3NH3 dissociation, the wall, the hand that leads a molecule
node tests/editor.test.mjs                  # the editor's model: building by hydrogens, bonds, erasing, closing a ring; what it builds relaxes to the library's energies
node tests/smiles.test.mjs                  # SMILES in and out, the geometries built from it (rings, cages) and their relaxation, PubChem records and answers (played back, nothing is asked)
node tests/adduct.test.mjs                  # BH3 + NH3, brought together "by mouse", form the adduct
node tests/experiments-screen.mjs           # experiment screening: some fifty pairs, brought together "by mouse" or heated, and what they yield
node tests/md-dissociation.mjs              # at what temperature and after how long B–N breaks (8 runs per temperature)
node engine/probe/compare.mjs               # WebAssembly against the native build: energy, forces, orbitals
node tests/run-browser.mjs render-check.html   # rendering in headless Chromium: frames drawn from the orbital grid and the body grid against the exact evaluation, the exact frame of a still picture, screenshots, frame time
node tests/run-browser.mjs plane-check.html    # cutting plane: accuracy of the field on the cut, snapping, mouse control, screenshots
node tests/run-browser.mjs md-check.html       # dynamics in the browser: start, pause, step, reset, steps per second, orbitals without flicker
node tests/run-browser.mjs hands-check.html    # a second molecule, moving and rotating, grabbing a nucleus with a spring, the third color, the energy plot
node tests/run-browser.mjs editor-check.html   # the molecule editor by mouse: build, undo, relax, keep under "My molecules", reload, delete; a molecule from SMILES and one from PubChem (answers played back)
node tests/run-browser.mjs large-check.html    # large molecules: built from SMILES, drawn, relaxed, run in dynamics, and what each step costs on this machine
node tests/run-browser.mjs site-check.html     # the built dist/ site: works from a subfolder, requests nothing external, licenses and sources in place
node tests/run-browser.mjs library-check.html  # every library molecule opens and runs in dynamics
node tests/run-browser.mjs experiments-check.html  # every experiment in the list loads; those that run hands-free give the promised products
node tests/run-browser.mjs speed-check.html    # frame time per molecule: from the orbital grid, exact, and filling the grid
node tests/screenshot.mjs "mol=bh3nh3&wobble=0" tests/out/app.png   # a screenshot of the whole page
```

The rendering checks run on the machine's GPU through ANGLE/OpenGL, as in an ordinary Chromium window on Linux. Flags for `tests/run-browser.mjs`: `--angle=vulkan` (a different backend), `--scale=1.875` (display scale), `--log` (the browser and GPU process log in `tests/out/browser.log`), `--software` (software renderer, very slow). It is worth running on both backends: a driver bug that crashed the GPU process showed up only on OpenGL.

```sh
node tests/run-browser.mjs context-check.html --scale=1.875   # the WebGL context is not lost on a large canvas
node tests/run-browser.mjs gl-bisect.html --log               # after which GL call the context is lost
```

URL parameters for the checks: `mol`, `orb` (`homo`, `lumo`, a number), `body`, `wobble`, `yaw`, `pitch`, `canvas=WxH`, `bodygrid=0` (the body is computed exactly at every point of the ray), `orbgrid=0` (the selected orbital as well: every frame is an exact one, as before the orbital grid), `plane` (`view`, `three:0,1,2`, `along:0,1`, `across:0,1`; the numbers are atom indices), `pmode` (`knife`, `sheet`), `pfield` (`psi`, `rho`), `poffset` (offset in Å), `second` (the second molecule), `experiment` (the experiment key), `md=run`, `edit` (`1`: the editor on the scene, `new`: on a new molecule), `temp` (K), `thermostat=0`, `speed` (fs/s, 0 for maximum), `dt` (fs), `seed`.

## Licenses

The `xtb_wasm.wasm` module statically includes the `axtb`, `axtb-d4` and `axtb-solvation` crates under **LGPL-3.0-or-later**. The `axtb-gxtb` crate (GPL-3.0-only) is not linked in and must not be.

What the LGPL requires when the module is distributed (that is, when the site is published), and how each point is covered:

1. State prominently that `axtb` is used and that it is under the LGPL: the "Licenses" paragraph in the "About the model" window.
2. Include the GPL-3.0 and LGPL-3.0 texts and the `NOTICE` listing the original works: `dist/licenses/`.
3. Provide the source code of the library in the version that is in the module: it is unmodified, and the licenses page links to the exact commit. If `axtb` is ever modified, the changes must be published under the LGPL.
4. Let the recipient rebuild the module with a modified version of `axtb`: the wrapper source and the instructions are in `dist/source/`.

The rest of the application code (JavaScript, shaders) talks to the module through its interface and does not fall under the LGPL. If the whole project repository is public, points 3 and 4 are satisfied automatically.

The other packages in the module (`faer`, `rayon` and their dependencies) are under MIT, Apache-2.0 and similar licenses; their texts are collected in `dist/licenses/third-party/`. Twelve MIT packages ship no license file of their own: for those, the standard MIT text is written out with the authors from the package manifest.

The project itself (JavaScript, shaders, the wrapper) is under the MIT license, see `LICENSE`. The MIT license does not restrict changes to `axtb` or the rebuilding of the module, so it meets the LGPL. The build script copies `LICENSE` into `dist/licenses/` and `dist/source/xtb-wasm/`.
