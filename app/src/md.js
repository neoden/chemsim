// Born–Oppenheimer molecular dynamics of the nuclei: classical nuclei on the forces of the
// quantum engine. No DOM, no worker: usable in the browser worker and in node tests alike.
//
// Integrator: velocity Verlet. With the thermostat on it becomes the BAOAB Langevin scheme
// (half kick, half drift, friction + noise, half drift, half kick), which is exactly
// velocity Verlet when the thermostat is off.
//
// Why Langevin: the systems here are single small molecules (3–12 atoms) in vacuum. A
// thermostat that rescales all velocities at once (Berendsen, Nosé–Hoover, stochastic
// rescaling) only fixes the total kinetic energy and relies on the molecule itself to share
// it between vibrations; with so few, nearly harmonic modes that sharing is slow or never
// happens, and energy can pile up in rotation. Langevin couples every atom to the bath
// separately, samples the right distribution for any size, and stays stable up to 10 000 K.
// Its price is a mild damping of coherent vibrations, set by the coupling time tau.
import { FS, KB, AMU, MASS_AMU } from './units.js';

function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export class Dynamics {
  // force(positions, wantOrbitals) -> { ok, energy, forces, ... } or { ok: false, converged, error }.
  // Everything in atomic units except the parameters, which are in the UI's units.
  constructor({ numbers, positions, force, seed = 1, wallRadius = Infinity }) {
    this.numbers = numbers;
    this.n = numbers.length;
    this.mass = Float64Array.from(numbers, (z) => MASS_AMU[z] * AMU);
    this.x = Float64Array.from(positions);
    this.v = new Float64Array(3 * this.n);
    this.f = new Float64Array(3 * this.n);
    this.force = force;
    this.uniform = mulberry32(seed);
    this.params = {
      dt: 0.5, // fs
      thermostat: true,
      temperature: 300, // K
      tau: 100, // fs, coupling time of the thermostat
      wallRadius, // bohr: beyond it a soft wall pushes atoms back towards the centre
      wallK: 0.1, // hartree / bohr^2
    };
    // An atom held by the mouse: { atom, target: [x, y, z] (bohr) }, or { atom, toAtom, rest }
    // for a target that follows another atom: the point `rest` bohr from it on the line to the
    // held one. Two things act while it is held.
    // A spring on the atom itself; its pull is capped so that it turns and bends the molecule
    // instead of tearing the atom off, and the atom's motion is damped so that it does not ring.
    // With grab.piece (the atoms of the molecule the held one belongs to), a carrier as well:
    // the same acceleration for every atom of the piece, steering the piece towards the speed
    // min(distance / carryTime, carrySpeed) in the direction of the target. It strains no bond,
    // so it may be strong, and a heavy molecule follows as readily as a light one.
    this.grab = null;
    this.grabK = 0.08; // hartree / bohr^2
    this.grabMaxForce = 0.03; // hartree / bohr (bonds here resist with 0.1 or more)
    this.grabDamping = 25; // fs
    this.carrySpeed = 0.04; // bohr / fs: about three thermal speeds of a small molecule at 300 K
    this.carryTime = 20; // fs: the approach slows down over the last carrySpeed * carryTime bohr
    this.carryRelax = 10; // fs: how fast the piece takes the wanted speed
    this.carryMaxForce = 0.12; // hartree / bohr on the whole piece: what it presses with on an obstacle
    this.time = 0; // fs
    this.steps = 0;
    this.epot = 0; // engine energy + wall energy
    this.result = null; // last engine result
    this.ready = false;
  }

  gaussian() {
    let u = 0;
    while (u === 0) u = this.uniform();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * this.uniform());
  }

  // engine forces plus the wall, at this.x; false if the engine failed
  evaluate(wantOrbitals) {
    const r = this.force(this.x, wantOrbitals);
    if (!r.ok) { this.error = r; return false; }
    this.result = r;
    this.f.set(r.forces);
    let wall = 0;
    const { wallRadius: R, wallK: k } = this.params;
    for (let i = 0; i < this.n; i++) {
      const d = Math.hypot(this.x[3 * i], this.x[3 * i + 1], this.x[3 * i + 2]);
      if (d > R) {
        wall += 0.5 * k * (d - R) ** 2;
        for (let c = 0; c < 3; c++) this.f[3 * i + c] -= k * (d - R) * this.x[3 * i + c] / d;
      }
    }
    // the hand is outside work: its energy is not counted in epot
    if (this.grab) {
      const g = this.grab, a = g.atom, d = this.grabVector(), len = Math.hypot(...d);
      const pull = Math.min(this.grabK * len, this.grabMaxForce);
      if (len > 1e-9) for (let c = 0; c < 3; c++) this.f[3 * a + c] += pull * d[c] / len;
      if (g.piece && len > 1e-9) {
        const { mass, velocity } = this.pieceMotion(g.piece), want = Math.min(len / this.carryTime, this.carrySpeed) / FS;
        const F = [0, 1, 2].map((c) => mass * (want * d[c] / len - velocity[c]) / (this.carryRelax * FS));
        const size = Math.hypot(...F), k = Math.min(1, this.carryMaxForce / Math.max(size, 1e-12)) / mass;
        for (const i of g.piece) for (let c = 0; c < 3; c++) this.f[3 * i + c] += k * this.mass[i] * F[c];
      }
    }
    this.epot = r.energy + wall;
    return true;
  }

  // from the held atom to where the hand wants it
  grabVector() {
    const g = this.grab, a = g.atom, x = this.x;
    if (g.toAtom === undefined) return [0, 1, 2].map((c) => g.target[c] - x[3 * a + c]);
    const u = [0, 1, 2].map((c) => x[3 * a + c] - x[3 * g.toAtom + c]), l = Math.hypot(...u) || 1;
    return u.map((v) => v * (g.rest / l - 1));
  }

  // mass and centre-of-mass velocity of a set of atoms
  pieceMotion(piece) {
    const velocity = [0, 0, 0];
    let mass = 0;
    for (const i of piece) { mass += this.mass[i]; for (let c = 0; c < 3; c++) velocity[c] += this.mass[i] * this.v[3 * i + c]; }
    return { mass, velocity: velocity.map((p) => p / mass) };
  }

  // Put the nuclei somewhere else by hand (velocities are kept); false if the engine failed.
  moveTo(positions, wantOrbitals = true) {
    const x0 = this.x.slice();
    this.x.set(positions);
    if (this.evaluate(wantOrbitals)) return true;
    this.x = x0;
    return false;
  }

  // the middle of the nuclei
  centre() {
    const c = [0, 0, 0];
    for (let k = 0; k < this.x.length; k++) c[k % 3] += this.x[k] / this.n;
    return c;
  }

  // Moves everything by one vector, so that the middle of the nuclei is at `point`. For a run
  // without a wall, where nothing is tied to a place: the forces and the electrons do not change.
  centreAt(point) {
    const by = this.centre().map((v, k) => point[k] - v);
    for (let k = 0; k < this.x.length; k++) this.x[k] += by[k % 3];
    if (this.grab?.target) this.grab.target = this.grab.target.map((v, k) => v + by[k]);
  }

  // The same run with only the atoms `keep` (indices): where they are, how they move, the
  // clock and the settings are carried over. `force` computes for the smaller set. The new
  // run still needs its init().
  kept(keep, force) {
    const pick = (a) => Float64Array.from({ length: 3 * keep.length }, (_, k) => a[3 * keep[(k / 3) | 0] + k % 3]);
    const d = new Dynamics({ numbers: keep.map((i) => this.numbers[i]), positions: pick(this.x), force });
    Object.assign(d.params, this.params);
    d.uniform = this.uniform;
    d.v = pick(this.v);
    d.time = this.time;
    d.steps = this.steps;
    return d;
  }

  // forces at the starting geometry; must succeed before step()
  init(wantOrbitals = true) {
    this.ready = this.evaluate(wantOrbitals);
    return this.ready;
  }

  removeDrift() {
    const p = [0, 0, 0];
    let m = 0;
    for (let i = 0; i < this.n; i++) { m += this.mass[i]; for (let c = 0; c < 3; c++) p[c] += this.mass[i] * this.v[3 * i + c]; }
    for (let i = 0; i < this.n; i++) for (let c = 0; c < 3; c++) this.v[3 * i + c] -= p[c] / m;
  }

  get dof() { return Math.max(1, 3 * this.n - 3); }
  get ekin() {
    let e = 0;
    for (let i = 0; i < this.n; i++) e += 0.5 * this.mass[i] * (this.v[3 * i] ** 2 + this.v[3 * i + 1] ** 2 + this.v[3 * i + 2] ** 2);
    return e;
  }
  get temperature() { return 2 * this.ekin / (KB * this.dof); }
  get etot() { return this.epot + this.ekin; }

  // Maxwell–Boltzmann velocities without drift, scaled to exactly `kelvin`
  thermalize(kelvin) {
    for (let i = 0; i < this.n; i++) {
      const s = Math.sqrt(KB * Math.max(kelvin, 0) / this.mass[i]);
      for (let c = 0; c < 3; c++) this.v[3 * i + c] = s * this.gaussian();
    }
    this.removeDrift();
    const t = this.temperature;
    if (t > 0) { const k = Math.sqrt(kelvin / t); for (let i = 0; i < this.v.length; i++) this.v[i] *= k; }
  }

  // One time step. Returns false if the electronic problem failed; the state is then left
  // as it was before the step.
  step(wantOrbitals = false) {
    const { dt, thermostat, temperature, tau } = this.params, h = dt * FS, n = this.n;
    const x0 = this.x.slice(), v0 = this.v.slice();
    for (let i = 0; i < n; i++) {
      const a = 0.5 * h / this.mass[i];
      for (let c = 0; c < 3; c++) { const k = 3 * i + c; this.v[k] += a * this.f[k]; this.x[k] += 0.5 * h * this.v[k]; }
    }
    if (thermostat) {
      const c1 = Math.exp(-dt / tau), c2 = Math.sqrt(1 - c1 * c1);
      for (let i = 0; i < n; i++) {
        const s = c2 * Math.sqrt(KB * temperature / this.mass[i]);
        for (let c = 0; c < 3; c++) this.v[3 * i + c] = c1 * this.v[3 * i + c] + s * this.gaussian();
      }
      this.removeDrift(); // the random kicks do not conserve momentum
    }
    for (let k = 0; k < 3 * n; k++) this.x[k] += 0.5 * h * this.v[k];
    const f0 = this.f.slice(), e0 = this.epot;
    if (!this.evaluate(wantOrbitals)) {
      this.x = x0; this.v = v0; this.f = f0; this.epot = e0;
      return false;
    }
    for (let i = 0; i < n; i++) {
      const a = 0.5 * h / this.mass[i];
      for (let c = 0; c < 3; c++) this.v[3 * i + c] += a * this.f[3 * i + c];
    }
    if (this.grab) {
      // with a carrier the atom is damped relative to its molecule, which is meant to move
      const k = Math.exp(-dt / this.grabDamping), a = this.grab.atom, base = this.grab.piece ? this.pieceMotion(this.grab.piece).velocity : [0, 0, 0];
      for (let c = 0; c < 3; c++) this.v[3 * a + c] = base[c] + k * (this.v[3 * a + c] - base[c]);
    }
    this.time += dt;
    this.steps++;
    return true;
  }
}
