//! Thin wrapper around axtb (GFN2-xTB) with a flat C ABI, usable both natively and from
//! WebAssembly without any JS binding generator.
//!
//! Units at this boundary are atomic: bohr, hartree, hartree/bohr.
//!
//! Result buffer layout (f64):
//!   [0] status: 0 ok, 1 SCF not converged, 2 other error (message via xtb_err_*)
//!   [1] energy  [2] nat  [3] nao  [4] nmo  [5] SCF iterations  [6] nshell  [7] has_gradient
//!   then gradient (3*nat, dE/dR), orbital energies (nmo), occupations (nmo, both spins
//!   together, 0..2), occupations of the alpha spin alone (nmo, 0..1),
//!   coefficients (nmo*nao, orbital by orbital), and per shell:
//!   atom, l, nprim, alpha[nprim], coeff[nprim].
//! AO order: atom by atom, shell by shell, m = -l..l (p: y, z, x).

use axtb::{Calculator, Error, Method, Options};
use std::sync::Mutex;

pub const HEADER: usize = 8;

static LAST_ERROR: Mutex<String> = Mutex::new(String::new());

fn set_error(msg: String) {
    if let Ok(mut e) = LAST_ERROR.lock() {
        *e = msg;
    }
}

pub struct Engine {
    calc: Calculator,
    nat: usize,
    pub out: Vec<f64>,
}

impl Engine {
    /// `unpaired`: number of unpaired electrons (alpha minus beta); `None` takes the lowest
    /// that the electron count allows (0 or 1).
    pub fn new(numbers: &[usize], positions: &[[f64; 3]], charge: f64, etemp: f64, unpaired: Option<usize>) -> Result<Self, Error> {
        let mut options = Options::default();
        options.scf.etemp = etemp;
        options.scf.threads = 1;
        let calc = Calculator::new(Method::Gfn2, numbers, positions, charge, unpaired, options)?;
        Ok(Self { calc, nat: numbers.len(), out: Vec::new() })
    }

    pub fn set_etemp(&mut self, etemp: f64) {
        self.calc.options.scf.etemp = etemp;
    }

    /// SCF accuracy multiplier (axtb: energy threshold 1e-6 * acc, residual 2e-5 * acc).
    pub fn set_accuracy(&mut self, accuracy: f64) {
        self.calc.options.scf.accuracy = accuracy;
    }

    pub fn reset_guess(&mut self) {
        self.calc.reset_guess();
    }

    /// Runs one calculation and fills `self.out`; returns the status code.
    pub fn singlepoint(&mut self, positions: &[[f64; 3]], gradient: bool, orbitals: bool) -> i32 {
        self.calc.options.scf.keep_orbitals = orbitals;
        self.out.clear();
        self.out.resize(HEADER, 0.0);
        let sp = match self.calc.singlepoint(positions, gradient) {
            Ok(sp) => sp,
            Err(e) => {
                let status = if matches!(e, Error::ScfNotConverged { .. }) { 1 } else { 2 };
                set_error(e.to_string());
                self.out[0] = status as f64;
                return status;
            }
        };
        let nat = self.nat;
        let out = &mut self.out;
        out[1] = sp.energy;
        out[2] = nat as f64;
        out[5] = sp.iterations as f64;
        match &sp.gradient {
            Some(g) => {
                out[7] = 1.0;
                out.extend(g.iter().flatten());
            }
            None => out.extend(std::iter::repeat(0.0).take(3 * nat)),
        }
        if let Some(orb) = &sp.orbitals {
            let set = &orb.sets[0];
            let nmo = set.energies.len();
            let nao = set.coefficients.len() / nmo.max(1);
            out[3] = nao as f64;
            out[4] = nmo as f64;
            out.extend(&set.energies);
            // an open shell comes as an alpha and a beta set of the same orbitals
            match orb.sets.get(1) {
                Some(beta) => {
                    out.extend(set.occupations.iter().zip(&beta.occupations).map(|(a, b)| a + b));
                    out.extend(&set.occupations);
                }
                None => {
                    out.extend(&set.occupations);
                    out.extend(set.occupations.iter().map(|o| 0.5 * o));
                }
            }
            out.extend(&set.coefficients);
            let mut nshell = 0;
            for (atom, shells) in orb.shells.iter().enumerate() {
                for sh in shells {
                    nshell += 1;
                    out.extend([atom as f64, sh.l as f64, sh.nprim as f64]);
                    out.extend(&sh.alpha[..sh.nprim]);
                    out.extend(&sh.coeff[..sh.nprim]);
                }
            }
            out[6] = nshell as f64;
        } else {
            let nmo = sp.orbital_energies.len();
            out[4] = nmo as f64;
            out.extend(&sp.orbital_energies);
            out.extend(&sp.occupations);
            out.extend(sp.occupations.iter().map(|o| 0.5 * o)); // not meaningful without the orbitals
        }
        0
    }
}

fn positions_from(ptr: *const f64, nat: usize) -> Vec<[f64; 3]> {
    let flat = unsafe { std::slice::from_raw_parts(ptr, 3 * nat) };
    flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn install_panic_hook() {
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        std::panic::set_hook(Box::new(|info| set_error(format!("panic: {info}"))));
    });
}

#[no_mangle]
pub extern "C" fn xtb_alloc(len: usize) -> *mut f64 {
    let mut v = vec![0.0f64; len];
    let p = v.as_mut_ptr();
    std::mem::forget(v);
    p
}

#[no_mangle]
pub unsafe extern "C" fn xtb_free(ptr: *mut f64, len: usize) {
    drop(Vec::from_raw_parts(ptr, len, len));
}

/// `numbers`: nat atomic numbers as f64; `positions`: 3*nat in bohr; `unpaired`: number of
/// unpaired electrons, negative for the lowest possible. Null on error.
#[no_mangle]
pub unsafe extern "C" fn xtb_create(numbers: *const f64, nat: usize, positions: *const f64, charge: f64, etemp: f64, unpaired: i32) -> *mut Engine {
    install_panic_hook();
    let numbers: Vec<usize> = std::slice::from_raw_parts(numbers, nat).iter().map(|&z| z as usize).collect();
    match Engine::new(&numbers, &positions_from(positions, nat), charge, etemp, usize::try_from(unpaired).ok()) {
        Ok(e) => Box::into_raw(Box::new(e)),
        Err(e) => {
            set_error(e.to_string());
            std::ptr::null_mut()
        }
    }
}

#[no_mangle]
pub unsafe extern "C" fn xtb_destroy(e: *mut Engine) {
    drop(Box::from_raw(e));
}

#[no_mangle]
pub unsafe extern "C" fn xtb_set_etemp(e: *mut Engine, etemp: f64) {
    (*e).set_etemp(etemp);
}

#[no_mangle]
pub unsafe extern "C" fn xtb_set_accuracy(e: *mut Engine, accuracy: f64) {
    (*e).set_accuracy(accuracy);
}

#[no_mangle]
pub unsafe extern "C" fn xtb_reset_guess(e: *mut Engine) {
    (*e).reset_guess();
}

#[no_mangle]
pub unsafe extern "C" fn xtb_singlepoint(e: *mut Engine, positions: *const f64, gradient: i32, orbitals: i32) -> i32 {
    let e = &mut *e;
    let pos = positions_from(positions, e.nat);
    e.singlepoint(&pos, gradient != 0, orbitals != 0)
}

#[no_mangle]
pub unsafe extern "C" fn xtb_out_ptr(e: *const Engine) -> *const f64 {
    (*e).out.as_ptr()
}

#[no_mangle]
pub unsafe extern "C" fn xtb_out_len(e: *const Engine) -> usize {
    (*e).out.len()
}

#[no_mangle]
pub extern "C" fn xtb_err_ptr() -> *const u8 {
    LAST_ERROR.lock().map(|e| e.as_ptr()).unwrap_or(std::ptr::null())
}

#[no_mangle]
pub extern "C" fn xtb_err_len() -> usize {
    LAST_ERROR.lock().map(|e| e.len()).unwrap_or(0)
}
