// Internal units are atomic (bohr, hartree). The UI shows angstrom, eV, kelvin, fs.
export const BOHR_ANGSTROM = 0.529177210903;
export const HARTREE_EV = 27.211386245988;
// Elements the app offers: the first 36. GFN2-xTB gives them s, p and (from Ne on, for most)
// d functions; the engine itself goes up to Rn.
const ELEMENTS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split(' ');
const table = (values) => Object.fromEntries(values.map((v, i) => [i + 1, v]));
export const SYMBOLS = table(ELEMENTS);
export const NUMBERS = Object.fromEntries(Object.entries(SYMBOLS).map(([z, el]) => [el, +z]));
export const FS = 41.341373335; // atomic units of time in a femtosecond
export const KB = 3.166811563e-6; // Boltzmann constant, hartree per kelvin
export const AMU = 1822.888486209; // electron masses in an atomic mass unit
export const MASS_AMU = table([1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948,
  39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798]);
// valence electrons of each element in GFN2-xTB (the d metals count their d electrons; zinc does not)
export const VALENCE = table([1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8]);
// covalent radii, angstrom (Cordero et al., 2008)
export const COVALENT = table([0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06,
  2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16]);
