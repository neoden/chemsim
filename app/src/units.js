// Internal units are atomic (bohr, hartree). The UI shows angstrom, eV, kelvin, fs.
export const BOHR_ANGSTROM = 0.529177210903;
export const HARTREE_EV = 27.211386245988;
// Elements the app offers: those for which GFN2-xTB uses only s and p functions (the renderer
// draws nothing else yet).
export const SYMBOLS = { 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 19: 'K' };
export const NUMBERS = Object.fromEntries(Object.entries(SYMBOLS).map(([z, el]) => [el, +z]));
export const FS = 41.341373335; // atomic units of time in a femtosecond
export const KB = 3.166811563e-6; // Boltzmann constant, hartree per kelvin
export const AMU = 1822.888486209; // electron masses in an atomic mass unit
export const MASS_AMU = { 1: 1.008, 2: 4.0026, 3: 6.94, 4: 9.0122, 5: 10.81, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998, 11: 22.990, 19: 39.098 };
// valence electrons of each element in GFN2-xTB
export const VALENCE = { 1: 1, 2: 2, 3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 11: 1, 19: 1 };
