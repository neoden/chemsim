// Writes the stage-0 test geometries (angstrom) into ../cases/*.xyz
import { writeFileSync } from 'node:fs';
const dir = new URL('../cases/', import.meta.url);
const save = (name, atoms, note) =>
  writeFileSync(new URL(name + '.xyz', dir),
    `${atoms.length}\n${note}\n` + atoms.map(([e, x, y, z]) => `${e} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}`).join('\n') + '\n');

save('h2o', [['O', 0, 0, 0.1173], ['H', 0, 0.7572, -0.4692], ['H', 0, -0.7572, -0.4692]], 'water');

// planar BH3 at the origin, NH3 on the z axis with its lone pair towards boron, B-N 3.0 A
const pair = [['B', 0, 0, 0]];
for (const a of [0, 120, 240]) pair.push(['H', 1.19 * Math.cos(a * Math.PI / 180), 1.19 * Math.sin(a * Math.PI / 180), 0]);
pair.push(['N', 0, 0, 3.0]);
for (const a of [60, 180, 300]) pair.push(['H', 0.938 * Math.cos(a * Math.PI / 180), 0.938 * Math.sin(a * Math.PI / 180), 3.38]);
save('bh3_nh3', pair, 'BH3 + NH3, B-N 3.0 A');

// all-trans n-hexane, 20 atoms, idealized geometry
const hex = [], dx = 1.268, dy = 0.856, n = 6;
for (let i = 0; i < n; i++) {
  const up = i % 2, x = i * dx, y = up * dy, s = up ? 1 : -1;
  hex.push(['C', x, y, 0]);
  hex.push(['H', x, y + s * 0.63, 0.89], ['H', x, y + s * 0.63, -0.89]);
  if (i === 0) hex.push(['H', x - 1.09 * 0.829, y - s * 1.09 * 0.559, 0]);
  if (i === n - 1) hex.push(['H', x + 1.09 * 0.829, y - s * 1.09 * 0.559, 0]);
}
save('hexane', hex, 'n-hexane, idealized');
