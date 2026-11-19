// The flat 2D plot of the cut: the same field as on the plane in the 3D view, with the
// nuclei, a scale bar and a caption.
import { BOHR_ANGSTROM } from './units.js';
import { dot, sub, planeCentre, planeUp, planeOffset } from './plane.js';

const NEAR = 0.3; // bohr: nuclei closer than this to the plane are drawn as lying in it

export function drawSlice(canvas, caption, renderer, plane, field, molecule) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const size = Math.max(64, Math.min(640, Math.round(canvas.clientWidth * dpr)));
  if (canvas.width !== size) { canvas.width = size; canvas.height = size; }
  const ctx = canvas.getContext('2d');
  const centre = planeCentre(plane), e2 = planeUp(plane), half = renderer.scene.planeRadius;
  const pixels = renderer.renderSlice(size, centre, plane.e1, e2, half, field);
  ctx.putImageData(new ImageData(pixels, size, size), 0, 0);

  const toPx = (u, v) => [(u / half * 0.5 + 0.5) * size, (0.5 - v / half * 0.5) * size];
  ctx.font = `${Math.round(11 * dpr)}px system-ui, sans-serif`;
  ctx.textBaseline = 'middle';
  const p = molecule.positions;
  let inPlane = 0;
  for (let k = 0; k < molecule.elements.length; k++) {
    const r = sub([p[3 * k], p[3 * k + 1], p[3 * k + 2]], centre);
    const [x, y] = toPx(dot(r, plane.e1), dot(r, e2)), h = dot(r, plane.normal);
    const near = Math.abs(h) < NEAR;
    if (near) inPlane++;
    ctx.beginPath();
    ctx.arc(x, y, (near ? 2.6 : 2.2) * dpr, 0, 2 * Math.PI);
    if (near) { ctx.fillStyle = '#fff'; ctx.fill(); }
    else { ctx.strokeStyle = 'rgba(255,255,255,.4)'; ctx.lineWidth = dpr; ctx.stroke(); }
    ctx.fillStyle = near ? 'rgba(255,255,255,.9)' : 'rgba(255,255,255,.4)';
    ctx.fillText(molecule.elements[k], x + 5 * dpr, y - 7 * dpr);
  }
  // scale bar: 1 angstrom
  const bar = (1 / BOHR_ANGSTROM) / (2 * half) * size, m = 10 * dpr;
  ctx.strokeStyle = 'rgba(255,255,255,.8)';
  ctx.lineWidth = 1.5 * dpr;
  ctx.beginPath(); ctx.moveTo(m, size - m); ctx.lineTo(m + bar, size - m); ctx.stroke();
  ctx.fillStyle = 'rgba(255,255,255,.8)';
  ctx.fillText('1 Å', m + bar + 5 * dpr, size - m);

  const off = planeOffset(plane) * BOHR_ANGSTROM + 0; // + 0: no negative zero
  caption.textContent =
    (field === 'rho' ? 'Square root of the total density.' : 'Wave function of the selected orbital; the white line is a node.') +
    ` Contour lines at every ×2. Nuclei in the plane: ${inPlane}, the rest shown as circles. ` +
    `Offset from the centre: ${(Math.abs(off) < 0.005 ? 0 : off).toFixed(2).replace('-', '−')} Å.`;
}
