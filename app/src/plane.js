// Geometry of the cutting plane. A plane is { origin, normal, e1 }: a point on it, the unit
// normal (in "knife" mode the side the normal points to is cut away) and a unit in-plane
// vector that is "right" on the flat 2D plot. All in bohr, molecule coordinates.
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const scale = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = (a) => Math.hypot(a[0], a[1], a[2]);
export const unit = (a) => scale(a, 1 / norm(a));

// component of v perpendicular to unit vector n, normalized; null if v is parallel to n
function perpendicular(v, n) {
  const w = sub(v, scale(n, dot(v, n)));
  return norm(w) > 1e-6 ? unit(w) : null;
}

// e1 in the plane, as close as possible to `hint` (then to the fallbacks)
function inPlane(normal, ...hints) {
  for (const h of [...hints, [1, 0, 0], [0, 1, 0], [0, 0, 1]]) {
    const e = h && perpendicular(h, normal);
    if (e) return e;
  }
}

// make the normal point towards the camera, so that the cut face is what one sees
const facing = (normal, toCamera) => (dot(normal, toCamera) < 0 ? scale(normal, -1) : normal);

// cam: { ro, fw, rt, up } from cameraBasis()
export function planeFacingCamera(cam) {
  const normal = unit(cam.ro);
  return { origin: [0, 0, 0], normal, e1: inPlane(normal, cam.rt) };
}

export function planeThroughThree(a, b, c, cam) {
  const n = cross(sub(b, a), sub(c, a));
  if (norm(n) < 1e-6) return null; // the three atoms are on one line
  const normal = facing(unit(n), cam.ro);
  return { origin: scale(add(add(a, b), c), 1 / 3), normal, e1: inPlane(normal, cam.rt, sub(b, a)) };
}

// contains the bond a-b; of all such planes, the one that faces the camera most
export function planeAlongBond(a, b, cam) {
  const axis = unit(sub(b, a));
  const normal = perpendicular(cam.ro, axis) ?? perpendicular(cam.up, axis);
  return { origin: scale(add(a, b), 0.5), normal, e1: dot(axis, cam.rt) < 0 ? scale(axis, -1) : axis };
}

// perpendicular to the bond a-b through its midpoint
export function planeAcrossBond(a, b, cam) {
  const normal = facing(unit(sub(b, a)), cam.ro);
  return { origin: scale(add(a, b), 0.5), normal, e1: inPlane(normal, cam.rt, cam.up) };
}

export const planeOffset = (p) => dot(p.normal, p.origin); // signed distance from the molecule centre
export const planeCentre = (p) => scale(p.normal, planeOffset(p)); // point of the plane nearest to the centre
export const planeUp = (p) => cross(p.normal, p.e1);

export function shifted(p, distance) {
  return { ...p, origin: add(p.origin, scale(p.normal, distance)) };
}

export function flipped(p) {
  return { ...p, normal: scale(p.normal, -1) };
}

function rotate(v, axis, angle) {
  const c = Math.cos(angle), s = Math.sin(angle);
  return add(add(scale(v, c), scale(cross(axis, v), s)), scale(axis, dot(axis, v) * (1 - c)));
}

// tilt about the plane's own origin (the point it was snapped to, carried along by shifts)
export function tilted(p, axis, angle) {
  const normal = unit(rotate(p.normal, axis, angle));
  return { origin: p.origin, normal, e1: inPlane(normal, rotate(p.e1, axis, angle)) };
}
