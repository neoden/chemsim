// Energy against time during the dynamics: total and potential energy, in eV relative to
// the potential energy at the start.
import { HARTREE_EV } from './units.js';

const TOTAL = '#6aa9ff', POTENTIAL = '#b9bfcc', GRID = 'rgba(255,255,255,.10)', TEXT = '#9299a8';

// history: { t: [fs], etot: [Eh], epot: [Eh] }
export function drawEnergy(canvas, history) {
  const dpr = Math.min(window.devicePixelRatio || 1, 2), w = Math.round(canvas.clientWidth * dpr), h = Math.round(canvas.clientHeight * dpr);
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
  const ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, w, h);
  const n = history.t.length;
  if (n < 2) return;
  const e0 = history.epot[0], ev = (x) => (x - e0) * HARTREE_EV;
  let lo = Infinity, hi = -Infinity;
  for (let i = 0; i < n; i++) { lo = Math.min(lo, ev(history.epot[i]), ev(history.etot[i])); hi = Math.max(hi, ev(history.epot[i]), ev(history.etot[i])); }
  if (hi - lo < 0.02) { const mid = (hi + lo) / 2; lo = mid - 0.01; hi = mid + 0.01; }
  const pad = { l: 34 * dpr, r: 6 * dpr, t: 18 * dpr, b: 14 * dpr };
  const t0 = history.t[0], t1 = history.t[n - 1];
  const X = (t) => pad.l + (t - t0) / (t1 - t0) * (w - pad.l - pad.r), Y = (v) => pad.t + (hi - v) / (hi - lo) * (h - pad.t - pad.b);
  ctx.font = `${Math.round(10 * dpr)}px system-ui, sans-serif`;
  ctx.fillStyle = TEXT;
  ctx.strokeStyle = GRID;
  ctx.lineWidth = dpr;
  ctx.textBaseline = 'middle';
  ctx.textAlign = 'right';
  for (const v of [lo, hi]) {
    ctx.beginPath(); ctx.moveTo(pad.l, Y(v)); ctx.lineTo(w - pad.r, Y(v)); ctx.stroke();
    ctx.fillText(v.toFixed(2).replace('-', '−'), pad.l - 4 * dpr, Y(v));
  }
  if (lo < 0 && hi > 0) { ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(w - pad.r, Y(0)); ctx.stroke(); }
  ctx.textBaseline = 'alphabetic';
  ctx.fillText(t1 >= 2000 ? `${(t1 / 1000).toFixed(2)} ps` : `${t1.toFixed(0)} fs`, w - pad.r, h - 2 * dpr);
  ctx.textAlign = 'left';
  const line = (series, colour, width) => {
    ctx.strokeStyle = colour; ctx.lineWidth = width * dpr; ctx.lineJoin = 'round';
    ctx.beginPath();
    for (let i = 0; i < n; i++) { const x = X(history.t[i]), y = Y(ev(series[i])); if (i) ctx.lineTo(x, y); else ctx.moveTo(x, y); }
    ctx.stroke();
  };
  line(history.epot, POTENTIAL, 1);
  line(history.etot, TOTAL, 1.6);
  // direct labels instead of a legend box
  ctx.fillStyle = TEXT; ctx.fillText('eV', 2 * dpr, pad.t - 4 * dpr);
  ctx.fillStyle = TOTAL; ctx.fillText('total', pad.l + 4 * dpr, pad.t - 4 * dpr);
  ctx.fillStyle = POTENTIAL; ctx.fillText('potential', pad.l + 52 * dpr, pad.t - 4 * dpr);
}
