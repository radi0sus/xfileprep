// e-statistics.js
// Wilson-plot scaling and E² (normalized structure factor) statistics.
//
// This answers a question systematic absences can never answer:
// Friedel's law makes every diffraction pattern APPEAR centrosymmetric
// (up to anomalous-scattering effects) regardless of whether the crystal
// STRUCTURE itself has a real inversion center — so absences alone
// cannot tell e.g. Cc from C2/c, or Pn from P21/n whenever the acentric
// option's missing symmetry element doesn't add any absence of its own.
// E-statistics instead looks at the INTENSITY DISTRIBUTION itself: a
// centrosymmetric structure's normalized intensities follow a different
// probability distribution (Gaussian-derived) than an acentric one, and
// that shows up in simple summary statistics.
//
// Reference values for ideal (untwinned, unbiased) data — Wilson (1949);
// Howells, Phillips & Rogers (1950) — are the standard comparison XPREP
// and PLATON both quote:
//   <|E^2 - 1|>    ≈ 0.968 (centrosymmetric)  vs  0.736 (acentric)
//   <|E|>          ≈ 0.798 (centrosymmetric)  vs  0.886 (acentric)
//   fraction |E|>3 ≈ 0.30%  (centrosymmetric) vs  0.01%  (acentric)
// Real data rarely sits exactly on either value — heavy atoms, NCS,
// twinning, and a poor B-factor fit all pull it around — so this is a
// diagnostic to weigh alongside everything else, not a verdict to trust
// blindly on its own.

const WILSON_REFERENCE = {
  centrosymmetric: { meanE2minus1: 0.968, meanE: 0.798, fractionEgt3: 0.003 },
  acentric: { meanE2minus1: 0.736, meanE: 0.886, fractionEgt3: 0.0001 },
};

// Ordinary least squares for a small set of (x, y) points.
function linearRegression(points) {
  const n = points.length;
  const sx = points.reduce((a, p) => a + p.x, 0);
  const sy = points.reduce((a, p) => a + p.y, 0);
  const sxx = points.reduce((a, p) => a + p.x * p.x, 0);
  const sxy = points.reduce((a, p) => a + p.x * p.y, 0);
  const denom = n * sxx - sx * sx;
  if (denom === 0) return { slope: 0, intercept: sy / n };
  const slope = (n * sxy - sx * sy) / denom;
  const intercept = (sy - slope * sx) / n;
  return { slope, intercept };
}

// Sum of f_j(s)^2 over the unit cell contents at one resolution point s.
// Returns the elements missing from the scattering-factor table (if any)
// alongside the sum, so the caller can report exactly what's missing
// rather than silently mis-scaling the whole plot.
function sumSquaredScatteringFactors(composition, s) {
  let sum = 0;
  const missing = [];
  for (const { element, count } of composition) {
    const f = scatteringFactor(element, s);
    if (f === null) { missing.push(element); continue; }
    sum += count * f * f;
  }
  return { sum, missing };
}

// Bins reflections into resolution shells with roughly equal reflection
// counts per shell — more robust than equal-width sin(theta)/lambda bins
// when reflections are unevenly spread across resolution (typical for
// real unmerged data, which thins out sharply at high resolution).
function binByResolution(reflections, sValues, targetShells = 20) {
  const indexed = reflections.map((r, i) => ({ r, s: sValues[i] })).sort((a, b) => a.s - b.s);
  const perShell = Math.max(20, Math.ceil(indexed.length / targetShells));
  const shells = [];
  for (let start = 0; start < indexed.length; start += perShell) {
    shells.push(indexed.slice(start, start + perShell));
  }
  return shells;
}

// Full computation: merges reflections under the given Laue-class point
// group (Wilson/E-statistics should run on symmetry-unique data, not
// raw redundant observations — see mergeReflections in laue-merge.js),
// fits the Wilson plot, then computes E² per unique reflection and the
// summary statistics used to judge centrosymmetric vs. acentric.
//
// `reflections` are RAW (possibly redundant) observations; `cell` is the
// direct unit cell; `composition` is [{element, count}] for the WHOLE
// unit cell (see ins-parser.js / the formula-times-Z path in index.html).
function computeEStatistics(reflections, cell, composition, groupMatrices) {
  const merged = mergeReflections(reflections, groupMatrices);
  if (merged.length < 60) {
    return { error: `Only ${merged.length} symmetry-unique reflections — too few for a meaningful Wilson plot (need at least ~60).` };
  }

  const recip = reciprocalCellFromDirect(cell);
  const sValues = merged.map(r => sinThetaOverLambda(r.h, r.k, r.l, recip));

  const shells = binByResolution(merged, sValues);
  const wilsonPoints = [];
  const missingElements = new Set();

  for (const shell of shells) {
    const meanF2 = shell.reduce((a, p) => a + p.r.f2, 0) / shell.length;
    const meanS = shell.reduce((a, p) => a + p.s, 0) / shell.length;
    const { sum: sumF2calc, missing } = sumSquaredScatteringFactors(composition, meanS);
    missing.forEach(m => missingElements.add(m));
    if (sumF2calc <= 0 || meanF2 <= 0) continue; // ln() undefined; skip this shell
    wilsonPoints.push({ x: meanS * meanS, y: Math.log(meanF2 / sumF2calc), s: meanS, meanF2, n: shell.length });
  }

  if (missingElements.size > 0) {
    return { error: `No scattering-factor data for: ${[...missingElements].join(', ')}. Check the element composition.`, missingElements: [...missingElements] };
  }
  if (wilsonPoints.length < 4) {
    return { error: 'Not enough usable resolution shells for a Wilson plot (too few positive-intensity shells).' };
  }

  // ln(<Fo^2>/sum f^2) = -ln(k) - 2*B*s^2  (see module doc for the
  // derivation) => slope = -2B, intercept = -ln(k)
  const { slope, intercept } = linearRegression(wilsonPoints.map(p => ({ x: p.x, y: p.y })));
  const B = -slope / 2;
  const k = Math.exp(-intercept); // scales observed F² onto an absolute (electron²) scale

  let usable = 0, sumE2 = 0, sumAbsE2minus1 = 0, sumAbsE = 0, countEgt3 = 0;
  for (let i = 0; i < merged.length; i++) {
    const s = sValues[i];
    const { sum: sumF2calc } = sumSquaredScatteringFactors(composition, s);
    if (sumF2calc <= 0) continue;
    const expected = sumF2calc * Math.exp(-2 * B * s * s);
    if (expected <= 0) continue;
    const e2 = (k * merged[i].f2) / expected;
    if (!Number.isFinite(e2) || e2 < 0) continue;
    usable++;
    sumE2 += e2;
    sumAbsE2minus1 += Math.abs(e2 - 1);
    const eAbs = Math.sqrt(e2);
    sumAbsE += eAbs;
    if (eAbs > 3) countEgt3++;
  }

  if (usable < 60) {
    return { error: `Only ${usable} reflections had a usable E² value — too few to judge reliably.` };
  }

  const meanE2 = sumE2 / usable;
  const meanE2minus1 = sumAbsE2minus1 / usable;
  const meanE = sumAbsE / usable;
  const fractionEgt3 = countEgt3 / usable;

  const distCentro = Math.abs(meanE2minus1 - WILSON_REFERENCE.centrosymmetric.meanE2minus1);
  const distAcentric = Math.abs(meanE2minus1 - WILSON_REFERENCE.acentric.meanE2minus1);
  const verdict = distCentro < distAcentric ? 'centrosymmetric' : 'acentric';
  // how much closer to the winning reference than the losing one — a
  // crude but honest confidence signal: near 0 means "right on the
  // fence", closer to 1 means unambiguous.
  const margin = Math.abs(distCentro - distAcentric) / (distCentro + distAcentric);

  return {
    uniqueReflections: merged.length, usableReflections: usable,
    B, k, wilsonPoints,
    meanE2, meanE2minus1, meanE, fractionEgt3,
    verdict, distCentro, distAcentric, margin,
  };
}

// --- Wilson plot visualization --------------------------------------------
// Canvas scatter of ln(<Fo^2>/sum f^2) per resolution shell against s^2,
// with the fitted line overlaid — lets a person sanity-check the fit
// itself (a heavy-atom "hook" at low resolution, obvious curvature, or a
// kink from twinning are all visible here, not just in the summary
// numbers). Same plain-canvas approach as absences-viz.js: no new
// dependency, file://-compatible.
let _lastWilsonPlotArgs = null;

function renderWilsonPlot(container, wilsonPoints, B, k) {
  _lastWilsonPlotArgs = [container, wilsonPoints, B, k];
  const colors = getPlotColors();
  const cssW = 480, cssH = 320;
  const dpr = window.devicePixelRatio || 1;
  const canvas = document.createElement('canvas');
  canvas.width = cssW * dpr;
  canvas.height = cssH * dpr;
  canvas.style.width = cssW + 'px';
  canvas.style.height = cssH + 'px';
  container.innerHTML = '';
  container.appendChild(canvas);
  const ctx = canvas.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, cssW, cssH);
  ctx.fillStyle = colors.bg;
  ctx.fillRect(0, 0, cssW, cssH);

  const margin = { left: 50, right: 15, top: 15, bottom: 40 };
  const plotW = cssW - margin.left - margin.right;
  const plotH = cssH - margin.top - margin.bottom;

  const xs = wilsonPoints.map(p => p.x), ys = wilsonPoints.map(p => p.y);
  const xMin = 0, xMax = Math.max(...xs) * 1.05;
  const yMin = Math.min(...ys) - 0.3, yMax = Math.max(...ys) + 0.3;
  const px = x => margin.left + (x / xMax) * plotW;
  const py = y => margin.top + plotH - ((y - yMin) / (yMax - yMin)) * plotH;

  // axes + a handful of ticks (floats, so a simple fixed tick count beats
  // trying to force "nice" integers the way absences-viz.js does)
  ctx.strokeStyle = colors.grid;
  ctx.fillStyle = colors.text;
  ctx.font = '10px sans-serif';
  ctx.lineWidth = 1;
  const xTicks = 5, yTicks = 5;
  for (let i = 0; i <= xTicks; i++) {
    const xv = xMin + (xMax - xMin) * (i / xTicks);
    const x = px(xv);
    ctx.beginPath(); ctx.moveTo(x, margin.top); ctx.lineTo(x, margin.top + plotH); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(xv.toFixed(2), x, cssH - margin.bottom + 14);
  }
  for (let i = 0; i <= yTicks; i++) {
    const yv = yMin + (yMax - yMin) * (i / yTicks);
    const y = py(yv);
    ctx.beginPath(); ctx.moveTo(margin.left, y); ctx.lineTo(margin.left + plotW, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(yv.toFixed(2), margin.left - 6, y + 3);
  }

  // fitted line: y = -ln(k) - 2*B*x
  const intercept = -Math.log(k);
  const slope = -2 * B;
  ctx.strokeStyle = colors.fitLine;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(px(xMin), py(intercept + slope * xMin));
  ctx.lineTo(px(xMax), py(intercept + slope * xMax));
  ctx.stroke();

  // shell points, sized by how many reflections went into that shell
  const maxN = Math.max(...wilsonPoints.map(p => p.n));
  for (const p of wilsonPoints) {
    const r = 2 + 3 * Math.sqrt(p.n / maxN);
    ctx.fillStyle = colors.point;
    ctx.beginPath();
    ctx.arc(px(p.x), py(p.y), r, 0, 2 * Math.PI);
    ctx.fill();
  }

  ctx.fillStyle = colors.text;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  ctx.fillText('sin²(θ)/λ² (Å⁻²)', margin.left + plotW / 2, cssH - 6);
  ctx.save();
  ctx.translate(12, margin.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('ln(⟨Fo²⟩ / Σf²)', 0, 0);
  ctx.restore();
}

onPlotThemeChange(() => {
  if (_lastWilsonPlotArgs) renderWilsonPlot(..._lastWilsonPlotArgs);
});
