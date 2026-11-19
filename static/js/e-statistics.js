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

// --- Reference cumulative distributions of z = E² (intensity) --------------
// N(z) = fraction of reflections with E² <= z (Howells, Phillips & Rogers
// 1950). Ideal, untwinned data:
//   acentric:                  N(z) = 1 - exp(-z)
//   centrosymmetric:           N(z) = erf(sqrt(z/2))
//   acentric, perfect twin:    N(z) = 1 - (1 + 2z) exp(-2z)   (hemihedral, α = 0.5)
// The twinned curve looks "more centrosymmetric" at high z but has far
// fewer weak reflections than the centric curve — that shape difference is
// what separates a twinned acentric crystal from a truly centric one.

// Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7) — plenty for a plot and a
// Kolmogorov distance; Math.erf does not exist in JavaScript.
function erfApprox(x) {
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * x);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return sign * y;
}
const NZ_REFERENCE = {
  centrosymmetric: z => erfApprox(Math.sqrt(Math.max(0, z) / 2)),
  acentric: z => 1 - Math.exp(-Math.max(0, z)),
  twinnedAcentric: z => 1 - (1 + 2 * Math.max(0, z)) * Math.exp(-2 * Math.max(0, z)),
};

// Largest vertical gap between the empirical step function of a SORTED
// array and a reference CDF (Kolmogorov–Smirnov distance D).
function ksDistance(sorted, cdf) {
  const n = sorted.length;
  let d = 0;
  for (let i = 0; i < n; i++) {
    const F = cdf(sorted[i]);
    d = Math.max(d, Math.abs((i + 1) / n - F), Math.abs(i / n - F));
  }
  return d;
}

// Builds the N(z) summary from the per-reflection E² values.
function computeNZ(e2Values, negClipped) {
  const sorted = e2Values.slice().sort((a, b) => a - b);
  const n = sorted.length;
  const dCentro = ksDistance(sorted, NZ_REFERENCE.centrosymmetric);
  const dAcentric = ksDistance(sorted, NZ_REFERENCE.acentric);
  const dTwin = ksDistance(sorted, NZ_REFERENCE.twinnedAcentric);
  // empirical curve on a fixed grid z = 0 ... NZ_ZMAX (pointer walk, O(n + grid))
  const curve = [];
  let idx = 0;
  const steps = 200;
  for (let g = 0; g <= steps; g++) {
    const z = NZ_ZMAX * g / steps;
    while (idx < n && sorted[idx] <= z) idx++;
    curve.push({ z, n: idx / n });
  }
  const best = Math.min(dCentro, dAcentric, dTwin);
  const closest = best === dTwin ? 'acentric, twinned?' : (best === dCentro ? 'centrosymmetric' : 'acentric');
  return { n, negClipped, dCentro, dAcentric, dTwin, curve, closest };
}
const NZ_ZMAX = 4;

// Systematic absences caused by the LATTICE (centering) are not a feature
// of the structure: those reflections are zero by construction, so mixing
// them into Wilson/E² statistics distorts every distribution (a C-centred
// acentric structure with the absent reflections left in came out as
// "centrosymmetric" in testing). The centering is read from the Hall symbol
// of the currently selected space-group candidate — only the lattice
// letter is used, never the (circular) screw/glide conditions.
function filterByLattice(reflections, hall) {
  if (!hall) return { reflections, note: 'no space-group candidate yet — all reflections used', removed: 0, test: null };
  const m = String(hall).trim().replace(/^-/, '').match(/^[A-Za-z]/);
  const letter = m ? m[0].toUpperCase() : 'P';
  const mod = (x, n) => ((x % n) + n) % n;
  let test = null, label = letter;
  switch (letter) {
    case 'A': test = r => mod(r.k + r.l, 2) === 0; break;
    case 'B': test = r => mod(r.h + r.l, 2) === 0; break;
    case 'C': test = r => mod(r.h + r.k, 2) === 0; break;
    case 'I': test = r => mod(r.h + r.k + r.l, 2) === 0; break;
    case 'F': test = r => mod(r.h + r.k, 2) === 0 && mod(r.h + r.l, 2) === 0; break;
    case 'R': {
      // obverse (-h+k+l = 3n) or reverse (h-k+l = 3n) setting: pick whichever
      // one the data actually obey (higher mean intensity of allowed reflections)
      const obv = r => mod(-r.h + r.k + r.l, 3) === 0;
      const rev = r => mod(r.h - r.k + r.l, 3) === 0;
      const meanOf = t => { let s = 0, c = 0; for (const r of reflections) if (t(r)) { s += r.f2; c++; } return c ? s / c : 0; };
      const mo = meanOf(obv), mr = meanOf(rev);
      if (mo > 1.5 * mr) { test = obv; label = 'R (obverse)'; }
      else if (mr > 1.5 * mo) { test = rev; label = 'R (reverse)'; }
      break;
    }
    default: break; // P (and anything unrecognised): nothing to remove
  }
  if (!test) return { reflections, note: letter === 'R' ? 'lattice R: obverse/reverse setting unclear — all reflections used' : 'primitive lattice — no centering absences to remove', removed: 0, test: null };
  const kept = reflections.filter(test);
  const removed = reflections.length - kept.length;
  return { reflections: kept, note: `lattice ${label}: ${removed} centering-absent reflection(s) excluded before merging`, removed, test };
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
// --- E² statistics per resolution shell --------------------------------------
// The overall <|E²−1|> mixes all resolutions. A dominant heavy atom (strong at
// low resolution only) or pseudosymmetry bends the distribution in part of
// the range, so the same number is also worked out shell by shell. Within each
// shell E² is rescaled to <E²> = 1 — otherwise the global Wilson fit (one
// straight line through all shells) would leak into the comparison. The raw
// <E²> of the global fit is kept to show where that fit misses.
const E_SHELL_MIN = 100;   // reflections per shell below which the scatter gets too large
const E_SHELL_MAX = 8;     // only used by the fallback (no theoretical shell list)

// `shellInfo` = result of computeShellBounds(): the same shell boundaries as
// the merging table (row n = same resolution range). Shells with fewer than
// E_SHELL_MIN usable reflections are still listed but marked `few` and left
// out of the agreement check. Without shellInfo the shells hold equal numbers
// of measured reflections instead.
function computeEShells(pairs, shellInfo = null) {
  const n = pairs.length;
  const sorted = pairs.slice().sort((a, b) => a.s - b.s);
  let chunks, sLo = [];
  if (shellInfo && shellInfo.haveTheory) {
    const { bounds, sMin } = shellInfo;
    chunks = bounds.map(() => []);
    let idx = 0;
    for (const p of sorted) {
      while (idx < bounds.length - 1 && p.s > bounds[idx]) idx++;
      chunks[idx].push(p);
    }
    sLo = bounds.map((b, i) => (i === 0 ? sMin : bounds[i - 1]));
  } else {
    const nShells = Math.min(E_SHELL_MAX, Math.floor(n / E_SHELL_MIN));
    if (nShells < 2) return { tooFew: true, n };
    const per = Math.ceil(n / nShells);
    chunks = [];
    for (let i = 0; i < n; i += per) chunks.push(sorted.slice(i, i + per));
  }
  const mid = (WILSON_REFERENCE.centrosymmetric.meanE2minus1 + WILSON_REFERENCE.acentric.meanE2minus1) / 2;
  const shells = [];
  chunks.forEach((chunk, i) => {
    const m = chunk.length;
    if (m === 0) return;
    const meanE2raw = chunk.reduce((a, p) => a + p.e2, 0) / m;
    if (!(meanE2raw > 0)) return;
    let sumAbs = 0, gt3 = 0, sumS2 = 0;
    for (const p of chunk) {
      const e2 = p.e2 / meanE2raw;
      sumAbs += Math.abs(e2 - 1);
      if (Math.sqrt(e2) > 3) gt3++;
      sumS2 += p.s * p.s;
    }
    const meanAbs = sumAbs / m;
    const lo = sLo.length ? sLo[i] : chunk[0].s;
    const hi = shellInfo && shellInfo.haveTheory ? shellInfo.bounds[i] : chunk[m - 1].s;
    shells.push({
      n: m, meanE2raw, meanAbs, fractionEgt3: gt3 / m,
      s2: sumS2 / m,
      dMax: lo > 0 ? 1 / (2 * lo) : null, dMin: 1 / (2 * hi),
      few: m < E_SHELL_MIN,
      verdict: meanAbs < mid ? 'acentric' : 'centrosymmetric',
    });
  });
  const used = shells.filter(sh => !sh.few);
  if (used.length < 2) return { tooFew: true, n };
  const verdicts = new Set(used.map(sh => sh.verdict));
  return { shells, consistent: verdicts.size === 1, nFew: shells.length - used.length, n };
}

// Plot: shell-wise <|E²−1|> against (sin θ/λ)², with the two reference values.
let _lastEShellPlotArgs = null;

function renderEShellPlot(container, eShells) {
  _lastEShellPlotArgs = [container, eShells];
  const canvas = document.createElement('canvas');
  container.innerHTML = '';
  container.appendChild(canvas);
  drawEShellPlot(canvas, eShells);
  attachPlotZoom(canvas, '<strong>E² statistics by resolution shell</strong>', 480,
    (cv, zoom) => drawEShellPlot(cv, eShells, zoom));
}

function drawEShellPlot(canvas, eShells, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 480, cssH = 300;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  const margin = { left: 50, right: 15, top: 15, bottom: 40 };
  const plotW = cssW - margin.left - margin.right;
  const plotH = cssH - margin.top - margin.bottom;
  const sh = eShells.shells;
  const xMax = Math.max(...sh.map(p => p.s2)) * 1.08;
  const yMin = 0.5, yMax = 1.2;
  const px = x => margin.left + (x / xMax) * plotW;
  const py = y => margin.top + plotH - ((y - yMin) / (yMax - yMin)) * plotH;

  ctx.strokeStyle = colors.grid; ctx.fillStyle = colors.text; ctx.font = '10px sans-serif'; ctx.lineWidth = 1;
  for (let y = 0.5; y <= 1.2001; y += 0.1) {
    ctx.beginPath(); ctx.moveTo(margin.left, py(y)); ctx.lineTo(margin.left + plotW, py(y)); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(y.toFixed(1), margin.left - 6, py(y) + 3);
  }
  const xStep = xMax > 0.3 ? 0.1 : 0.05;
  for (let x = 0; x <= xMax; x += xStep) {
    ctx.beginPath(); ctx.moveTo(px(x), margin.top); ctx.lineTo(px(x), margin.top + plotH); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(x.toFixed(2), px(x), cssH - margin.bottom + 14);
  }

  const refLine = (y, color, label) => {
    ctx.strokeStyle = color; ctx.lineWidth = 1.5; ctx.setLineDash([5, 3]);
    ctx.beginPath(); ctx.moveTo(margin.left, py(y)); ctx.lineTo(margin.left + plotW, py(y)); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = color; ctx.textAlign = 'right';
    ctx.fillText(label, margin.left + plotW - 4, py(y) - 4);
  };
  refLine(WILSON_REFERENCE.centrosymmetric.meanE2minus1, colors.present, 'centrosymmetric 0.968');
  refLine(WILSON_REFERENCE.acentric.meanE2minus1, colors.fitLine, 'acentric 0.736');

  const y = v => Math.max(yMin, Math.min(yMax, v)); // keep outliers on the plot
  ctx.strokeStyle = colors.point; ctx.lineWidth = 1.5;
  ctx.beginPath();
  sh.filter(p => !p.few).forEach((p, i) => { i === 0 ? ctx.moveTo(px(p.s2), py(y(p.meanAbs))) : ctx.lineTo(px(p.s2), py(y(p.meanAbs))); });
  ctx.stroke();
  ctx.fillStyle = colors.point;
  sh.forEach(p => {   // hollow circle = fewer than E_SHELL_MIN reflections (not used for the comparison)
    ctx.beginPath(); ctx.arc(px(p.s2), py(y(p.meanAbs)), 4, 0, 2 * Math.PI);
    if (p.few) ctx.stroke(); else ctx.fill();
  });

  ctx.fillStyle = colors.text; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('(sin θ/λ)² (Å⁻²), shell mean', margin.left + plotW / 2, cssH - 6);
  ctx.save(); ctx.translate(12, margin.top + plotH / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('⟨|E²−1|⟩ per shell', 0, 0); ctx.restore();
}

function computeEStatistics(reflections, cell, composition, groupMatrices, latticeTest = null) {
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
  const nzValues = []; // E² per reflection for the N(z) test
  const eBySin = [];   // (s, E²) of the usable reflections, for the per-shell statistics
  let negClipped = 0;
  for (let i = 0; i < merged.length; i++) {
    const s = sValues[i];
    const { sum: sumF2calc } = sumSquaredScatteringFactors(composition, s);
    if (sumF2calc <= 0) continue;
    const expected = sumF2calc * Math.exp(-2 * B * s * s);
    if (expected <= 0) continue;
    const e2 = (k * merged[i].f2) / expected;
    if (!Number.isFinite(e2)) continue;
    // Negative net intensities are real WEAK reflections: for N(z) they count
    // as z = 0 instead of being dropped (dropping them would remove exactly
    // the low-z end of the distribution that the test relies on). The other
    // statistics below keep their original behaviour and skip them.
    if (e2 < 0) { nzValues.push(0); negClipped++; continue; }
    nzValues.push(e2);
    eBySin.push({ s, e2 });
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
    nz: computeNZ(nzValues, negClipped),
    eShells: computeEShells(eBySin, typeof computeShellBounds === 'function'
      ? computeShellBounds(sValues.filter(v => v > 0).sort((a, b) => a - b), cell, groupMatrices, latticeTest) : null),
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
  const canvas = document.createElement('canvas');
  container.innerHTML = '';
  container.appendChild(canvas);
  drawWilsonPlot(canvas, wilsonPoints, B, k);
  // click-to-enlarge (plot-zoom.js): re-draws the same plot at a larger zoom
  attachPlotZoom(canvas, '<strong>Wilson plot</strong>', 480,
    (cv, zoom) => drawWilsonPlot(cv, wilsonPoints, B, k, zoom));
}

// `zoom` > 1 draws the same plot larger (see setupCanvasForDPR in absences-viz.js).
function drawWilsonPlot(canvas, wilsonPoints, B, k, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 480, cssH = 320;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);

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
  ctx.fillText('(sin θ/λ)² (Å⁻²)', margin.left + plotW / 2, cssH - 6);
  ctx.save();
  ctx.translate(12, margin.top + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText('ln(⟨Fo²⟩ / Σf²), unscaled', 0, 0);
  ctx.restore();
}

// --- N(z) plot --------------------------------------------------------------
// Cumulative distribution of E² against the three ideal curves. Observed
// curve = solid blue; centrosymmetric = green; acentric = red; perfectly
// twinned acentric = dashed grey.
let _lastNZPlotArgs = null;

function renderNZPlot(container, nz) {
  _lastNZPlotArgs = [container, nz];
  const canvas = document.createElement('canvas');
  container.innerHTML = '';
  container.appendChild(canvas);
  drawNZPlot(canvas, nz);
  attachPlotZoom(canvas, '<strong>N(z) test</strong>', 480,
    (cv, zoom) => drawNZPlot(cv, nz, zoom));
}

function drawNZPlot(canvas, nz, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 480, cssH = 320;
  const ctx = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  const margin = { left: 50, right: 15, top: 15, bottom: 40 };
  const plotW = cssW - margin.left - margin.right;
  const plotH = cssH - margin.top - margin.bottom;
  const px = z => margin.left + (z / NZ_ZMAX) * plotW;
  const py = y => margin.top + plotH - y * plotH;

  ctx.strokeStyle = colors.grid; ctx.fillStyle = colors.text; ctx.font = '10px sans-serif'; ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i++) {            // x: z = 0 ... 4 in steps of 1
    const x = px(i);
    ctx.beginPath(); ctx.moveTo(x, margin.top); ctx.lineTo(x, margin.top + plotH); ctx.stroke();
    ctx.textAlign = 'center'; ctx.fillText(String(i), x, cssH - margin.bottom + 14);
  }
  for (let i = 0; i <= 5; i++) {            // y: N(z) = 0 ... 1 in steps of 0.2
    const y = py(i / 5);
    ctx.beginPath(); ctx.moveTo(margin.left, y); ctx.lineTo(margin.left + plotW, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText((i / 5).toFixed(1), margin.left - 6, y + 3);
  }

  const drawCurve = (fn, color, width, dash) => {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash || []);
    ctx.beginPath();
    for (let g = 0; g <= 200; g++) {
      const z = NZ_ZMAX * g / 200;
      g === 0 ? ctx.moveTo(px(z), py(fn(z))) : ctx.lineTo(px(z), py(fn(z)));
    }
    ctx.stroke(); ctx.setLineDash([]);
  };
  drawCurve(NZ_REFERENCE.centrosymmetric, colors.present, 1.5);
  drawCurve(NZ_REFERENCE.acentric, colors.fitLine, 1.5);
  drawCurve(NZ_REFERENCE.twinnedAcentric, colors.axis, 1.5, [4, 3]);

  ctx.strokeStyle = colors.point; ctx.lineWidth = 2.5;
  ctx.beginPath();
  nz.curve.forEach((p, i) => { i === 0 ? ctx.moveTo(px(p.z), py(p.n)) : ctx.lineTo(px(p.z), py(p.n)); });
  ctx.stroke();

  // legend, bottom right (the curves all sit in the upper left)
  const items = [['observed', colors.point, 2.5, null], ['centrosymmetric', colors.present, 1.5, null],
                 ['acentric', colors.fitLine, 1.5, null], ['acentric, perfect twin', colors.axis, 1.5, [4, 3]]];
  ctx.font = '10px sans-serif'; ctx.textAlign = 'left';
  items.forEach(([label, col, w, dash], i) => {
    const y = margin.top + plotH - 12 - (items.length - 1 - i) * 14;
    const x0 = margin.left + plotW - 135;
    ctx.strokeStyle = col; ctx.lineWidth = w; ctx.setLineDash(dash || []);
    ctx.beginPath(); ctx.moveTo(x0, y - 3); ctx.lineTo(x0 + 22, y - 3); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = colors.text; ctx.fillText(label, x0 + 28, y);
  });

  ctx.fillStyle = colors.text; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
  ctx.fillText('z = E²', margin.left + plotW / 2, cssH - 6);
  ctx.save(); ctx.translate(12, margin.top + plotH / 2); ctx.rotate(-Math.PI / 2);
  ctx.fillText('N(z) = fraction with E² ≤ z', 0, 0); ctx.restore();
}

onPlotThemeChange(() => {
  if (_lastWilsonPlotArgs) renderWilsonPlot(..._lastWilsonPlotArgs);
  if (_lastNZPlotArgs) renderNZPlot(..._lastNZPlotArgs);
  if (_lastEShellPlotArgs) renderEShellPlot(..._lastEShellPlotArgs);
});
