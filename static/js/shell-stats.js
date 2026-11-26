// shell-stats.js
// Merging statistics per resolution shell (shown on the "Laue Symmetry" tab):
// unique reflections, redundancy, <I/σ>, share of I > 2σ, and Rint / Rmeas /
// Rpim. Purely informational — it helps to judge where the data stop being
// useful (a resolution limit), but nothing here is applied to the .ins or
// the .hkl.
//
// Definitions (merged over symmetry equivalents under the given point group,
// see groupReflectionsByCanonicalHKL in laue-merge.js; only groups with at
// least two members enter the R values):
//   Rint  = Σ|I − <I>| / Σ I
//   Rmeas = Σ sqrt(n/(n−1)) Σ|I − <I>| / Σ I   (redundancy-independent)
//   Rpim  = Σ sqrt(1/(n−1)) Σ|I − <I>| / Σ I   (precision of the merged value)
// <I/σ> is the mean of I/σ over the MERGED reflections, with σ the standard
// error sqrt(Σσᵢ²)/n — the same σ that mergeReflections() reports.
// CC½ by the σ–τ method (Assmann, Brehm & Diederichs 2016, J. Appl. Cryst. 49,
// 1021; as in XDSCC12 and CrystFEL's partialator; see also the XDSwiki page
// "CC1/2"): with σ²_y = variance of the averaged intensities of the unique
// reflections of the shell and σ²_ε = mean over the reflections of the
// variance of their merged intensity (sample variance of the observations,
// divided by n/2),
//     CC½ = (σ²_y − ½ σ²_ε) / (σ²_y + ½ σ²_ε).
// No random split into two halves: the result is deterministic and unbiased,
// and it can become slightly negative when only noise is left. Unweighted
// (the σ values of the single observations are not used), only reflections
// measured at least twice enter.

function accumulateShell(chunk, sLo, sHi, possible, nAllowed = null) {
  let nTotal = 0, sumISig = 0, nISig = 0, nGt2 = 0;
  let numInt = 0, numMeas = 0, numPim = 0, den = 0;
  let sumMeanI = 0, sumMeanSig = 0; // R_sigma = sum sigma(<I>) / sum <I> over the merged reflections
  const ccLists = []; // intensities of each reflection measured at least twice, for CC½
  for (const { members } of chunk) {
    const n = members.length;
    nTotal += n;
    const sum = members.reduce((a, r) => a + r.f2, 0);
    const mean = sum / n;
    const sig = Math.sqrt(members.reduce((a, r) => a + r.sigma * r.sigma, 0)) / n;
    sumMeanI += mean; sumMeanSig += sig;
    if (sig > 0) {
      const x = mean / sig;
      sumISig += x; nISig++;
      if (x > 2) nGt2++;
    }
    if (n >= 2) {
      let dev = 0;
      for (const r of members) dev += Math.abs(r.f2 - mean);
      numInt += dev;
      numMeas += Math.sqrt(n / (n - 1)) * dev;
      numPim += Math.sqrt(1 / (n - 1)) * dev;
      den += sum;
      ccLists.push(members.map(r => r.f2));
    }
  }
  const nUnique = chunk.length;
  return {
    dMax: sLo > 0 ? 1 / (2 * sLo) : null, dMin: 1 / (2 * sHi),
    nUnique, nTotal, redundancy: nUnique ? nTotal / nUnique : null,
    possible: possible === null ? null : possible,
    completeness: possible ? Math.min(1, (nAllowed === null ? nUnique : nAllowed) / possible) : null,
    meanISigma: nISig ? sumISig / nISig : null,
    fractionGt2: nISig ? nGt2 / nISig : null,
    rInt: den > 0 ? numInt / den : null,
    rMeas: den > 0 ? numMeas / den : null,
    rPim: den > 0 ? numPim / den : null,
    rSigma: sumMeanI > 0 ? sumMeanSig / sumMeanI : null,
    cc12: ccLists.length >= 10 ? cc12SigmaTau(ccLists) : null,
    nCC: ccLists.length,
  };
}

// σ–τ CC½ from the lists of observations of the unique reflections (each with
// n ≥ 2 observations). Returns null if it cannot be formed (fewer than two
// reflections, or a non-positive denominator).
function cc12SigmaTau(lists) {
  const N = lists.length;
  if (N < 2) return null;
  let sx = 0, sxx = 0, sEps = 0;
  for (const xs of lists) {
    const n = xs.length;
    let sum = 0, sum2 = 0;
    for (const x of xs) { sum += x; sum2 += x * x; }
    const xbar = sum / n;
    sx += xbar; sxx += xbar * xbar;
    sEps += ((sum2 - sum * sum / n) / (n - 1)) / (n / 2);   // unbiased variance of the mean of the two half-data sets
  }
  const sigY2 = (sxx - sx * sx / N) / (N - 1);
  const sigEps2 = sEps / N;
  const den = sigY2 + 0.5 * sigEps2;
  return den > 0 ? (sigY2 - 0.5 * sigEps2) / den : null;
}

// All symmetry-unique hkl up to sMax that the lattice allows — the denominator
// of the completeness. A canonical hkl is its own representative, so counting
// the hkl that map onto themselves under the Laue group counts each unique
// reflection exactly once, without storing anything. Returns the sorted
// s = sinθ/λ values, or null if the cell/resolution is so large that the loop
// would take noticeably long.
function theoreticalUniqueS(cell, recip, sMax, groupMatrices, latticeTest, maxLoop = 2e7) {
  const hMax = Math.floor(2 * sMax * cell.a), kMax = Math.floor(2 * sMax * cell.b), lMax = Math.floor(2 * sMax * cell.c);
  if ((2 * hMax + 1) * (2 * kMax + 1) * (2 * lMax + 1) > maxLoop) return null;
  const out = [];
  const lim = sMax + 1e-9;
  for (let h = -hMax; h <= hMax; h++) {
    for (let k = -kMax; k <= kMax; k++) {
      for (let l = -lMax; l <= lMax; l++) {
        if (h === 0 && k === 0 && l === 0) continue;
        const s = sinThetaOverLambda(h, k, l, recip);
        if (s > lim) continue;
        if (latticeTest && !latticeTest({ h, k, l })) continue;
        const c = canonicalHKL(h, k, l, groupMatrices);
        if (c[0] === h && c[1] === k && c[2] === l) out.push(s);
      }
    }
  }
  return out.sort((a, b) => a - b);
}

// Shell boundaries shared by the merging table and the E² shell table, so
// that row n means the same resolution range in both. `sSorted` = s = sinθ/λ
// of the measured unique reflections, ascending. Shells have equal width in
// d*³, i.e. the same number of POSSIBLE unique reflections each. Without the
// theoretical list (huge cell) the boundaries fall back to equal numbers of
// MEASURED reflections.
// Returns { bounds (upper s of each shell), theory (sorted s list or null),
//           haveTheory, sMin, sMax }.
function computeShellBounds(sSorted, cell, groupMatrices, latticeTest = null, targetShells = 10) {
  const sMin = sSorted[0], sMax = sSorted[sSorted.length - 1];
  const recip = reciprocalCellFromDirect(cell);
  const theory = theoreticalUniqueS(cell, recip, sMax, groupMatrices, latticeTest);
  const haveTheory = !!theory && theory.length >= 40;
  const basis = haveTheory ? theory : sSorted;
  const nShells = Math.max(1, Math.min(targetShells, Math.floor(basis.length / 20)));
  const per = Math.ceil(basis.length / nShells);
  const bounds = [];
  for (let i = 0; i < nShells; i++) bounds.push(basis[Math.min((i + 1) * per, basis.length) - 1]);
  bounds[nShells - 1] = Math.max(bounds[nShells - 1], sMax);
  return { bounds, theory: haveTheory ? theory : null, haveTheory, sMin, sMax };
}

// `reflections` raw (possibly redundant); `groupMatrices` from generateGroup();
// `latticeTest` (optional, from filterByLattice) = which hkl the lattice
// centering allows. The Completeness column shows where the data have gaps.
function computeShellStatistics(reflections, cell, groupMatrices, targetShells = 10, latticeTest = null) {
  const groups = groupReflectionsByCanonicalHKL(reflections, groupMatrices);
  const recip = reciprocalCellFromDirect(cell);
  const uniques = [];
  for (const g of groups.values()) {
    const s = sinThetaOverLambda(g.h, g.k, g.l, recip);
    if (s > 0) uniques.push({ s, members: g.members, allowed: !latticeTest || latticeTest({ h: g.h, k: g.k, l: g.l }) });
  }
  if (uniques.length < 40) {
    return { error: `Only ${uniques.length} unique reflections — too few for resolution shells.` };
  }
  uniques.sort((a, b) => a.s - b.s);
  const info = computeShellBounds(uniques.map(u => u.s), cell, groupMatrices, latticeTest, targetShells);
  const { bounds, theory, haveTheory, sMax } = info;

  const chunks = bounds.map(() => []);
  let idx = 0;
  for (const u of uniques) {
    while (idx < bounds.length - 1 && u.s > bounds[idx]) idx++;
    chunks[idx].push(u);
  }
  const possibleIn = (lo, hi) => {
    if (!haveTheory) return null;
    let c = 0;
    for (const s of theory) { if (s > hi) break; if (s > lo) c++; }
    return c;
  };
  const shells = chunks.map((chunk, i) => {
    const lo = i === 0 ? 0 : bounds[i - 1];
    const sLo = i === 0 ? uniques[0].s : lo;
    return accumulateShell(chunk, sLo, bounds[i], possibleIn(i === 0 ? -1 : lo, bounds[i]), chunk.filter(u => u.allowed).length);
  });
  const overall = accumulateShell(uniques, uniques[0].s, sMax, haveTheory ? theory.length : null, uniques.filter(u => u.allowed).length);
  const cutoffIndex = shells.findIndex(sh => sh.meanISigma !== null && sh.meanISigma < 2);
  return { shells, overall, cutoffIndex, haveTheory };
}

function renderShellStatistics(container, stats, info) {
  if (stats.error) { container.innerHTML = `<p><small>${stats.error}</small></p>`; return; }
  const pct = x => (x === null ? '—' : (x * 100).toFixed(1) + ' %');
  const num = (x, d = 1) => (x === null ? '—' : x.toFixed(d));
  const compCell = sh => {
    if (sh.completeness === null) return '<td>—</td>';
    const low = sh.completeness < 0.95;
    return `<td${low ? ' style="color:#b60;font-weight:600;"' : ''} title="${sh.nUnique} of ${sh.possible} possible unique reflections">${(sh.completeness * 100).toFixed(1)} %</td>`;
  };
  const ccCell = sh => sh.cc12 === null ? '<td>—</td>'
    : `<td${sh.cc12 < 0.3 ? ' style="color:#b60;font-weight:600;"' : ''} title="${sh.nCC} reflections measured at least twice">${(sh.cc12 * 100).toFixed(1)} %</td>`;
  const row = (label, sh, cls) => `
    <tr${cls ? ` class="${cls}"` : ''}>
      <td>${label}</td><td>${sh.nUnique}</td>${compCell(sh)}<td>${sh.redundancy === null ? '—' : sh.redundancy.toFixed(2)}</td>
      <td>${num(sh.meanISigma)}</td><td>${pct(sh.fractionGt2)}</td>${ccCell(sh)}
      <td>${pct(sh.rInt)}</td><td>${pct(sh.rMeas)}</td><td>${pct(sh.rPim)}</td><td>${pct(sh.rSigma)}</td>
    </tr>`;
  const dLabel = sh => `${sh.dMax === null ? '∞' : sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)} Å`;
  const rows = stats.shells.map((sh, i) =>
    row(dLabel(sh), sh, i === stats.cutoffIndex ? 'shellCut' : '')).join('');
  const total = row(`all (${dLabel(stats.overall)})`, stats.overall, 'shellTotal');

  let cutNote;
  if (stats.cutoffIndex >= 0) {
    const sh = stats.shells[stats.cutoffIndex];
    cutNote = `<p><small>⟨I/σ⟩ drops below 2 in the highlighted shell (${dLabel(sh)}). A common rule of thumb puts the resolution limit at or just before such a shell.
      This is shown for orientation only — nothing is cut or written to the .ins/.hkl.</small></p>`;
  } else {
    cutNote = `<p><small>⟨I/σ⟩ stays above 2 in every shell — no obvious resolution cut-off in this data set. Nothing is cut or written to the .ins/.hkl.</small></p>`;
  }
  const ccCutIdx = stats.shells.findIndex(sh => sh.cc12 !== null && sh.cc12 < 0.3);
  if (ccCutIdx >= 0) {
    cutNote += `<p><small>CC<sub>½</sub> falls below 30 % from the shell ${dLabel(stats.shells[ccCutIdx])} on — there the two halves of the data no longer agree
      with each other, whatever the σ values say. Many programs use this (about 0.3 – 0.5) rather than ⟨I/σ⟩ ≈ 2 as the resolution criterion. Orange cells in the table: CC<sub>½</sub> &lt; 30 %.</small></p>`;
  }
  const redNote = stats.overall.rInt === null
    ? `<p style="color:#a60;"><small>⚠ No repeated observations: R values are undefined (data already merged to a unique set).</small></p>` : '';

  container.innerHTML = `
    <p><small>Merged under Laue class ${formatPG(info.laue)}${info.source === 'cif' ? ' (from the CIF space group)' : info.source === 'merging' ? ' (confirmed by merging)' : ' (assumed — highest-symmetry candidate for this cell, not confirmed)'}; ${info.lattice}.
      Shells have equal width in d*³ (the same number of possible unique reflections each), from low to high resolution.
      Where the mean intensity approaches zero the R values grow without limit — use ⟨I/σ⟩ there.</small></p>
    <table class="shellTable">
      <tr><th>d range</th><th>Unique</th><th>Complete</th><th>Redundancy</th><th>⟨I/σ⟩</th><th>I &gt; 2σ</th><th>CC<sub>½</sub></th>
          <th>R<sub>int</sub></th><th>R<sub>meas</sub></th><th>R<sub>pim</sub></th><th>R<sub>σ</sub></th></tr>
      ${rows}${total}
    </table>
    ${redNote}${cutNote}
    <div class="plotContainer shellCurves"><div id="shellCurvesA"></div><div id="shellCurvesB"></div></div>
    <details class="shellLegend">
      <summary>What do these columns mean?</summary>
      <ul>
        <li><strong>d range</strong> — resolution shell in Å (low to high resolution).</li>
        <li><strong>Unique</strong> / <strong>Redundancy</strong> — symmetry-unique reflections measured in the shell, and how often each was measured on average.</li>
        <li><strong>Complete</strong> — measured unique reflections divided by the unique reflections the cell, lattice centering and Laue class allow in that shell (hover for the counts). Below 95 % is shown in orange: gaps in the data (beamstop, detector gaps, short scan). The bottom row is the completeness up to the highest measured resolution. Reflections that the lattice centering forbids, or that the absence scan found to be systematically absent (screw axes, glides), count neither as possible nor as measured — so files already cleaned of absences (deposited CIF/.hkl data) are not penalised. If a systematic absence was NOT recognised by the scan, the value is too low.</li>
        <li><strong>⟨I/σ⟩</strong> — mean signal-to-noise of the merged reflections (intensity divided by its standard uncertainty σ). Higher is better; a common rule of thumb puts the resolution limit where it falls below about 2.</li>
        <li><strong>I &gt; 2σ</strong> — share of merged reflections whose intensity is more than twice its uncertainty ("observed" reflections).</li>
        <li><strong>CC<sub>½</sub></strong> — how well the data agree with themselves: the variance of the averaged intensities against the variance of their merging error, CC<sub>½</sub> = (σ²<sub>y</sub> − ½σ²<sub>ε</sub>) / (σ²<sub>y</sub> + ½σ²<sub>ε</sub>) (σ–τ method, Assmann et al. 2016, as in XDSCC12 and CrystFEL; no random split). Only reflections measured at least twice. Close to 100 % = reproducible; it falls towards 0 where only noise is left (it can even turn slightly negative). Does not use the σ values. Orange below 30 %.</li>
        <li><strong>R<sub>int</sub></strong> — spread of equivalent measurements around their mean, Σ|I − ⟨I⟩| / ΣI. It grows with redundancy even when the data get better.</li>
        <li><strong>R<sub>meas</sub></strong> — R<sub>int</sub> corrected for redundancy (factor √(n/(n−1))): the quality of a single measurement, comparable between data sets.</li>
        <li><strong>R<sub>pim</sub></strong> — R<sub>int</sub> corrected with √(1/(n−1)): the precision of the merged (averaged) intensities that go into refinement; always smaller than R<sub>meas</sub>.</li>
        <li><strong>R<sub>σ</sub></strong> — Σσ(⟨I⟩) / Σ⟨I⟩ over the merged reflections: the noise level of the data from the σ values alone (SHELXL's R(sigma)). Defined for merged data as well.</li>
      </ul>
      <p><small>R values are only defined for reflections measured at least twice, and they diverge where the mean intensity approaches zero — rely on ⟨I/σ⟩ there.</small></p>
    </details>`;
  renderShellCurves(container, stats);
}


// --- Curves under the table ---------------------------------------------------
// Plot A (percent scale): completeness, I > 2σ, CC½ and R_meas (clipped at
// 100 %) against (sin θ/λ)² of the shell. Plot B: ⟨I/σ⟩ on a log axis with the
// ⟨I/σ⟩ = 2 line. Gaps in a series (R, CC½ undefined) just break the line.
const SHELL_CURVE_COLORS = { complete: '#2a7bd4', gt2: '#2a9d5c', cc12: '#d9822b', rMeas: '#c0392b', isig: '#2a7bd4' };
let _lastShellCurves = null;

// The two curve plots (also drawn by the report, report-writer.js).
function shellCurveSpecs(stats) {
  const clip = x => (x === null ? null : Math.min(100, x * 100));
  const specA = {
    caption: '<strong>Shell statistics (percent)</strong>',
    yMin: 0, yMax: 100, log: false, yLabel: '%',
    series: [
      { label: 'completeness', color: SHELL_CURVE_COLORS.complete, vals: sh => clip(sh.completeness) },
      { label: 'I > 2σ', color: SHELL_CURVE_COLORS.gt2, vals: sh => clip(sh.fractionGt2) },
      { label: 'CC½', color: SHELL_CURVE_COLORS.cc12, vals: sh => clip(sh.cc12) },
      { label: 'Rmeas (≤ 100)', color: SHELL_CURVE_COLORS.rMeas, vals: sh => clip(sh.rMeas) },
    ],
  };
  const maxI = Math.max(...stats.shells.map(sh => sh.meanISigma || 0), 4);
  const specB = {
    caption: '<strong>⟨I/σ⟩ per shell</strong>',
    yMin: 0.5, yMax: Math.pow(10, Math.ceil(Math.log10(maxI))), log: true, yLabel: '⟨I/σ⟩',
    refLine: { y: 2, label: '⟨I/σ⟩ = 2' },
    series: [{ label: '⟨I/σ⟩', color: SHELL_CURVE_COLORS.isig, vals: sh => (sh.meanISigma === null ? null : Math.max(0.5, sh.meanISigma)) }],
  };
  return [specA, specB];
}

function renderShellCurves(container, stats) {
  _lastShellCurves = [container, stats];
  const a = container.querySelector('#shellCurvesA'), b = container.querySelector('#shellCurvesB');
  if (!a || !b || !stats.shells || stats.shells.length < 2) return;
  const [specA, specB] = shellCurveSpecs(stats);
  for (const [el, spec] of [[a, specA], [b, specB]]) {
    const canvas = document.createElement('canvas');
    el.innerHTML = ''; el.appendChild(canvas);
    drawShellCurve(canvas, stats, spec);
    attachPlotZoom(canvas, spec.caption, 480, (cv, zoom) => drawShellCurve(cv, stats, spec, zoom));
  }
}

function drawShellCurve(canvas, stats, spec, zoom = 1) {
  const colors = getPlotColors();
  const cssW = 480, cssH = 260;
  const g = setupCanvasForDPR(canvas, cssW, cssH, colors, zoom);
  const m = { left: 50, right: 15, top: 14, bottom: 40 };
  const pw = cssW - m.left - m.right, ph = cssH - m.top - m.bottom;
  const xs = stats.shells.map(sh => {
    const sHi = 1 / (2 * sh.dMin), sLo = sh.dMax === null ? 0 : 1 / (2 * sh.dMax);
    const mid = (sLo + sHi) / 2;
    return mid * mid;
  });
  const xMax = Math.max(...xs) * 1.05;
  const px = x => m.left + (x / xMax) * pw;
  const ty = v => (spec.log ? Math.log10(v) : v);
  const lo = ty(spec.yMin), hi = ty(spec.yMax);
  const py = v => m.top + ph - ((ty(v) - lo) / (hi - lo)) * ph;

  g.lineWidth = 1; g.strokeStyle = colors.grid; g.fillStyle = colors.text; g.font = '10px sans-serif';
  const yTicks = spec.log ? [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500].filter(v => v >= spec.yMin && v <= spec.yMax) : [0, 20, 40, 60, 80, 100];
  for (const v of yTicks) {
    g.beginPath(); g.moveTo(m.left, py(v)); g.lineTo(m.left + pw, py(v)); g.stroke();
    g.textAlign = 'right'; g.fillText(String(v), m.left - 6, py(v) + 3);
  }
  const xStep = xMax > 0.3 ? 0.1 : 0.05;
  for (let x = 0; x <= xMax; x += xStep) {
    g.beginPath(); g.moveTo(px(x), m.top); g.lineTo(px(x), m.top + ph); g.stroke();
    g.textAlign = 'center'; g.fillText(x.toFixed(2), px(x), cssH - m.bottom + 14);
  }
  if (spec.refLine) {
    g.strokeStyle = colors.axis; g.setLineDash([5, 3]);
    g.beginPath(); g.moveTo(m.left, py(spec.refLine.y)); g.lineTo(m.left + pw, py(spec.refLine.y)); g.stroke(); g.setLineDash([]);
    g.fillStyle = colors.text; g.textAlign = 'right'; g.fillText(spec.refLine.label, m.left + pw - 4, py(spec.refLine.y) - 4);
  }
  for (const ser of spec.series) {
    g.strokeStyle = ser.color; g.fillStyle = ser.color; g.lineWidth = 1.8;
    let pen = false;
    stats.shells.forEach((sh, i) => {
      const v = ser.vals(sh);
      if (v === null) { pen = false; return; }
      if (!pen) { g.beginPath(); g.moveTo(px(xs[i]), py(v)); pen = true; } else { g.lineTo(px(xs[i]), py(v)); }
      g.stroke();
      g.beginPath(); g.arc(px(xs[i]), py(v), 3, 0, 2 * Math.PI); g.fill();
      g.beginPath(); g.moveTo(px(xs[i]), py(v));
    });
  }
  if (spec.series.length > 1) {   // legend, bottom left
    g.font = '10px sans-serif'; g.textAlign = 'left';
    spec.series.forEach((ser, i) => {
      const y = m.top + ph - 8 - (spec.series.length - 1 - i) * 13;
      g.fillStyle = ser.color; g.fillRect(m.left + 8, y - 7, 14, 3);
      g.fillStyle = colors.text; g.fillText(ser.label, m.left + 28, y - 2);
    });
  }
  g.fillStyle = colors.text; g.font = '11px sans-serif'; g.textAlign = 'center';
  g.fillText('(sin θ/λ)² (Å⁻²), shell mid', m.left + pw / 2, cssH - 6);
  g.save(); g.translate(12, m.top + ph / 2); g.rotate(-Math.PI / 2); g.fillText(spec.yLabel, 0, 0); g.restore();
}

onPlotThemeChange(() => { if (_lastShellCurves) renderShellCurves(..._lastShellCurves); });
