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

function accumulateShell(chunk, sLo, sHi, possible, nAllowed = null) {
  let nTotal = 0, sumISig = 0, nISig = 0, nGt2 = 0;
  let numInt = 0, numMeas = 0, numPim = 0, den = 0;
  for (const { members } of chunk) {
    const n = members.length;
    nTotal += n;
    const sum = members.reduce((a, r) => a + r.f2, 0);
    const mean = sum / n;
    const sig = Math.sqrt(members.reduce((a, r) => a + r.sigma * r.sigma, 0)) / n;
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
  };
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
  const row = (label, sh, cls) => `
    <tr${cls ? ` class="${cls}"` : ''}>
      <td>${label}</td><td>${sh.nUnique}</td>${compCell(sh)}<td>${sh.redundancy === null ? '—' : sh.redundancy.toFixed(2)}</td>
      <td>${num(sh.meanISigma)}</td><td>${pct(sh.fractionGt2)}</td>
      <td>${pct(sh.rInt)}</td><td>${pct(sh.rMeas)}</td><td>${pct(sh.rPim)}</td>
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
  const redNote = stats.overall.rInt === null
    ? `<p style="color:#a60;"><small>⚠ No repeated observations: R values are undefined (data already merged to a unique set).</small></p>` : '';

  container.innerHTML = `
    <p><small>Merged under Laue class ${info.laue}${info.source === 'cif' ? ' (from the CIF space group)' : info.source === 'merging' ? ' (confirmed by merging)' : ' (assumed — highest-symmetry candidate for this cell, not confirmed)'}; ${info.lattice}.
      Shells have equal width in d*³ (the same number of possible unique reflections each), from low to high resolution.
      Where the mean intensity approaches zero the R values grow without limit — use ⟨I/σ⟩ there.</small></p>
    <table class="shellTable">
      <tr><th>d range</th><th>Unique</th><th>Complete</th><th>Redundancy</th><th>⟨I/σ⟩</th><th>I &gt; 2σ</th>
          <th>R<sub>int</sub></th><th>R<sub>meas</sub></th><th>R<sub>pim</sub></th></tr>
      ${rows}${total}
    </table>
    ${redNote}${cutNote}
    <details class="shellLegend">
      <summary>What do these columns mean?</summary>
      <ul>
        <li><strong>d range</strong> — resolution shell in Å (low to high resolution).</li>
        <li><strong>Unique</strong> / <strong>Redundancy</strong> — symmetry-unique reflections measured in the shell, and how often each was measured on average.</li>
        <li><strong>Complete</strong> — measured unique reflections divided by the unique reflections the cell, lattice centering and Laue class allow in that shell (hover for the counts). Below 95 % is shown in orange: gaps in the data (beamstop, detector gaps, short scan). The bottom row is the completeness up to the highest measured resolution. Reflections that the lattice centering forbids, or that the absence scan found to be systematically absent (screw axes, glides), count neither as possible nor as measured — so files already cleaned of absences (deposited CIF/.hkl data) are not penalised. If a systematic absence was NOT recognised by the scan, the value is too low.</li>
        <li><strong>⟨I/σ⟩</strong> — mean signal-to-noise of the merged reflections (intensity divided by its standard uncertainty σ). Higher is better; a common rule of thumb puts the resolution limit where it falls below about 2.</li>
        <li><strong>I &gt; 2σ</strong> — share of merged reflections whose intensity is more than twice its uncertainty ("observed" reflections).</li>
        <li><strong>R<sub>int</sub></strong> — spread of equivalent measurements around their mean, Σ|I − ⟨I⟩| / ΣI. It grows with redundancy even when the data get better.</li>
        <li><strong>R<sub>meas</sub></strong> — R<sub>int</sub> corrected for redundancy (factor √(n/(n−1))): the quality of a single measurement, comparable between data sets.</li>
        <li><strong>R<sub>pim</sub></strong> — R<sub>int</sub> corrected with √(1/(n−1)): the precision of the merged (averaged) intensities that go into refinement; always smaller than R<sub>meas</sub>.</li>
      </ul>
      <p><small>R values are only defined for reflections measured at least twice, and they diverge where the mean intensity approaches zero — rely on ⟨I/σ⟩ there.</small></p>
    </details>`;
}
