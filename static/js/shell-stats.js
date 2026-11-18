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

function accumulateShell(chunk) {
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
  const sMin = chunk[0].s, sMax = chunk[chunk.length - 1].s;
  return {
    dMax: 1 / (2 * sMin), dMin: 1 / (2 * sMax),
    nUnique: chunk.length, nTotal, redundancy: nTotal / chunk.length,
    meanISigma: nISig ? sumISig / nISig : null,
    fractionGt2: nISig ? nGt2 / nISig : null,
    rInt: den > 0 ? numInt / den : null,
    rMeas: den > 0 ? numMeas / den : null,
    rPim: den > 0 ? numPim / den : null,
  };
}

// `reflections` raw (possibly redundant); `groupMatrices` from generateGroup().
// Shells hold (roughly) equal numbers of unique reflections.
function computeShellStatistics(reflections, cell, groupMatrices, targetShells = 10) {
  const groups = groupReflectionsByCanonicalHKL(reflections, groupMatrices);
  const recip = reciprocalCellFromDirect(cell);
  const uniques = [];
  for (const g of groups.values()) {
    const s = sinThetaOverLambda(g.h, g.k, g.l, recip);
    if (s > 0) uniques.push({ s, members: g.members });
  }
  if (uniques.length < 40) {
    return { error: `Only ${uniques.length} unique reflections — too few for resolution shells.` };
  }
  uniques.sort((a, b) => a.s - b.s);
  const nShells = Math.max(1, Math.min(targetShells, Math.floor(uniques.length / 20)));
  const per = Math.ceil(uniques.length / nShells);
  const shells = [];
  for (let i = 0; i < uniques.length; i += per) shells.push(accumulateShell(uniques.slice(i, i + per)));
  const overall = accumulateShell(uniques);
  // first shell (from low to high resolution) where <I/σ> falls below 2: a
  // common rule of thumb puts the resolution limit at or just before it
  const cutoffIndex = shells.findIndex(sh => sh.meanISigma !== null && sh.meanISigma < 2);
  return { shells, overall, cutoffIndex };
}

function renderShellStatistics(container, stats, info) {
  if (stats.error) { container.innerHTML = `<p><small>${stats.error}</small></p>`; return; }
  const pct = x => (x === null ? '—' : (x * 100).toFixed(1) + ' %');
  const num = (x, d = 1) => (x === null ? '—' : x.toFixed(d));
  const row = (label, sh, cls) => `
    <tr${cls ? ` class="${cls}"` : ''}>
      <td>${label}</td><td>${sh.nUnique}</td><td>${sh.redundancy.toFixed(2)}</td>
      <td>${num(sh.meanISigma)}</td><td>${pct(sh.fractionGt2)}</td>
      <td>${pct(sh.rInt)}</td><td>${pct(sh.rMeas)}</td><td>${pct(sh.rPim)}</td>
    </tr>`;
  const rows = stats.shells.map((sh, i) =>
    row(`${sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)} Å`, sh, i === stats.cutoffIndex ? 'shellCut' : '')).join('');
  const total = row(`all (${stats.overall.dMax.toFixed(2)} – ${stats.overall.dMin.toFixed(2)} Å)`, stats.overall, 'shellTotal');

  let cutNote;
  if (stats.cutoffIndex >= 0) {
    const sh = stats.shells[stats.cutoffIndex];
    cutNote = `<p><small>⟨I/σ⟩ drops below 2 in the highlighted shell (${sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)} Å). A common rule of thumb puts the resolution limit at or just before such a shell.
      This is shown for orientation only — nothing is cut or written to the .ins/.hkl.</small></p>`;
  } else {
    cutNote = `<p><small>⟨I/σ⟩ stays above 2 in every shell — no obvious resolution cut-off in this data set. Nothing is cut or written to the .ins/.hkl.</small></p>`;
  }
  const redNote = stats.overall.rInt === null
    ? `<p style="color:#a60;"><small>⚠ No repeated observations: R values are undefined (data already merged to a unique set).</small></p>` : '';

  container.innerHTML = `
    <p><small>Merged under Laue class ${info.laue}${info.assumed ? ' (assumed, not confirmed by merging)' : ''}; ${info.lattice}.
      Shells hold about the same number of unique reflections each, from low to high resolution.
      Where the mean intensity approaches zero the R values grow without limit — use ⟨I/σ⟩ there.</small></p>
    <table class="shellTable">
      <tr><th>d range</th><th>Unique</th><th>Redundancy</th><th>⟨I/σ⟩</th><th>I &gt; 2σ</th>
          <th>R<sub>int</sub></th><th>R<sub>meas</sub></th><th>R<sub>pim</sub></th></tr>
      ${rows}${total}
    </table>
    ${redNote}${cutNote}
    <details class="shellLegend">
      <summary>What do these columns mean?</summary>
      <ul>
        <li><strong>d range</strong> — resolution shell in Å (low to high resolution).</li>
        <li><strong>Unique</strong> / <strong>Redundancy</strong> — symmetry-unique reflections in the shell, and how often each was measured on average.</li>
        <li><strong>⟨I/σ⟩</strong> — mean signal-to-noise of the merged reflections (intensity divided by its standard uncertainty σ). Higher is better; a common rule of thumb puts the resolution limit where it falls below about 2.</li>
        <li><strong>I &gt; 2σ</strong> — share of merged reflections whose intensity is more than twice its uncertainty ("observed" reflections).</li>
        <li><strong>R<sub>int</sub></strong> — spread of equivalent measurements around their mean, Σ|I − ⟨I⟩| / ΣI. It grows with redundancy even when the data get better.</li>
        <li><strong>R<sub>meas</sub></strong> — R<sub>int</sub> corrected for redundancy (factor √(n/(n−1))): the quality of a single measurement, comparable between data sets.</li>
        <li><strong>R<sub>pim</sub></strong> — R<sub>int</sub> corrected with √(1/(n−1)): the precision of the merged (averaged) intensities that go into refinement; always smaller than R<sub>meas</sub>.</li>
      </ul>
      <p><small>R values are only defined for reflections measured at least twice, and they diverge where the mean intensity approaches zero — rely on ⟨I/σ⟩ there.</small></p>
    </details>`;
}
