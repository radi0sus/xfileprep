    const dropZone = document.getElementById('dropZone');
    const input = document.getElementById('hklInput');
    const sourceInfoEl = document.getElementById('sourceInfo');
    const cellFields = {
      a: document.getElementById('cellA'),
      b: document.getElementById('cellB'),
      c: document.getElementById('cellC'),
      alpha: document.getElementById('cellAlpha'),
      beta: document.getElementById('cellBeta'),
      gamma: document.getElementById('cellGamma'),
    };
    const cellNoteEl = document.getElementById('cellNote');
    const cellExtraEl = document.getElementById('cellExtra');
    const cellEsdEl = document.getElementById('cellEsd');
    const cellMetricResultEl = document.getElementById('cellMetricResult');
    const formulaInput = document.getElementById('formulaInput');
    const applyFormulaBtn = document.getElementById('applyFormulaBtn');
    const formulaUpdateNoteEl = document.getElementById('formulaUpdateNote');
    const laueBtn = document.getElementById('determineLaueBtn');
    const laueResultEl = document.getElementById('laueResult');
    let lastReflections = null;
    let determinedLaueClass = null;
    let lastRawText = null;
    let lastCellExtras = null;
    let determinedLaueClassSource = null; // 'merging' | null
    let lastMatchResults = null;
    let selectedCandidateIndex = 0;
    let selectionSource = 'default'; // 'default' | 'manual' | 'auto-e-stat'
    let lastInsComposition = null; // exact whole-unit-cell composition, from a loaded .ins's SFAC/UNIT
    let lastCellEsd = null; // {a,b,c,alpha,beta,gamma} esds from whichever source (CIF/.ins/.p4p) provided them
    // Tracks whether the person has typed into the formula field themselves
    // since the last file load. When true, a manually-typed formula (× Z)
    // is used for the .ins SFAC/UNIT instead of an exact whole-cell
    // composition from a previously loaded .ins/.p4p — otherwise editing
    // the formula field had no visible effect whenever a composition was
    // already known, since composition always took priority (see Bug 8).
    let formulaManuallyEdited = false;
    let lastLoadedFileBaseName = null; // filename (no extension) of whichever file was loaded most recently, used as a TITL prefix

    // Renders screw-axis notation (21, 31, 32, 41, 42, 43, 61-65) with a
    // proper Unicode subscript for display, e.g. "P 1 21 1" -> "P 1 2₁ 1".
    // Only for on-screen text — the SHELX .ins TITL line deliberately keeps
    // plain ASCII digits for maximum compatibility with older software.
    const SUBSCRIPT_DIGITS = { '0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄', '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉' };
    const SCREW_AXIS_TOKENS = /\b(21|31|32|41|42|43|61|62|63|64|65)\b/g;
    function formatHM(hm) {
      return hm.replace(SCREW_AXIS_TOKENS, m => m[0] + m.slice(1).split('').map(d => SUBSCRIPT_DIGITS[d]).join(''));
    }

    const absencesBtn = document.getElementById('scanAbsencesBtn');
    const absencesResultEl = document.getElementById('absencesResult');

    const VERDICT_STYLE = {
      holds: 'color:#080', 
      violated: 'color:#aaa',
      inconclusive: 'color:#b90',
    };
    function verdictCell(r) {
      const style = VERDICT_STYLE[r.verdict] || '';
      const title = r.reason ? ` title="${r.reason.replace(/"/g, '&quot;')}"` : '';
      return `<td style="${style}"${title}>${r.verdict} <small>(${r.violationPercent.toFixed(1)}%)</small></td>`;
    }

    function runAbsencesScan() {
      return new Promise((resolve) => {
        if (!lastReflections) { resolve(false); return; }
        absencesResultEl.innerHTML = `<p><em>Scanning…</em></p>`;

        setTimeout(() => {
        const scan = scanAllAbsences(lastReflections);

        const centeringRows = Object.entries(scan.centering)
          .map(([name, r]) => `<tr><td>${name}</td>${verdictCell(r)}</tr>`).join('');

        const axialRows = Object.entries(scan.axial).map(([axis, results]) => {
          const cells = results.map(r => verdictCell(r)).join('');
          return `<tr><td>${axis}</td><td>${results[0].reflectionCount}</td>${cells}</tr>`;
        }).join('');

        const zonalRows = Object.entries(scan.zonal).flatMap(([zone, results]) =>
          results.map((r, i) => `<tr><td>${i === 0 ? zone : ''}</td><td>${r.name}</td>${verdictCell(r)}</tr>`)
        ).join('');

        absencesResultEl.innerHTML = `
          <h2>Systematic absences</h2>

          <h3>Lattice centering (whole dataset)</h3>
          <table><tr><th>Type</th><th>Verdict</th></tr>${centeringRows}</table>

          <h3>Axial conditions (screw axes)</h3>
          <table>
            <tr><th>Axis</th><th>n refl.</th><th>period 2</th><th>period 3</th><th>period 4</th><th>period 6</th></tr>
            ${axialRows}
          </table>

          <h3>Zonal conditions (glide planes)</h3>
          <table><tr><th>Zone</th><th>Condition</th><th>Verdict</th></tr>${zonalRows}</table>

          <p><small>Percentage = share of "forbidden" reflections that are still significant (I/σ ≥ 3) — 0% is the clean case.
          A zonal condition consistent with the lattice centering alone (e.g. an apparent glide in a zone that intersects the
          centering condition) doesn't imply an independent glide plane — check the centering result first.
          <strong>"Inconclusive"</strong> means the test itself isn't trustworthy either way — hover a cell for the specific reason.
          It happens for one of two causes: too few reflections in the allowed (reference) class to say anything, or the allowed
          class itself is mostly weak too (common for a condition that only really shows up at higher resolution) — in both cases
          the percentage next to it isn't meaningful evidence, whether it's 0% or high.</small></p>`;

        document.getElementById('absenceVizSection').style.display = 'block';
        renderAbsenceVisuals(
          document.getElementById('zonePlotsContainer'),
          document.getElementById('axialStripsContainer'),
          lastReflections
        );
        resolve(true);
        }, 10);
      });
    }
    absencesBtn.addEventListener('click', () => { runAbsencesScan(); });


    // Reads the current cell fields and (re-)shows the metric classification.
    // Runs both after auto-fill from CIF and on any manual edit of the fields.
    function updateCellMetricClassification() {
      const values = {
        a: parseFloat(cellFields.a.value),
        b: parseFloat(cellFields.b.value),
        c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value),
        beta: parseFloat(cellFields.beta.value),
        gamma: parseFloat(cellFields.gamma.value),
      };
      if (!Object.values(values).every(Number.isFinite)) {
        cellMetricResultEl.innerHTML = '';
        return;
      }
      const classification = classifyCellMetric(values);
      cellMetricResultEl.innerHTML = `
        <h2>Cell metric symmetry</h2>
        <p><strong>Crystal system:</strong> ${classification.system}</p>
        <p><strong>Laue class candidates:</strong> ${classification.laueClasses.join(', ')}</p>
        <p><small>Based on cell metric only — the actual Laue class still needs to be confirmed by comparing merging statistics for each candidate (next step).</small></p>`;
    }
    Object.values(cellFields).forEach(field => field.addEventListener('input', updateCellMetricClassification));
    // For raw HKL (no cell in the file): once all six values are filled in
    // by hand and the field loses focus, run the same auto-pipeline that
    // fires automatically after a CIF load — 'change' (not 'input') so it
    // doesn't refire on every keystroke.
    Object.values(cellFields).forEach(field => field.addEventListener('change', () => autoRunPipeline()));
    // Same idea for the formula field: E² statistics needs a composition,
    // which for raw HKL / no-.ins data usually only shows up once the
    // person types a formula by hand — pick that up automatically too,
    // rather than requiring a manual "Compute E² statistics" click. Also
    // wired to an explicit "Apply formula" button/click, since the results
    // (Wilson plot, .ins) live on other tabs and a silent background
    // recompute is easy to miss — this button gives a deliberate action
    // plus visible confirmation right next to the field.
    function applyFormulaChange() {
      let did = [];
      if (lastReflections && lastMatchResults) { runEStatistics(true); did.push('E² statistics'); }
      if (lastMatchResults) { generateIns(); did.push('the SHELX .ins'); } // formula/Z feeds directly into the .ins too, regardless of E-stat outcome
      formulaUpdateNoteEl.innerHTML = did.length > 0
        ? `<p style="color:#080"><small>✓ Recalculated ${did.join(' and ')} with the current formula — check the corresponding tab.</small></p>`
        : `<p style="color:#b90"><small>No space-group match yet to recalculate — load reflection data first.</small></p>`;
    }
    formulaInput.addEventListener('change', applyFormulaChange);
    applyFormulaBtn.addEventListener('click', applyFormulaChange);
    formulaInput.addEventListener('input', () => { formulaManuallyEdited = true; });
    const summaryEl = document.getElementById('summary');
    const warningsEl = document.getElementById('warnings');

    const SOURCE_LABELS = {
      raw_hkl: 'Raw HKL file',
      embedded_shelx_hkl_file: 'Embedded HKL from CIF (_shelx_hkl_file)',
      refln_loop: 'CIF _refln loop',
    };

    // Shows the six cell esds (standard uncertainties) as given by whichever
    // source file provided them — CIF parenthesis notation, a SHELX .ins
    // ZERR line, or a .p4p CELLSD line. Clears the display when nothing
    // (or only an all-zero placeholder) was available.
    function renderCellEsd(esd) {
      lastCellEsd = esd;
      if (!esd) { cellEsdEl.innerHTML = ''; return; }
      const fmt = (v) => (v === null || v === undefined || !Number.isFinite(v)) ? '—' : v;
      cellEsdEl.innerHTML = `<p><small>Standard uncertainties (as given by the source file):
        a = ${fmt(esd.a)} · b = ${fmt(esd.b)} · c = ${fmt(esd.c)} Å ·
        α = ${fmt(esd.alpha)} · β = ${fmt(esd.beta)} · γ = ${fmt(esd.gamma)} °</small></p>`;
    }

    // Shared by both CIF-loading paths (a full CIF with reflections, and a
    // CIF with no reflection data at all — see handleFile) so cell/formula
    // get filled in identically either way.
    function applyCellFromCIF(cell, text) {
      lastCellExtras = cell;
      cellFields.a.value = cell.a;
      cellFields.b.value = cell.b;
      cellFields.c.value = cell.c;
      cellFields.alpha.value = cell.alpha;
      cellFields.beta.value = cell.beta;
      cellFields.gamma.value = cell.gamma;
      cellNoteEl.innerHTML = `<p><em>Filled in from CIF — edit if needed.</em></p>`;

      const extras = [];
      if (cell.z !== null) extras.push(`Z = ${cell.z}`);
      if (cell.wavelength !== null) extras.push(`λ = ${cell.wavelength} Å`);
      if (cell.reportedSpaceGroup !== null) extras.push(`reported space group: ${cell.reportedSpaceGroup}`);
      cellExtraEl.innerHTML = extras.length > 0 ? `<p>${extras.join(' · ')}</p>` : '';

      const formula = parseChemicalFormula(text);
      formulaInput.value = formula ? formula.map(f => f.element + f.count).join(' ') : '';
      formulaManuallyEdited = false; // fresh load
      renderCellEsd(cell.esd);
    }

    // A .ins carries cell + wavelength + Z + exact whole-cell composition
    // (SFAC/UNIT) but no reflections — fills in what it can, leaves any
    // already-loaded reflection data untouched.
    function handleInsFile(text) {
      const insData = parseInsFile(text);
      if (!insData) {
        sourceInfoEl.innerHTML = `<p style="color:#b00">Couldn't parse this .ins — need CELL, SFAC, and UNIT lines.</p>`;
        return;
      }

      cellFields.a.value = insData.cell.a;
      cellFields.b.value = insData.cell.b;
      cellFields.c.value = insData.cell.c;
      cellFields.alpha.value = insData.cell.alpha;
      cellFields.beta.value = insData.cell.beta;
      cellFields.gamma.value = insData.cell.gamma;
      cellNoteEl.innerHTML = `<p><em>Cell filled in from .ins — edit if needed.</em></p>`;

      lastInsComposition = insData.composition;
      lastCellExtras = { ...(lastCellExtras || {}), wavelength: insData.wavelength, z: insData.z, esd: insData.esd };
      formulaManuallyEdited = false; // fresh load — the exact composition above should win until the person edits it themselves

      // The formula field is documented as "per formula unit" — a loaded
      // .ins only gives the whole-cell composition (SFAC/UNIT) directly,
      // so divide by Z to show the same thing in the same convention as a
      // CIF-derived formula. Whole-cell counts stay in lastInsComposition
      // regardless and are what .ins generation actually uses (exact,
      // untouched by this division) — this is purely for display/parity
      // with the CIF path, plus a manual edit starting point if needed.
      if (insData.z && insData.z > 0) {
        formulaInput.value = insData.composition
          .map(c => {
            const perUnit = c.count / insData.z;
            const rounded = Math.round(perUnit);
            const isExact = Math.abs(perUnit - rounded) < 1e-6;
            return c.element + (isExact ? rounded : perUnit.toFixed(2));
          })
          .join(' ');
      }

      const extras = [];
      if (insData.z !== null) extras.push(`Z = ${insData.z}`);
      if (insData.wavelength !== null) extras.push(`λ = ${insData.wavelength} Å`);
      extras.push(`composition (whole cell, from SFAC/UNIT): ${insData.composition.map(c => c.element + c.count).join(' ')}`);
      cellExtraEl.innerHTML = `<p>${extras.join(' · ')}</p>`;
      renderCellEsd(insData.esd);

      updateCellMetricClassification();
      sourceInfoEl.innerHTML = `<p><strong>Source:</strong> SHELX .ins (cell + composition only — drop the matching .hkl for reflection data too, if you haven't already).</p>`;

      // If an .hkl was already dropped first, this only just now completed
      // the cell — pick up the auto-run from here rather than requiring
      // the person to re-drop the .hkl or click through manually.
      autoRunPipeline();
    }

    // A .p4p (CrysAlisPro/Bruker orientation-matrix export) carries cell +
    // esds + wavelength + formula but, like a .ins, no reflections — same
    // "fill in what it can" contract as handleInsFile.
    function handleP4PFile(text) {
      const p4pData = parseP4PFile(text);
      if (!p4pData) {
        sourceInfoEl.innerHTML = `<p style="color:#b00">Couldn't parse this .p4p — need a CELL line with six cell constants.</p>`;
        return;
      }

      cellFields.a.value = p4pData.cell.a;
      cellFields.b.value = p4pData.cell.b;
      cellFields.c.value = p4pData.cell.c;
      cellFields.alpha.value = p4pData.cell.alpha;
      cellFields.beta.value = p4pData.cell.beta;
      cellFields.gamma.value = p4pData.cell.gamma;
      cellNoteEl.innerHTML = `<p><em>Cell filled in from .p4p — edit if needed.</em></p>`;

      // .p4p gives a formula (per formula unit, like a CIF), not an exact
      // whole-cell SFAC/UNIT count like a .ins — don't let a stale exact
      // composition from an earlier .ins linger and silently take priority.
      lastInsComposition = null;
      formulaManuallyEdited = false; // fresh load
      lastCellExtras = { ...(lastCellExtras || {}), wavelength: p4pData.wavelength, z: lastCellExtras?.z ?? null };

      const extras = [];
      if (p4pData.wavelength !== null) extras.push(`λ = ${p4pData.wavelength} Å`);
      if (p4pData.volume !== null) extras.push(`V = ${p4pData.volume} Å³`);
      cellExtraEl.innerHTML = extras.length > 0 ? `<p>${extras.join(' · ')}</p>` : '';
      renderCellEsd(p4pData.esd);

      if (p4pData.formula) {
        formulaInput.value = p4pData.formula.map(f => f.element + f.count).join(' ');
      }

      updateCellMetricClassification();
      sourceInfoEl.innerHTML = `<p><strong>Source:</strong> .p4p (cell${p4pData.wavelength !== null ? ' + wavelength' : ''}${p4pData.formula ? ' + formula' : ''} —
        drop the matching .hkl for reflection data too, if you haven't already).</p>`;

      autoRunPipeline();
    }

    async function handleFile(file) {
      if (!file) return;
      const text = await file.text();
      if (file.name) lastLoadedFileBaseName = file.name.replace(/\.[^./\\]+$/, ''); // strip extension (.hkl/.cif/.ins/.p4p/...)

      // A .p4p (Bruker/CrysAlisPro orientation-matrix export) carries cell +
      // esds + wavelength + formula but no reflections at all — same
      // reasoning as the .ins branch below, handle it separately rather
      // than forcing it through the reflection-parsing path.
      if (looksLikeP4P(text)) {
        handleP4PFile(text);
        return;
      }

      // A .ins carries cell + exact whole-cell composition (SFAC/UNIT), but
      // no reflections at all — handle it separately rather than forcing
      // it through the reflection-parsing path, which would just report
      // "no reflection data found".
      if (looksLikeSHELXIns(text)) {
        handleInsFile(text);
        return;
      }

      const result = parseReflectionInput(text);

      if (result.reflections.length === 0) {
        // Sniffed as CIF, but no _refln_index_h loop or embedded
        // _shelx_hkl_file block — many deposited/reduced CIFs only carry
        // cell + formula, with reflections shipped separately. Rather than
        // just erroring out, grab whatever the CIF does have (same as the
        // .ins/.p4p "cell-only" files) and wait for a matching HKL to be
        // dropped afterwards — in either order, since a raw HKL with no
        // cell info already keeps whatever cell was established earlier.
        // NB: this branch has to trigger regardless of result.skipped —
        // every non-blank CIF line (cell/formula tags etc.) fails the HKL
        // fallback parser and lands in `skipped`, so requiring an empty
        // skipped list here would never fire for exactly the CIFs this is
        // meant to catch, and the file would fall through into the normal
        // "we have reflections" branch below with an empty array instead
        // (0 reflections, but not rejected) — silently running the whole
        // pipeline on nothing and rendering "undefined" ranges.
        if (looksLikeCIF(text)) {
          const cellOnly = extractCellFromCIF(text);
          if (cellOnly) {
            lastRawText = text;
            applyCellFromCIF(cellOnly, text);
            sourceInfoEl.innerHTML = `<p style="color:#b90">This CIF has no reflection data (no <code>_refln_index_h</code> loop or embedded HKL) —
              cell and formula were loaded from it below. Drop the matching .hkl file to continue.</p>`;
            summaryEl.innerHTML = '';
            warningsEl.innerHTML = '';
            updateCellMetricClassification();
            autoRunPipeline(); // in case reflections were already loaded earlier and were just waiting on this cell
            return;
          }
        }
        sourceInfoEl.innerHTML = `<p style="color:#b00">${result.skipped.length > 0
          ? `No valid reflections found — ${result.skipped.length} line(s) in this file couldn't be parsed as HKL data.`
          : 'No reflection data found in this file.'}</p>`;
        summaryEl.innerHTML = '';
        warningsEl.innerHTML = '';
        return;
      }

      const { reflections, skipped, source } = result;
      lastReflections = reflections;
      determinedLaueClass = null;
      determinedLaueClassSource = null;
      laueBtn.style.display = 'inline-block';
      laueResultEl.innerHTML = '';
      absencesBtn.style.display = 'inline-block';
      absencesResultEl.innerHTML = '';
      spaceGroupBtn.style.display = 'inline-block';
      spaceGroupResultEl.innerHTML = '';
      eStatBtn.style.display = 'inline-block';
      eStatResultEl.innerHTML = '';
      const stats = summarizeHKL(reflections);

      sourceInfoEl.innerHTML = `<p><strong>Source:</strong> ${SOURCE_LABELS[source] || source}</p>`;

      // For CIF input, cell constants are already present in the file —
      // pre-fill them instead of asking again. For raw HKL, there's no
      // cell info at all, so fields stay empty for manual entry — UNLESS
      // an .ins/.p4p was already dropped first and already established the
      // cell (+ wavelength/Z/composition/esds/formula); in that case this
      // reflection-only file must not wipe any of that out.
      const cell = looksLikeCIF(text) ? extractCellFromCIF(text) : null;
      lastRawText = text;
      if (cell) {
        // A full CIF (this one has reflections, unlike the no-HKL branch
        // above) is always authoritative for its own cell/formula — even
        // if a cell-only CIF or an .ins/.p4p was dropped earlier, this
        // replaces it outright, same as dropping any fresh CIF always has.
        applyCellFromCIF(cell, text);
      } else if (lastCellExtras || lastInsComposition) {
        // No cell info in THIS file, but an .ins/.p4p already established
        // one earlier in this session (in either drop order) — keep it,
        // and the cell fields, exactly as they were. This file only
        // contributes reflections.
        cellNoteEl.innerHTML = `<p><em>No cell info in this file — keeping the cell already filled in from an earlier .ins/.p4p.</em></p>`;
      } else {
        cellNoteEl.innerHTML = looksLikeCIF(text)
          ? `<p style="color:#b00"><em>No cell constants found in this CIF — please enter manually.</em></p>`
          : `<p><em>Raw HKL has no cell info — please enter manually.</em></p>`;
        cellExtraEl.innerHTML = '';
        formulaInput.value = '';
        formulaManuallyEdited = false; // fresh load
        renderCellEsd(null);
      }
      updateCellMetricClassification();

      summaryEl.innerHTML = `
        <h2>Summary</h2>
        <table>
          <tr><td>Reflection count</td><td>${stats.count}</td></tr>
          <tr><td>h range</td><td>${stats.hRange?.join(' … ')}</td></tr>
          <tr><td>k range</td><td>${stats.kRange?.join(' … ')}</td></tr>
          <tr><td>l range</td><td>${stats.lRange?.join(' … ')}</td></tr>
          <tr><td>F² range</td><td>${stats.f2Range?.map(v => v.toFixed(2)).join(' … ')}</td></tr>
          <tr><td>Mean F²</td><td>${stats.meanF2?.toFixed(2)}</td></tr>
        </table>`;

      if (skipped.length > 0) {
        warningsEl.innerHTML = `
          <h2 style="color:#b00">⚠ ${skipped.length} line(s) could not be parsed</h2>
          <pre>${skipped.slice(0, 5).map(s => `Line ${s.lineNumber ?? '?'}: "${s.content}"`).join('\n')}</pre>`;
      } else {
        warningsEl.innerHTML = '';
      }

      // Cell is already known (e.g. CIF autofill) — no need to make the
      // person click through Laue class → absences → space group match
      // one button at a time; run the whole chain straight to the
      // candidate table. If cell constants are still missing (raw HKL),
      // this quietly does nothing until they're filled in by hand — see
      // the cellFields 'change' listener below, which picks up from there.
      autoRunPipeline();
    }

    // Runs determineLaue → scanAbsences → matchSpaceGroup back to back and
    // lands the person directly at the space-group candidate table,
    // instead of requiring three separate button clicks. Safe to call
    // repeatedly (e.g. after editing the cell by hand) — each step is a
    // no-op if its own prerequisites aren't met yet.
    async function autoRunPipeline() {
      if (!lastReflections) return;
      const values = {
        a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
      };
      if (!Object.values(values).every(Number.isFinite)) return;

      await runLaueDetermination();
      await runAbsencesScan();
      const matched = await runSpaceGroupMatch();
      if (matched) {
        // Wilson-plot/E² statistics only if a composition is already known
        // (exact .ins SFAC/UNIT, or a filled-in formula + Z) — otherwise
        // leave it for a manual click once the formula is entered.
        await runEStatistics(true);
      }
      if (matched && typeof window.switchTab === 'function') {
        window.switchTab('tab-spacegroup');
      }
    }

    function runLaueDetermination() {
      return new Promise((resolve) => {
        if (!lastReflections) { resolve(false); return; }

        const values = {
          a: parseFloat(cellFields.a.value),
          b: parseFloat(cellFields.b.value),
          c: parseFloat(cellFields.c.value),
          alpha: parseFloat(cellFields.alpha.value),
          beta: parseFloat(cellFields.beta.value),
          gamma: parseFloat(cellFields.gamma.value),
        };
        if (!Object.values(values).every(Number.isFinite)) {
          laueResultEl.innerHTML = `<p style="color:#b00">Please fill in all six cell constants first.</p>`;
          resolve(false);
          return;
        }

        const classification = classifyCellMetric(values);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';

        laueResultEl.innerHTML = `<p><em>Merging reflections under ${classification.laueClasses.length} candidate(s) — this may take a moment for large datasets…</em></p>`;

        // deferred so the "merging…" message actually paints before the (synchronous) work runs
        setTimeout(() => {
        const { chosen, jumpDetectedAt, evaluations, redundancyInformative } =
          determineLaueClass(lastReflections, classification.laueClasses, { monoclinicAxis });

        // With no redundant (duplicate) observations under ANY candidate,
        // Rint is undefined everywhere and this method has no signal at
        // all — "chosen" would otherwise silently default to the lowest
        // candidate, which looks like a confident answer but isn't one.
        // Don't treat that as determined; leave it unset so the
        // space-group matching step keeps warning the user honestly
        // instead of quietly trusting it.
        determinedLaueClass = redundancyInformative ? chosen : null;
        determinedLaueClassSource = redundancyInformative ? 'merging' : null;

        const rows = evaluations.map(r => `
          <tr ${r.laueClass === chosen ? 'style="font-weight:bold;background:#f3fbf6;"' : ''}>
            <td>${r.laueClass}${r.laueClass === chosen ? ' ✓' : ''}</td>
            <td>${r.groupOrder}</td>
            <td>${r.uniqueReflections}</td>
            <td>${r.redundancy.toFixed(2)}</td>
            <td>${r.rInt !== null ? (r.rInt * 100).toFixed(2) + ' %' : '—'}</td>
          </tr>`).join('');

        const jumpNote = jumpDetectedAt
          ? `<p><small>Stopped before <strong>${jumpDetectedAt}</strong> — R<sub>int</sub> jumped there by more than the 50% tolerance, so that candidate's extra symmetry doesn't seem to be real.</small></p>`
          : '';

        const verdictBlock = redundancyInformative
          ? `<p><strong>Auto-picked: ${chosen}</strong> — highest symmetry candidate with no R<sub>int</sub> jump from the baseline.</p>
             ${jumpNote}
             <p><small>This is a heuristic (50% relative R<sub>int</sub> increase tolerated) — cross-check against the table above, especially in borderline cases.</small></p>`
          : `<p style="color:#b00"><strong>⚠ No redundant (duplicate) reflections found under ANY candidate symmetry — this method has no signal to pick between them.</strong>
             This happens when the reflection list is already merged to symmetry-unique data (e.g. a deposited CIF <code>_refln</code> loop) rather than
             raw, multiply-measured data straight off the diffractometer. Merging statistics fundamentally need redundancy to work with; there's nothing
             here to compare. The table above is shown for transparency but its "auto-pick" would just be the lowest-symmetry candidate by default (<strong>${chosen}</strong>) —
             don't trust that. If you have the original unmerged HKL (with duplicate/batch observations), use that instead; otherwise the systematic-absence
             pattern (zone plots below, once you've scanned) is your best remaining evidence, read by eye.</p>`;

        laueResultEl.innerHTML = `
          <h2>Laue class candidates (merging statistics)</h2>
          <table>
            <tr><th>Laue class</th><th>Order</th><th>Unique refl.</th><th>Redundancy</th><th>R<sub>int</sub></th></tr>
            ${rows}
          </table>
          ${verdictBlock}`;
        resolve(true);
        }, 10);
      });
    }
    laueBtn.addEventListener('click', () => { runLaueDetermination(); });

    const spaceGroupBtn = document.getElementById('matchSpaceGroupBtn');
    const spaceGroupResultEl = document.getElementById('spaceGroupResult');

    const insBtn = document.getElementById('generateInsBtn');
    const insResultEl = document.getElementById('insResult');
    const insActionsEl = document.getElementById('insActions');
    const insCopyBtn = document.getElementById('insCopyBtn');
    const insSaveBtn = document.getElementById('insSaveBtn');
    const insCopyNoteEl = document.getElementById('insCopyNote');
    let lastGeneratedIns = null; // raw .ins text of whatever's currently shown, for the Copy/Save buttons

    // Builds the SHELX .ins for whichever candidate is currently selected
    // (selectedCandidateIndex) and renders it. Called automatically —
    // right after a space-group match, whenever the radio selection
    // changes, whenever E² statistics resolves a centro/acentric tie, and
    // whenever cell/formula fields change — as well as from the manual
    // button, so the .ins is always kept in sync with the current pick
    // rather than requiring a click to (re)generate it.
    function generateIns() {
      if (!lastMatchResults || lastMatchResults.length === 0) {
        insResultEl.innerHTML = '';
        insActionsEl.style.display = 'none';
        lastGeneratedIns = null;
        return;
      }

      const top = lastMatchResults[selectedCandidateIndex].entry;

      // If the SELECTED candidate is tied with another that disagrees on
      // centrosymmetry, absences alone couldn't break that tie on their
      // own — note that plainly, but phrase it according to how this
      // particular selection actually came about (see selectionSource).
      const topScore = lastMatchResults[selectedCandidateIndex].score;
      const topIsCentro = hallSymbolIsCentrosymmetric(top.hall);
      const tiedOpposite = lastMatchResults.find(m =>
        m.entry !== top && m.score === topScore && hallSymbolIsCentrosymmetric(m.entry.hall) !== topIsCentro);

      const values = {
        a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
      };
      const formula = parseFormulaString(formulaInput.value);
      const wavelength = lastCellExtras?.wavelength;
      let z = lastCellExtras?.z;
      let zWasEstimated = false;

      if (!z && formula.length > 0) {
        const volume = cellVolume(values.a, values.b, values.c, values.alpha, values.beta, values.gamma);
        z = estimateZ(volume, formula);
        zWasEstimated = z !== null;
      }

      const spaceGroupTag = top.hm.replace(/\s+/g, '');
      const ins = generateInsFile({
        title: lastLoadedFileBaseName ? `${lastLoadedFileBaseName} in ${spaceGroupTag}` : spaceGroupTag,
        cell: values,
        wavelength,
        hallSymbol: top.hall,
        z,
        formula: formula.length > 0 ? formula : null,
        composition: formulaManuallyEdited ? null : lastInsComposition, // exact whole-cell counts, if a .ins/.res was loaded — takes priority over formula×Z UNLESS the person has since typed their own formula
        esd: lastCellEsd, // cell esds, if the source file provided them (CIF parens / .ins ZERR / .p4p CELLSD)
      });

      const hasComposition = (lastInsComposition && lastInsComposition.length > 0) || formula.length > 0;
      const missingParts = [!wavelength && 'wavelength', !hasComposition && 'chemical formula'].filter(Boolean);
      const missingNote = missingParts.length > 0
        ? `<p style="color:#b00"><small>Missing: ${missingParts.join(', ')} — placeholder(s) inserted, fill in manually before use.</small></p>`
        : '';
      const zNote = zWasEstimated
        ? `<p style="color:#b90"><small>Z = ${z} was estimated from cell volume ÷ (~18 Å³ × non-H atoms) — a starting guess, not a measured value. Sanity-check it (Z is usually 1, 2, 3, 4, 6, or 8).</small></p>`
        : '';
      // If the SELECTED candidate is tied with another that disagrees on
      // centrosymmetry, absences alone couldn't break that tie on their own.
      // Only worth flagging for an unreviewed default pick or to confirm an
      // E²-statistics-backed choice — a deliberate manual selection doesn't
      // need to be second-guessed here.
      const tieNote = tiedOpposite && selectionSource !== 'manual'
        ? (selectionSource === 'auto-e-stat'
            ? `<p style="color:#080"><strong>✓ ${formatHM(top.hm)} was picked over the tied ${formatHM(tiedOpposite.entry.hm)} (${tiedOpposite.entry.hall}) based on the E² statistics result above</strong> —
               absences alone couldn't distinguish them (both ${(topScore * 100).toFixed(0)}%), but the intensity-distribution check can, and did.</p>`
            : `<p style="color:#b00"><strong>⚠ This is tied at ${(topScore * 100).toFixed(0)}% with ${formatHM(tiedOpposite.entry.hm)} (${tiedOpposite.entry.hall}), which is ${topIsCentro ? 'acentric' : 'centrosymmetric'}
               where this one is ${topIsCentro ? 'centrosymmetric' : 'acentric'}.</strong> Systematic absences alone cannot break a centrosymmetric/acentric tie — Friedel's law makes both look
               identical in a diffraction pattern regardless of which one the crystal actually is. This is just the top-scoring row, not a real determination between the two — run
               <strong>Compute E² statistics</strong> on the Space Group tab, then select whichever row it points to.</p>`)
        : '';

      insResultEl.innerHTML = `
        <h2>SHELX .ins (space group ${top.number}${top.qualifier ? ':' + top.qualifier : ''}, ${formatHM(top.hm)})</h2>
        <p><small>Kept in sync automatically with whichever candidate is selected on the Space Group tab — edit the cell/formula fields or the selection there and this updates on its own.</small></p>
        ${tieNote}
        ${missingNote}
        ${zNote}
        <pre style="background:#f7f7f7; padding:0.8rem; border:1px solid #ddd; overflow-x:auto;">${ins}</pre>`;
      lastGeneratedIns = ins;
      insActionsEl.style.display = 'flex';
      insCopyNoteEl.textContent = '';
    }
    insBtn.addEventListener('click', generateIns);

    // Clipboard API needs a secure context and can silently be unavailable
    // on a file:// page (this app has no server) — fall back to the classic
    // hidden-textarea + execCommand('copy') trick rather than just failing.
    async function copyInsToClipboard() {
      if (!lastGeneratedIns) return;
      try {
        if (!navigator.clipboard || !navigator.clipboard.writeText) throw new Error('Clipboard API unavailable');
        await navigator.clipboard.writeText(lastGeneratedIns);
        insCopyNoteEl.textContent = '✓ Copied';
      } catch (e) {
        const ta = document.createElement('textarea');
        ta.value = lastGeneratedIns;
        ta.style.position = 'fixed';
        ta.style.opacity = '0';
        document.body.appendChild(ta);
        ta.select();
        try {
          document.execCommand('copy');
          insCopyNoteEl.textContent = '✓ Copied';
        } catch (e2) {
          insCopyNoteEl.textContent = 'Copy failed — please select and copy manually.';
        }
        document.body.removeChild(ta);
      }
      setTimeout(() => { insCopyNoteEl.textContent = ''; }, 2000);
    }
    insCopyBtn.addEventListener('click', copyInsToClipboard);

    insSaveBtn.addEventListener('click', () => {
      if (!lastGeneratedIns) return;
      const blob = new Blob([lastGeneratedIns], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${lastLoadedFileBaseName || 'structure'}.ins`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    });

    function runSpaceGroupMatch() {
      return new Promise((resolve) => {
        if (!lastReflections) { resolve(false); return; }

        const values = {
          a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
          alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
        };
        if (!Object.values(values).every(Number.isFinite)) {
          spaceGroupResultEl.innerHTML = `<p style="color:#b00">Please fill in all six cell constants first.</p>`;
          resolve(false);
          return;
        }

        spaceGroupResultEl.innerHTML = `<p><em>Matching against the 530-entry Hall symbol table…</em></p>`;

        setTimeout(() => {
        const classification = classifyCellMetric(values);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const targetLaueClass = determinedLaueClass || classification.laueClasses[0];
        const laueClassNote = determinedLaueClassSource === 'merging'
          ? `confirmed via merging statistics`
          : `assumed — highest metric candidate; run "Determine Laue class" above first to confirm it properly`;

        const scan = scanAllAbsences(lastReflections);
        const results = matchSpaceGroupCandidates(HALL_SYMBOLS_TABLE, targetLaueClass, scan, { monoclinicAxis });

        if (results.length === 0) {
          spaceGroupResultEl.innerHTML = `<p style="color:#b00">No candidates found for Laue class ${targetLaueClass} with the observed centering. Double-check the Laue class table above — the true Laue class may be lower symmetry than assumed here.</p>`;
          resolve(false);
          return;
        }

        lastMatchResults = results;
        selectedCandidateIndex = 0;
        selectionSource = 'default';
        insBtn.style.display = 'inline-block';
        generateIns(); // pre-populate the .ins tab with the top candidate right away

        const rows = results.slice(0, 10).map((r, i) => `
          <tr>
            <td><input type="radio" name="candidateSelect" value="${i}" ${i === 0 ? 'checked' : ''}></td>
            <td>${r.entry.number}${r.entry.qualifier ? ':' + r.entry.qualifier : ''}</td>
            <td>${formatHM(r.entry.hm)}</td>
            <td>${r.entry.hall}</td>
            <td>${(r.score * 100).toFixed(0)}%</td>
          </tr>`).join('');

        const topScore = results[0].score;
        const topIsCentro = hallSymbolIsCentrosymmetric(results[0].entry.hall);
        const tiedOpposite = results.find(r =>
          r !== results[0] && r.score === topScore && hallSymbolIsCentrosymmetric(r.entry.hall) !== topIsCentro);
        const tieWarning = tiedOpposite
          ? `<p style="color:#b00"><strong>⚠ The top ${(topScore * 100).toFixed(0)}% match is tied between an acentric and a centrosymmetric option</strong>
             (${formatHM(results[0].entry.hm)} vs. ${formatHM(tiedOpposite.entry.hm)}) — systematic absences alone can never break that particular tie, by Friedel's law.
             Run <strong>Compute E² statistics</strong> below, then pick whichever one it points to using the radio buttons — "Generate .ins" uses
             whichever row is selected, not just the top one.</p>`
          : '';

        spaceGroupResultEl.innerHTML = `
          <h2>Space group candidates</h2>
          <p><small>Assumed Laue class: ${targetLaueClass} (${laueClassNote})</small></p>
          ${tieWarning}
          <table>
            <tr><th>Use</th><th>No.</th><th>H-M symbol</th><th>Hall symbol</th><th>Match</th></tr>
            ${rows}
          </table>
          <p><small>Pick a different row any time before generating the .ins — it always uses whichever one is selected here, not just the top-scoring one.</small></p>
          <p><small>Match % = share of calculated reflection conditions that agree with the observed scan, weighted by how many reflections
          actually informed each check — a solid signal from thousands of reflections outweighs a shaky one from a sparse axial row (a check
          where the data was inconclusive — too few reflections either way — doesn't count for or against at all). A less-than-100% match on an
          otherwise clear top candidate is often just the same small-sample noise the weighting already discounts, not a wrong space group.
          Acentric candidates are included alongside centrosymmetric ones — Friedel's law makes diffraction always look centrosymmetric
          regardless of the true crystal symmetry, so absences alone usually can't rule acentric options out; a close acentric/centrosymmetric
          tie here is exactly the situation E² statistics is for.</small></p>`;

        spaceGroupResultEl.querySelectorAll('input[name="candidateSelect"]').forEach(input => {
          input.addEventListener('change', (e) => {
            selectedCandidateIndex = parseInt(e.target.value, 10);
            selectionSource = 'manual';
            generateIns(); // rebuild the .ins for the newly selected candidate right away
          });
        });
        resolve(true);
        }, 10);
      });
    }
    spaceGroupBtn.addEventListener('click', () => { runSpaceGroupMatch(); });

    const eStatBtn = document.getElementById('computeEStatBtn');
    const eStatResultEl = document.getElementById('eStatResult');

    // Composition for Wilson/E-statistics: prefer an exact whole-cell
    // count from a loaded .ins (SFAC/UNIT), otherwise derive it from the
    // "chemical formula per formula unit" field × Z (same Z logic already
    // used for .ins generation above — reused rather than re-estimated).
    function getUnitCellComposition() {
      if (lastInsComposition && !formulaManuallyEdited) {
        return { composition: lastInsComposition, source: 'exact — SFAC/UNIT from the loaded .ins' };
      }
      const formula = parseFormulaString(formulaInput.value);
      if (formula.length === 0) return null;

      let z = lastCellExtras?.z;
      let zWasEstimated = false;
      if (!z) {
        const values = {
          a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
          alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
        };
        if (!Object.values(values).every(Number.isFinite)) return null;
        const volume = cellVolume(values.a, values.b, values.c, values.alpha, values.beta, values.gamma);
        z = estimateZ(volume, formula);
        zWasEstimated = z !== null;
      }
      if (!z) return null;

      return {
        composition: formula.map(f => ({ element: f.element, count: f.count * z })),
        source: `formula × Z=${z}${zWasEstimated ? ' (Z estimated from cell volume — check it)' : ''}`,
      };
    }

    eStatBtn.addEventListener('click', () => { runEStatistics(); });

    // Runs the Wilson-plot / E² statistics step. Returns a Promise<boolean>
    // (true if it actually computed something). When called automatically
    // from autoRunPipeline (silent = true), missing prerequisites (cell,
    // composition) are just a quiet no-op — the person hasn't necessarily
    // gotten around to entering a formula yet, that's not an error. A
    // manual button click (silent = false) still explains what's missing.
    function runEStatistics(silent = false) {
      return new Promise((resolve) => {
        if (!lastReflections) { resolve(false); return; }

        const values = {
          a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
          alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
        };
        if (!Object.values(values).every(Number.isFinite)) {
          if (!silent) eStatResultEl.innerHTML = `<p style="color:#b00">Please fill in all six cell constants first.</p>`;
          resolve(false);
          return;
        }

        const compInfo = getUnitCellComposition();
        if (!compInfo) {
          if (!silent) eStatResultEl.innerHTML = `<p style="color:#b00">Need the unit cell's element composition first — either drop the matching
            SHELX .ins (exact SFAC/UNIT counts), or fill in the chemical formula field above (also needs Z, e.g. from a loaded CIF).</p>`;
          resolve(false);
          return;
        }

        eStatResultEl.innerHTML = `<p><em>Merging reflections and fitting the Wilson plot…</em></p>`;

        setTimeout(() => {
        const classification = classifyCellMetric(values);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const targetLaueClass = determinedLaueClass || classification.laueClasses[0];
        const group = generateGroup(getLaueClassGenerators(targetLaueClass, monoclinicAxis));

        const stats = computeEStatistics(lastReflections, values, compInfo.composition, group);

        if (stats.error) {
          eStatResultEl.innerHTML = `<p style="color:#b00">${stats.error}</p>`;
          resolve(false);
          return;
        }

        const ref = WILSON_REFERENCE;
        const verdictColor = stats.verdict === 'centrosymmetric' ? '#080' : '#06c';

        // If the space-group candidates included a centrosymmetric/acentric
        // tie at the top score, this verdict IS the evidence needed to
        // break it — apply it directly instead of leaving the person to
        // notice and click a radio button themselves.
        let autoSelectNote = '';
        if (lastMatchResults && lastMatchResults.length > 0) {
          const topScore = lastMatchResults[0].score;
          const tiedGroup = lastMatchResults.filter(m => m.score === topScore);
          const matchesVerdict = m => hallSymbolIsCentrosymmetric(m.entry.hall) === (stats.verdict === 'centrosymmetric');
          const resolved = tiedGroup.find(matchesVerdict);
          const wasActuallyTied = tiedGroup.some(m => !matchesVerdict(m));
          if (resolved && wasActuallyTied) {
            const newIndex = lastMatchResults.indexOf(resolved);
            selectionSource = 'auto-e-stat'; // this pick is now backed by E² statistics either way
            if (newIndex !== selectedCandidateIndex) {
              selectedCandidateIndex = newIndex;
              generateIns(); // rebuild the .ins for the now-resolved candidate right away
              const radio = document.querySelector(`input[name="candidateSelect"][value="${newIndex}"]`);
              if (radio) radio.checked = true;
              autoSelectNote = `<p style="color:#080"><strong>✓ Candidate selection updated above:</strong> switched the space-group pick to
                ${formatHM(resolved.entry.hm)} (${resolved.entry.hall}) to match this ${stats.verdict} result — the SHELX .ins tab has been updated to match.</p>`;
            } else {
              autoSelectNote = `<p><small>This confirms the already-selected candidate (${formatHM(resolved.entry.hm)}).</small></p>`;
            }
          }
        }

        eStatResultEl.innerHTML = `
          <h2>E² statistics (Wilson plot)</h2>
          <p><small>Composition: ${compInfo.source} · merged under Laue class ${targetLaueClass}${determinedLaueClassSource === 'merging' ? '' : ' (assumed, not yet confirmed above)'}
          · ${stats.uniqueReflections} unique reflections (${stats.usableReflections} usable)</small></p>
          <table>
            <tr><th></th><th>this data</th><th>ref. centrosymmetric</th><th>ref. acentric</th></tr>
            <tr><td>&lt;|E²−1|&gt;</td><td><strong>${stats.meanE2minus1.toFixed(3)}</strong></td><td>${ref.centrosymmetric.meanE2minus1}</td><td>${ref.acentric.meanE2minus1}</td></tr>
            <tr><td>&lt;|E|&gt;</td><td>${stats.meanE.toFixed(3)}</td><td>${ref.centrosymmetric.meanE}</td><td>${ref.acentric.meanE}</td></tr>
            <tr><td>fraction |E|&gt;3</td><td>${(stats.fractionEgt3 * 100).toFixed(3)}%</td><td>${(ref.centrosymmetric.fractionEgt3 * 100).toFixed(2)}%</td><td>${(ref.acentric.fractionEgt3 * 100).toFixed(2)}%</td></tr>
            <tr><td>&lt;E²&gt;</td><td>${stats.meanE2.toFixed(3)}</td><td colspan="2"><small>should be ≈1.000 if the Wilson scaling fit well</small></td></tr>
          </table>
          <p style="color:${verdictColor}"><strong>Best match: ${stats.verdict}</strong> <small>(margin ${(stats.margin * 100).toFixed(0)}% between the two reference distances — near 0% means right on the fence)</small></p>
          ${autoSelectNote}
          <p><small>Wilson fit: B = ${stats.B.toFixed(2)} Å² (overall isotropic falloff), scale k = ${stats.k.toFixed(3)}. If &lt;E²&gt; above is far from 1.000, treat the verdict cautiously —
          common real causes are a dominant heavy atom (breaks the "randomly distributed atoms" Wilson assumption, especially at low resolution), twinning, or non-ideal data scaling; none of those are
          errors in this calculation, they're properties of the data itself. This is a different question from the absences scan above — Friedel's law makes every diffraction pattern look centrosymmetric
          regardless of the crystal's true symmetry, so this is the one of the two that can actually distinguish an acentric structure from a centrosymmetric one.</small></p>
          <div id="wilsonPlotContainer"></div>`;

        renderWilsonPlot(document.getElementById('wilsonPlotContainer'), stats.wilsonPoints, stats.B, stats.k);
        resolve(true);
        }, 10);
      });
    }

    input.addEventListener('change', (e) => handleFile(e.target.files[0]));

    dropZone.addEventListener('dragover', (e) => {
      e.preventDefault();
      dropZone.classList.add('dragover');
    });
    dropZone.addEventListener('dragleave', () => {
      dropZone.classList.remove('dragover');
    });
    dropZone.addEventListener('drop', (e) => {
      e.preventDefault();
      dropZone.classList.remove('dragover');
      const file = e.dataTransfer.files[0];
      handleFile(file);
    });
