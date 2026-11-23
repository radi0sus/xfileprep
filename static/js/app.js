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
    const laueResultEl = document.getElementById('laueResult');
    let lastReflections = null;
    let determinedLaueClass = null;
    let lastRawText = null;
    let lastCellExtras = null;
    let determinedLaueClassSource = null; // 'cif' | 'merging' | null
    let lastMatchResults = null;
    let selectedCandidateIndex = 0;
    let selectionSource = 'default'; // 'default' | 'manual' | 'auto-e-stat' | 'cif'
    let lastInsComposition = null; // exact whole-unit-cell composition, from a loaded .ins's SFAC/UNIT
    let lastCellEsd = null; // {a,b,c,alpha,beta,gamma} esds from whichever source (CIF/.ins/.p4p) provided them
    // Tracks whether the person has typed into the formula field themselves
    // since the last file load. When true, a manually-typed formula (× Z)
    // is used for the .ins SFAC/UNIT instead of an exact whole-cell
    // composition from a previously loaded .ins/.p4p — otherwise editing
    // the formula field had no visible effect whenever a composition was
    // already known, since composition always took priority (see Bug 8).
    let formulaManuallyEdited = false;
    let lastHklExport = null; // {text, count, source} — only set when reflections came out of a CIF (never for a raw SHELX .hkl)
    // Which file supplied what — shown in a strip above the tabs so it stays
    // visible on every tab (the app otherwise "forgets" earlier drops visually).
    const loadedFiles = { reflections: null, cell: null };
    let awaitingJump = false; // pipeline is ready, but we stayed on tab 1 because the formula is still missing
    let batchDepth = 0;      // >0 while several files are processed in one go
    let pendingRun = false;  // pipeline run requested during a batch, executed once at its end
    const loadedFilesEl = document.getElementById('loadedFiles');
    function noteLoaded(kind, name, extra) {
      loadedFiles[kind] = { name: name || '(unnamed)', extra: extra || '' };
      renderLoadedFiles();
    }
    function renderLoadedFiles() {
      if (!loadedFilesEl) return;
      const chip = (label, item) => item
        ? `<span class="fileChip ok">✔ ${label}: <strong>${item.name.replace(/[<>&]/g, '')}</strong>${item.extra ? ` <small>${item.extra}</small>` : ''}</span>`
        : `<span class="fileChip missing">○ ${label}: not loaded yet</span>`;
      loadedFilesEl.innerHTML = chip('Reflections', loadedFiles.reflections) + chip('Cell', loadedFiles.cell)
        + (awaitingJump
          ? ` <span class="fileChip hint">Results ready. No formula yet — drop the CIF/.ins/.p4p with it, or
              <button type="button" data-go="tab-spacegroup">continue to Space Group →</button></span>`
          : '');
    }
    if (loadedFilesEl) loadedFilesEl.addEventListener('click', (e) => {
      const go = e.target.closest && e.target.closest('[data-go]');
      if (!go) return;
      awaitingJump = false;
      renderLoadedFiles();
      window.switchTab(go.dataset.go);
    });
    renderLoadedFiles();
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

    // Purely informational note for the Space Group candidate table: does
    // THIS row's cell/hkl match what was actually measured, or would
    // selecting it need a transform first? Never affects the .ins output.
    //
    // IMPORTANT: "needs a transform" is NOT the same question as "is this
    // row's qualifier non-null". The candidate list can (and often does)
    // contain several rows for the SAME space-group number — e.g. both
    // 45:cab (I 2 c b) and 45:null (I b a 2) — each scored against the
    // SAME as-measured reflections. Only whichever of them actually scores
    // highest is the one the raw, untransformed cell/hkl natively matches;
    // ALL other rows for that number — including the qualifier-null
    // "standard" one, if it isn't the top scorer — need the transform
    // applied just as much as any differently-labeled row does. Treating
    // "qualifier === null" as "no transform needed" (the earlier version
    // of this function) was exactly this bug: picking the standard-named
    // but lower-scoring row showed no matrix at all, when a real dataset
    // (I b a 2 scoring below I 2 c b for the same reflections) needs
    // exactly the same transform to be written correctly as the non-
    // standard row would.
    // Index (in `results`) of the space group the loaded CIF reports, or -1.
    // Hall symbol first; else H-M + IT number. A setting collapsed away by
    // deduplicateByConditionSignature is matched by number + identical conditions.
    function cifReportedCandidateIndex(results) {
      const cell = lastCellExtras;
      if (!cell) return -1;
      let entry = null;
      if (cell.reportedHallSymbol) entry = HALL_SYMBOLS_TABLE.find(e => e.hall === cell.reportedHallSymbol.trim());
      if (!entry && cell.reportedSpaceGroupNumber != null && cell.reportedSpaceGroup) {
        const norm = x => String(x).toLowerCase().replace(/\s+/g, '').replace(/:/g, '');
        entry = HALL_SYMBOLS_TABLE.find(e => e.number === cell.reportedSpaceGroupNumber && norm(e.hm) === norm(cell.reportedSpaceGroup));
      }
      if (!entry) return -1;
      const direct = results.findIndex(r => r.entry.hall === entry.hall);
      if (direct >= 0) return direct;
      const sig = conditionSignature(entry.hall);
      return results.findIndex(r => r.entry.number === entry.number && conditionSignature(r.entry.hall) === sig);
    }

    function currentCellValues() {
      return {
        a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
      };
    }

    // Monoclinic only: if the target cell choice is clearly more oblique
    // than the measured one (e.g. P 21/n with beta 95 -> P 21/c with beta
    // 125), say so — IUCr/ITA recommend the least oblique cell choice, so
    // the as-measured setting is normally the better one to refine in.
    function obliquenessNote(entry, matrix, cell) {
      if (!(entry.number >= 3 && entry.number <= 15) || !cellIsUsable(cell)) return '';
      const before = monoclinicAngle(cell);
      const after = monoclinicAngle(transformCellParameters(cell, matrix));
      if (Math.abs(after - 90) <= Math.abs(before - 90) + 0.5) return '';
      return `<br><small style="color:#b90" title="IUCr/ITA convention: choose the monoclinic cell choice with the least oblique angle. The space group is the same either way; only the description differs.">⚠ more oblique: ${after.toFixed(2)}° instead of ${before.toFixed(2)}° — the measured setting is the recommended one</small>`;
    }

    function describeSetting(entry) {
      if (!lastMatchResults) return '';
      const sameNumber = lastMatchResults.filter(r => r.entry.number === entry.number);
      const native = sameNumber.reduce((best, r) => (r.score > best.score ? r : best), sameNumber[0]).entry;

      if (native.hall === entry.hall) {
        return '<span style="color:#888">matches as measured</span>';
      }
      const cell = currentCellValues();
      const matrix = deriveTransformMatrix(native.hall, entry.hall, cell);
      if (!matrix) {
        return '<span style="color:#b90" title="No matrix in the known pool reproduces this row\'s operators from the as-measured description — e.g. hexagonal vs. rhombohedral axes, not supported">different axis convention — matrix not found</span>';
      }
      if (isIdentity3(matrix)) {
        // Origin choice 1 vs 2 (F d -3 m, P n n n, I 41/a m d, ...): the same
        // lattice and indexing, only the origin of the symmetry operators moves.
        return `<span style="color:#888" title="Differs from ${formatHM(native.hm)} only by an origin shift: same cell, same hkl, no HKLF matrix. The .ins gets this row's symmetry operators directly.">matches as measured <small>(other origin choice)</small></span>`;
      }
      const esdNote = isPurePermutation(matrix)
        ? ''
        : '<br><small style="color:#b90">⚠ cell-choice shear — esds not recalculated</small>';
      return `<span style="color:#069" title="Your cell/hkl as measured natively match ${formatHM(native.hm)}. To use this row (${formatHM(entry.hm)}) instead, transform with this matrix.">as measured (${formatHM(native.hm)}) → <strong>${formatHM(entry.hm)}</strong><br><code style="font-size:0.85em">${formatHKLFMatrix(matrix)}</code>${esdNote}${obliquenessNote(entry, matrix, cell)}</span>`;
    }

    const absencesResultEl = document.getElementById('absencesResult');

    const VERDICT_STYLE = {
      holds: 'color:#080', 
      violated: 'color:#aaa',
      inconclusive: 'color:#b90',
    };
    // Display wording only — the internal verdict values stay holds/violated/inconclusive.
    const VERDICT_LABEL = {
      holds: 'satisfied',
      violated: 'violated',
      inconclusive: 'undetermined',
    };
    function verdictCell(r) {
      const style = VERDICT_STYLE[r.verdict] || '';
      const title = r.reason ? ` title="${r.reason.replace(/"/g, '&quot;')}"` : '';
      const label = VERDICT_LABEL[r.verdict] || r.verdict;
      return `<td style="${style}"${title}>${label} <small title="share of forbidden reflections with I/σ ≥ 3">(${r.violationPercent.toFixed(1)}%)</small></td>`;
    }

    // Current cell fields as numbers, or null if any is missing (plots then fall back to an index grid).
    function readCellOrNull() {
      const v = {};
      for (const k of ['a', 'b', 'c', 'alpha', 'beta', 'gamma']) v[k] = parseFloat(cellFields[k].value);
      return Object.values(v).every(Number.isFinite) ? v : null;
    }

    function renderLaueTab() {
      if (!lastReflections) return;
      document.getElementById('laueVizSection').style.display = 'block';
      document.getElementById('laueVizEmpty').style.display = 'none';
      // Laue group for the optional "fill in symmetry equivalents" display
      // (same class and unique axis as the layer viewer and the shell table).
      const cellForGroup = readCellOrNull();
      let laueGroup = null, laueForFill = determinedLaueClass;
      try {
        if (cellForGroup) {
          const classification = classifyCellMetric(cellForGroup);
          const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
          laueForFill = determinedLaueClass || classification.laueClasses[0];
          laueGroup = generateGroup(getLaueClassGenerators(laueForFill, monoclinicAxis));
        }
      } catch (err) { laueGroup = null; }
      renderLaueOverview(document.getElementById('laueVizContainer'), document.getElementById('laueVizNote'),
        lastReflections, cellForGroup, determinedLaueClass, determinedLaueClassSource, laueGroup);
      renderShellStats();
      renderTwinCheck();
      renderLayerViewer();
    }

    // Drops the reflections that the OBSERVED lattice centering forbids (zero by
    // construction, they would distort every statistic). Independent of the
    // space-group match, so the E² and shell statistics do not wait for it.
    // Absence scan, computed once per loaded data set (the E² tab, the shell
    // table and the layer viewer all need it).
    let _scanCache = { reflections: null, scan: null };
    function cachedAbsenceScan() {
      if (_scanCache.reflections !== lastReflections) {
        _scanCache = { reflections: lastReflections, scan: scanAllAbsences(lastReflections) };
      }
      return _scanCache.scan;
    }

    // Drops the reflections that the OBSERVED lattice centering forbids (zero by
    // construction, they would distort every statistic). Independent of the
    // space-group match, so the E² and shell statistics do not wait for it.
    // `allowedTest` additionally excludes positions that the absence scan found
    // to be systematically forbidden (screw axes, glides): it defines which
    // reflections COUNT as possible for the completeness and the shell
    // boundaries. Data that were already cleaned of absences (deposited CIF
    // /.hkl files) would otherwise look incomplete.
    function observedLatticeFilter() {
      const scan = cachedAbsenceScan();
      const centName = observedCenteringName(scan.centering);
      const lat = filterByLattice(lastReflections, centName[0]);
      const allowedTest = r => (!lat.test || lat.test(r)) && !isForbiddenByHeldConditions(scan, r.h, r.k, r.l);
      return { ...lat, allowedTest, scan };
    }

    // Layer viewer (layer-viewer.js): same Laue class / observed lattice as the
    // shell table; the absence scan keeps systematically forbidden positions
    // from being drawn as gaps.
    function renderLayerViewer() {
      const el = document.getElementById('layerViewerContainer');
      if (!el) return;
      if (!lastReflections) { el.innerHTML = ''; return; }
      const cell = readCellOrNull();
      if (!cell) { el.innerHTML = ''; return; }
      try {
        const classification = classifyCellMetric(cell);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const laue = determinedLaueClass || classification.laueClasses[0];
        const group = generateGroup(getLaueClassGenerators(laue, monoclinicAxis));
        const lat = observedLatticeFilter();
        initLayerViewer(el, { reflections: lastReflections, cell, group, latticeTest: lat.allowedTest, scan: lat.scan, laue, source: determinedLaueClassSource });
      } catch (err) {
        el.innerHTML = `<p style="color:#b00"><small>Could not build the layer viewer: ${err.message}</small></p>`;
      }
    }

    // Resolution-shell merging statistics (shell-stats.js) under the Laue
    // class from merging (or the highest-symmetry candidate if none was
    // confirmed). Centering-absent reflections of the OBSERVED lattice are
    // dropped first (they are zero by construction and would drag <I/σ> and
    // the R values down). Information only: nothing is written to the .ins.
    function renderShellStats() {
      const el = document.getElementById('shellStatsContainer');
      if (!el) return;
      if (!lastReflections) { el.innerHTML = ''; return; }
      const cell = readCellOrNull();
      if (!cell) { el.innerHTML = '<p><small>Needs all six cell constants.</small></p>'; return; }
      try {
        const classification = classifyCellMetric(cell);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const laue = determinedLaueClass || classification.laueClasses[0];
        const group = generateGroup(getLaueClassGenerators(laue, monoclinicAxis));
        const lat = observedLatticeFilter();
        const stats = computeShellStatistics(lat.reflections, cell, group, 10, lat.allowedTest);
        renderShellStatistics(el, stats, { laue, source: determinedLaueClassSource, lattice: lat.note });
      } catch (err) {
        el.innerHTML = `<p style="color:#b00"><small>Could not compute shell statistics: ${err.message}</small></p>`;
      }
    }

    // Twin check (twin-check.js) under the current Laue class. The suggested
    // TWIN/BASF lines go into the .ins only if the person ticks the box.
    let lastTwinCheck = null;      // result of checkMerohedralTwinning, or null
    let twinWriteToIns = false;    // checkbox state, survives re-renders of the section
    // Hint on the E² tab: twinning narrows the intensity distribution (N(z),
    // <|E²-1|>), so the centro/acentric verdict there is not reliable.
    function renderE2TwinNote() {
      const el = document.getElementById('eStatTwinNote');
      if (!el) return;
      const laws = lastTwinCheck && lastTwinCheck.applicable ? lastTwinCheck.laws.filter(L => L.verdict === 'twinned' || L.verdict.startsWith('perfect')) : [];
      if (!laws.length) { el.innerHTML = ''; return; }
      const best = laws.sort((a, b) => b.alpha - a.alpha)[0];
      el.innerHTML = `<p style="color:#b60"><strong>⚠ Twinning found</strong> (Twinning tab: <code>TWIN ${best.shelxMatrix.flat().join(' ')}</code>, α ≈ ${best.alpha.toFixed(2)}).
        Twinning averages the intensities of twin-related reflections and narrows their distribution: ⟨|E²−1|⟩ drops below the acentric value
        and N(z) moves towards the twinned-acentric curve. The centrosymmetric/acentric verdict below is therefore unreliable — an "acentric" can
        come from the twin alone, while a "centrosymmetric" despite the twin points to pseudo-symmetry (e.g. Z′ > 1, heavy atoms on special positions).</p>`;
    }

    // What the twin checks found worth flagging: a twinned law and/or an
    // inversion twin (x >= 0.1 with enough expected signal). Null if nothing.
    function twinFindings() {
      const out = [];
      if (lastTwinCheck && lastTwinCheck.applicable) {
        const best = lastTwinCheck.laws.filter(L => L.verdict === 'twinned' || L.verdict.startsWith('perfect')).sort((x, y) => y.alpha - x.alpha)[0];
        if (best) out.push({ kind: 'law', text: `<code>TWIN ${best.shelxMatrix.flat().join(' ')}</code>, α ≈ ${best.alpha.toFixed(2)}` });
      }
      const t = lastAnomalousCheck && lastAnomalousCheck.twin;
      const sel = lastMatchResults && lastMatchResults[selectedCandidateIndex] ? lastMatchResults[selectedCandidateIndex].entry : null;
      // an inversion twin only exists in an acentric space group
      if (t && t.x >= 0.1 && t.expectedChi2 >= 1.5 && sel && !hallSymbolIsCentrosymmetric(sel.hall)) out.push({ kind: 'inversion', text: `inversion twin, x ≈ ${t.x.toFixed(2)}${t.x >= 0.35 ? ' (or centrosymmetric)' : ''}` });
      return out.length ? out : null;
    }

    // Smallest integer direction with the first non-zero component positive
    // ([0 -2 0] -> [0 1 0]; centring halves are doubled first).
    function reducedDirection(v) {
      let w = v.map(x => Math.round(x * 6)); // halves and thirds (centring vectors) become integers
      const gcd = (a, b) => (b ? gcd(b, a % b) : Math.abs(a));
      const g = w.reduce((a, x) => gcd(a, x), 0) || 1;
      w = w.map(x => x / g);
      const first = w.find(x => x !== 0);
      return first < 0 ? w.map(x => -x || 0) : w;
    }

    // Tab badge + banner under the tab bar, visible on every tab.
    function updateTwinAlerts() {
      const f = twinFindings();
      const btn = document.getElementById('twinTabBtn');
      if (btn) {
        btn.classList.toggle('twinAlert', !!f);
        btn.textContent = f ? '3. Twinning ⚠' : '3. Twinning';
      }
      const banner = document.getElementById('twinBanner');
      if (!banner) return;
      if (!f) { banner.style.display = 'none'; banner.innerHTML = ''; return; }
      banner.innerHTML = `<strong>⚠ Twinning found:</strong> ${f.map(x => x.text).join(' · ')} — <a id="twinBannerLink">see the Twinning tab</a>.`;
      banner.style.display = 'block';
      const link = document.getElementById('twinBannerLink');
      if (link) link.addEventListener('click', () => window.switchTab && window.switchTab('tab-twin'));
    }

    function renderTwinCheck() {
      const el = document.getElementById('twinCheckContainer');
      if (!el) return;
      lastTwinCheck = null;
      try { renderTwinCheckInner(el); } finally { renderE2TwinNote(); updateTwinAlerts(); }
    }

    function renderTwinCheckInner(el) {
      lastAnomalousCheck = null;
      if (!lastReflections) { el.innerHTML = ''; return; }
      const cell = readCellOrNull();
      if (!cell) { el.innerHTML = '<p><small>Needs all six cell constants.</small></p>'; return; }
      let html = '';
      try {
        const classification = classifyCellMetric(cell);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const laue = determinedLaueClass || classification.laueClasses[0];
        const res = checkMerohedralTwinning(lastReflections, cell, laue, monoclinicAxis, observedCenteringName(cachedAbsenceScan().centering));
        lastTwinCheck = res;
        // lattice 2-fold axes (Le Page) — where the twin laws come from
        const folds = res.folds || [];
        html += `<h3>Lattice symmetry (Le Page)</h3>
          <p><small>2-fold axes of the lattice, found independent of the cell choice: a lattice direction [uvw] (nearly) perpendicular to a lattice plane (hkl);
            the obliquity δ is the angle between them. Laue class ${res.laue}; lattice symmetry <strong>${res.holohedry}</strong> (${res.holohedryNote}).</small></p>`;
        html += folds.length
          ? `<table class="shellTable"><thead><tr><th>Axis [uvw]</th><th>⊥ plane (hkl)</th><th>Obliquity δ</th><th>In the Laue class?</th></tr></thead><tbody>${
              folds.map(f => `<tr><td>[${reducedDirection(f.t).join(' ')}]</td><td>(${reducedDirection(f.h).join(' ')})</td><td>${f.delta.toFixed(3)}°</td>
                <td>${f.inLaue ? 'yes' : (f.used ? '<strong>no — possible twin law</strong>' : 'no (inconsistent with the others, not used)')}</td></tr>`).join('')}</tbody></table>`
          : '<p><small>No 2-fold axis within the tolerance.</small></p>';
        html += '<h3>Twin laws beyond the Laue class (H test)</h3>';
        if (!res.applicable) {
          html += `<p><small>No twin law to test: ${res.reason}.</small></p>`;
        } else {
          const fmtM = M => M.flat().join(' ');
          const rows = res.laws.map(L => {
            const color = L.verdict === 'twinned' ? '#b60' : L.verdict.startsWith('perfect') ? '#b00' : '#080';
            const dev = L.deviation ? `${L.deviation.angle.toFixed(2)}°${L.deviation.length > 0.0005 ? ` / ${(100 * L.deviation.length).toFixed(2)} %` : ''}` : '—';
            return `<tr><td><code>TWIN ${fmtM(L.shelxMatrix)}</code></td><td>${dev}</td><td>${L.pairs}</td>
              <td>${L.meanH === null ? '—' : L.meanH.toFixed(3)}</td>
              <td>${L.alpha === null ? '—' : `${L.alpha.toFixed(2)} <small>± ${L.alphaSe.toFixed(3)}</small>`}</td>
              <td style="color:${color}">${L.verdict}</td></tr>`;
          }).join('');
          const best = res.laws.filter(L => L.verdict === 'twinned' || L.verdict.startsWith('perfect')).sort((x, y) => y.alpha - x.alpha)[0];
          html += `
            <table class="shellTable"><thead><tr><th>Twin law (SHELX)</th><th title="Le Page obliquity: angle between the lattice direction and the reciprocal-lattice normal of the 2-fold axis — 0 = exact (merohedry), small = pseudo-merohedry">Obliquity δ</th>
              <th title="Acentric twin-related pairs, I1+I2 > ${res.laws[0] ? res.laws[0].selection : 10}σ">Pairs</th><th title="Untwinned: 0.5, perfect twin: 0">⟨H⟩</th><th>Twin fraction α</th><th>Verdict</th></tr></thead>
              <tbody>${rows}</tbody></table>
            <p><small>α = ½ − ⟨H⟩ is an estimate: pseudo-symmetry raises it, noise lowers it — use it as start value for BASF and refine.
              α ≥ 0.40 cannot be told apart from a crystal that simply has the higher symmetry: compare R<sub>int</sub> of ${res.holohedry} on the Laue tab.</small></p>
            <div id="nhPlotContainer"></div>
            <p><small>N(H): straight lines are ideal twins with fraction α; an observed curve that bends away from them points to pseudo-symmetry rather than a clean twin.</small></p>
            ${best ? `<p style="color:#b60"><small>Twinning narrows the intensity distribution, so the centrosymmetric/acentric verdict on the E² tab is unreliable for this data set.</small></p>` : ''}
            ${best ? `<p><label><input type="checkbox" id="twinToIns" ${twinWriteToIns ? 'checked' : ''}> <span>write <code>TWIN ${fmtM(best.shelxMatrix)}</code> with <code>BASF ${best.alpha.toFixed(2)}</code> to the .ins</span></label></p>` : ''}`;
        }
      } catch (err) {
        html += `<p style="color:#b00"><small>Could not run the twin-law check: ${err.message}</small></p>`;
      }
      html += renderInversionTwinHtml(cell);
      el.innerHTML = html;
      const nhEl = document.getElementById('nhPlotContainer');
      if (nhEl && lastTwinCheck && lastTwinCheck.applicable) renderNHPlot(nhEl, lastTwinCheck.laws);
      const cb = document.getElementById('twinToIns');
      if (cb) cb.addEventListener('change', () => { twinWriteToIns = cb.checked; if (lastMatchResults) generateIns(); });
      const cbInv = document.getElementById('invTwinToIns');
      if (cbInv) cbInv.addEventListener('change', () => { invTwinWriteToIns = cbInv.checked; if (lastMatchResults) generateIns(); });
      if (lastMatchResults) generateIns();
    }

    // Inversion (racemic) twinning from the Bijvoet pairs (anomalous-check.js).
    // Needs an acentric point group to merge under: the selected candidate if
    // acentric, else the best-scoring acentric one (e.g. P 21 when P 21/m is selected).
    let lastAnomalousCheck = null;
    let invTwinWriteToIns = false;
    // f″ for the Bijvoet check: values from a loaded CIF take priority (they
    // belong to the refinement); elements the CIF does not list — or all of
    // them, without a CIF — come from the FFAST table for the selected anode
    // (Mo/Cu/Ag Kα). Returns {fpp: {El: f″} | null, source, note}.
    function effectiveFpp(composition) {
      const anode = anodeForWavelength(currentWavelength().value);
      const fromCif = lastFpp ? { ...lastFpp } : null;
      if (!anode) {
        return fromCif
          ? { fpp: fromCif, source: 'from the CIF' }
          : { fpp: null, source: '', note: 'the wavelength is not Mo, Cu or Ag Kα' };
      }
      const fpp = fromCif || {};
      const filled = [];
      for (const c of composition || []) {
        if (fpp[c.element] !== undefined) continue;
        const info = elementInfo(c.element);
        if (info && info.abs) { fpp[c.element] = info.abs[anode][1]; filled.push(c.element); }
      }
      if (!fromCif) return { fpp: filled.length ? fpp : null, source: `from NIST FFAST at ${anode} Kα`, note: filled.length ? '' : 'no composition yet' };
      return { fpp, source: filled.length ? `from the CIF; ${filled.join(', ')} from NIST FFAST at ${anode} Kα` : 'from the CIF' };
    }
    function renderInversionTwinHtml(cell) {
      let html = '<h3>Inversion (racemic) twinning — Bijvoet pairs</h3>';
      if (!lastMatchResults || !lastMatchResults.length) return html + '<p><small>Available once the space-group candidates are matched.</small></p>';
      const sel = lastMatchResults[selectedCandidateIndex];
      const acentric = r => !hallSymbolIsCentrosymmetric(r.entry.hall);
      const ref = acentric(sel) ? sel
        : (lastMatchResults.find(r => r.score === sel.score && acentric(r)) || lastMatchResults.find(acentric));
      if (!ref) return html + '<p><small>No acentric candidate — inversion twinning does not apply.</small></p>';
      try {
        const ops = generateSpaceGroupOperators(ref.entry.hall);
        const pointGroup = [...new Map(ops.map(o => [o.R.map(r => r.join(',')).join(';'), o.R])).values()];
        const comp = getUnitCellComposition();
        const fppInfo = effectiveFpp(comp ? comp.composition : null);
        const res = checkAnomalousSignal(lastReflections, cell, pointGroup, fppInfo.fpp, comp ? comp.composition : null);
        lastAnomalousCheck = res;
        if (!res.enough) return html + `<p><small>Only ${res.pairs} Bijvoet pairs (Friedel mates measured separately) — not enough. Merged data with Friedel pairs averaged carry no anomalous information.</small></p>`;
        const pgNote = ref === sel ? `point group of ${formatHM(ref.entry.hm)}` : `point group of the acentric candidate ${formatHM(ref.entry.hm)} (${formatHM(sel.entry.hm)} is centrosymmetric)`;
        const sigmaNote = res.calibrated
          ? (res.sigmaScale > 1.05 ? `σ scaled ×${Math.sqrt(res.sigmaScale).toFixed(2)} to match the scatter of repeated measurements` : 'σ consistent with repeated measurements')
          : 'σ not calibrated (merged data)';
        html += `<p>${res.pairs} Bijvoet pairs (${pgNote}; ${sigmaNote}).<br>
          Observed anomalous signal: χ² = <strong>${res.chi2.toFixed(2)}</strong> <small>(1 = noise only)</small>`;
        if (res.twin) {
          const t = res.twin;
          const color = t.verdict.startsWith('enantiopure') ? '#080' : t.verdict.startsWith('partial') ? '#b60' : t.verdict.startsWith('racemic') ? '#b00' : '#888';
          html += ` · expected without inversion twinning: χ² ≈ ${t.expectedChi2.toFixed(1)}<br>
            Estimated inversion-twin fraction: <strong>x ≈ ${t.x.toFixed(2)}</strong> <small>± ${t.xSe.toFixed(2)}</small> —
            <span style="color:${color}">${t.verdict}</span></p>
            <p><small>From observed/expected Bijvoet differences ≈ |1 − 2x| (f″ ${fppInfo.source}, composition ${comp ? comp.source : ''}). Like the Flack parameter but without a model, so without sign:
            x and 1 − x look the same, and a 50:50 racemic twin is indistinguishable from a centrosymmetric structure.
            ${lastTwinCheck && lastTwinCheck.applicable && lastTwinCheck.laws.some(L => L.verdict === 'twinned') ? ' Combined with the twin law above the two effects mix — rough indication only.' : ''}
            Refine the Flack parameter (or BASF) to get the real value.</small></p>`;
          if (t.x >= 0.1 && t.expectedChi2 >= 1.5) {
            html += ref === sel
              ? `<p><label><input type="checkbox" id="invTwinToIns" ${invTwinWriteToIns ? 'checked' : ''}> <span>write the inversion twin (<code>TWIN</code>, <code>BASF ${t.x.toFixed(2)}</code>) to the .ins</span></label></p>`
              : `<p><small>An inversion twin only exists for an acentric space group — select ${formatHM(ref.entry.hm)} on the Space Group tab to write it to the .ins.</small></p>`;
          }
        } else {
          const why = !fppInfo.fpp ? `no f″ values — they come from a CIF (<code>_atom_type_scat_dispersion_imag</code>) or, for Mo/Cu/Ag Kα, from the built-in table${fppInfo.note ? ` (${fppInfo.note})` : ''}`
            : !comp ? 'no composition (formula + Z or a .ins)'
            : `f″ missing for ${res.missingFpp.join(', ') || 'the anomalous scatterers'}`;
          html += `</p><p><small>The expected signal cannot be computed (${why}), so a racemic twin cannot be told apart from a structure with only light atoms.
            ${res.acentricEvidence === 'strong' ? 'The clear anomalous signal shows the structure is acentric and not a 50:50 racemic twin.' : ''}</small></p>`;
        }
        if (res.acentricEvidence === 'strong' && hallSymbolIsCentrosymmetric(sel.entry.hall)) {
          html += `<p style="color:#b00"><strong>⚠ Friedel's law is clearly broken (χ² = ${res.chi2.toFixed(1)}), which no centrosymmetric structure can do — the selected ${formatHM(sel.entry.hm)} is probably wrong;
            consider ${formatHM(ref.entry.hm)}.</strong> <small>(Unless absorption or scaling errors are large.)</small></p>`;
        }
      } catch (err) {
        html += `<p style="color:#b00"><small>Could not run the Bijvoet check: ${err.message}</small></p>`;
      }
      return html;
    }

    // Twin instructions for the .ins (ticked boxes, significant results), or null:
    // a twin law with BASF α, the inversion twin with BASF x, or both combined
    // as TWIN R -4 with the three component fractions.
    function twinForIns() {
      let law = null, inv = null;
      if (twinWriteToIns && lastTwinCheck && lastTwinCheck.applicable) {
        const best = lastTwinCheck.laws.filter(L => L.verdict === 'twinned' || L.verdict.startsWith('perfect')).sort((x, y) => y.alpha - x.alpha)[0];
        if (best) law = best;
      }
      const selEntry = lastMatchResults && lastMatchResults[selectedCandidateIndex] ? lastMatchResults[selectedCandidateIndex].entry : null;
      if (invTwinWriteToIns && selEntry && !hallSymbolIsCentrosymmetric(selEntry.hall)
          && lastAnomalousCheck && lastAnomalousCheck.twin && lastAnomalousCheck.twin.x >= 0.1) inv = lastAnomalousCheck.twin.x;
      const minus1 = [[-1, 0, 0], [0, -1, 0], [0, 0, -1]];
      if (law && inv !== null) {
        const a = law.alpha, x = inv; // components: R (a(1-x)), -1 (x(1-a)), R·-1 (a·x)
        return { matrix: law.shelxMatrix, n: -4, basf: [a * (1 - x), x * (1 - a), a * x] };
      }
      if (law) return { matrix: law.shelxMatrix, n: 2, basf: [law.alpha] };
      if (inv !== null) return { matrix: minus1, n: 2, basf: [inv] };
      return null;
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
          results.map((r, i) => `<tr><td>${i === 0 ? (zone === 'hhl' ? 'hhl (diagonal plane)' : zone) : ''}</td><td>${r.name}</td>${verdictCell(r)}</tr>`)
        ).join('');

        absencesResultEl.innerHTML = `
          <h2>Reflection conditions (from systematic absences)</h2>

          <h3>Lattice centering (whole dataset)</h3>
          <table><tr><th>Type</th><th>Status</th></tr>${centeringRows}</table>

          <h3>Axial conditions (screw axes)</h3>
          <table>
            <tr><th>Axis</th><th>n refl.</th><th title="reflections are allowed only if the index along this axis is a multiple of 2">index = 2n</th><th title="reflections are allowed only if the index along this axis is a multiple of 3">index = 3n</th><th title="reflections are allowed only if the index along this axis is a multiple of 4">index = 4n</th><th title="reflections are allowed only if the index along this axis is a multiple of 6">index = 6n</th></tr>
            ${axialRows}
          </table>

          <h3>Zonal conditions (glide planes)</h3>
          <table><tr><th>Zone</th><th>Condition</th><th>Status</th></tr>${zonalRows}</table>

          <details class="explain"><summary>How to read these results</summary><p><small>Percentage = share of "forbidden" reflections that are still significant (I/σ ≥ 3) — 0% is the clean case.
          <strong>d-glide</strong> = quarter translation (the condition is on a sum modulo 4, e.g. k+l = 4n). <strong>hhl</strong> is the plane through the diagonal (k = h): there the glide/mirror planes
          of tetragonal, trigonal/hexagonal and cubic groups show up; in other crystal systems that row only repeats the centering.
          A zonal condition consistent with the lattice centering alone (e.g. an apparent glide in a zone that intersects the
          centering condition) doesn't imply an independent glide plane — check the centering result first.
          <strong>"Undetermined"</strong> means the test itself isn't trustworthy either way — hover a cell for the specific reason.
          It happens for one of two causes: too few reflections in the allowed (reference) class to say anything, or the allowed
          class itself is mostly weak too (common for a condition that only really shows up at higher resolution) — in both cases
          the percentage next to it isn't meaningful evidence, whether it's 0% or high.
          A condition can also read <strong>"satisfied"</strong> despite a non-zero percentage: if the <em>median</em> I/σ of the forbidden class is below 3
          (hover the cell), the significant ones are just the upper tail of weak residual intensity. The plots below grey out cells forbidden by such accepted conditions.</small></p></details>`;

        document.getElementById('absenceVizSection').style.display = 'block';
        renderAbsenceVisuals(
          document.getElementById('zonePlotsContainer'),
          document.getElementById('axialStripsContainer'),
          lastReflections,
          scan,
          readCellOrNull()
        );
        resolve(true);
        }, 10);
      });
    }


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
        renderCrystalData();
        return;
      }
      renderCrystalData();
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
    // rather than requiring a manual step. Also
    // wired to an explicit "Apply formula" button/click, since the results
    // (Wilson plot, .ins) live on other tabs and a silent background
    // recompute is easy to miss — this button gives a deliberate action
    // plus visible confirmation right next to the field.
    function applyFormulaChange() {
      let did = [];
      if (lastReflections && lastMatchResults) { runEStatistics(true); did.push('E² statistics'); }
      if (lastMatchResults) { generateIns(); did.push('the SHELX .ins'); } // formula/Z feeds directly into the .ins too, regardless of E-stat outcome
      if (lastMatchResults) renderTwinCheck(); // the expected Bijvoet signal depends on the composition
      formulaUpdateNoteEl.innerHTML = did.length > 0
        ? `<p style="color:#080"><small>✓ Recalculated ${did.join(' and ')} with the current formula — check the corresponding tab.</small></p>`
        : `<p style="color:#b90"><small>No space-group match yet to recalculate — load reflection data first.</small></p>`;
    }
    formulaInput.addEventListener('change', applyFormulaChange);
    applyFormulaBtn.addEventListener('click', applyFormulaChange);
    formulaInput.addEventListener('input', () => { formulaManuallyEdited = true; renderCrystalData(); });

    // ---- Wavelength, Z, crystal size → "Contents & crystal data" ----------
    // λ and Z used to live only in lastCellExtras (whatever the CIF/.ins/.p4p
    // said); now they are visible, editable fields. A freshly loaded file that
    // carries a value sets the field (like the formula field); a later file
    // without one leaves it alone. With no λ from any file, Mo Kα is assumed
    // and marked as such — the .ins tab says so too.
    const lambdaSelect = document.getElementById('lambdaSelect');
    const lambdaCustom = document.getElementById('lambdaCustom');
    const lambdaSourceEl = document.getElementById('lambdaSource');
    const zInput = document.getElementById('zInput');
    const zSourceEl = document.getElementById('zSource');
    const sizeInputs = ['sizeD1', 'sizeD2', 'sizeD3'].map(id => document.getElementById(id));
    const sizeSourceEl = document.getElementById('sizeSource');
    // Crystal size: CIF _exptl_crystal_size_max/mid/min or a SHELX SIZE line
    // fill the three fields; editing any of them makes it "set by hand".
    function setSizeFromFile(size, label) {
      if (!size) return false;
      const sorted = [...size].sort((x, y) => y - x); // largest first, as in the CIF
      sizeInputs.forEach((inp, i) => { inp.value = sorted[i]; });
      sizeSourceEl.textContent = `(from ${label})`;
      return true;
    }
    function resetSize() {
      sizeInputs.forEach(inp => { inp.value = ''; });
      sizeSourceEl.textContent = '';
    }
    const crystalDataResultEl = document.getElementById('crystalDataResult');
    const zPrimeBoxEl = document.getElementById('zPrimeBox');
    let lambdaState = { source: 'assumed', label: null }; // source: 'file' | 'manual' | 'assumed'
    let zState = { source: null, label: null };            // source: 'file' | 'manual' | null (empty → estimated)
    // Same idea as formulaManuallyEdited: a Z typed by hand means formula × Z
    // should win over an exact SFAC/UNIT composition from a loaded .ins.
    let zManuallyEdited = false;

    function renderLambdaSource() {
      lambdaCustom.style.display = lambdaSelect.value === 'other' ? '' : 'none';
      const s = lambdaState;
      lambdaSourceEl.className = 'fieldSource' + (s.source === 'assumed' ? ' assumed' : '');
      lambdaSourceEl.textContent = s.source === 'file' ? `(from ${s.label}${s.exactFor && s.exactFor !== 'other' && Math.abs(s.exact - ANODE_WAVELENGTHS[s.exactFor]) > 1e-6 ? `: ${s.exact} Å` : ''})`
        : s.source === 'manual' ? '(set by hand)'
        : 'assumed — no wavelength in the loaded files, please check';
    }
    function setWavelengthFromFile(lambda, label) {
      if (!Number.isFinite(lambda) || lambda <= 0) return false;
      const anode = anodeForWavelength(lambda);
      lambdaSelect.value = anode || 'other';
      lambdaCustom.value = anode ? '' : lambda;
      // keep the file's exact value (e.g. 1.54178 from an older .ins) as long as
      // the matching anode stays selected — the .ins then writes it unchanged
      lambdaState = { source: 'file', label, exact: lambda, exactFor: lambdaSelect.value };
      renderLambdaSource();
      return true;
    }
    function resetWavelength() {
      lambdaSelect.value = 'Mo';
      lambdaCustom.value = '';
      lambdaState = { source: 'assumed', label: null };
      renderLambdaSource();
    }
    function currentWavelength() {
      const keepExact = lambdaState.source === 'file' && lambdaState.exactFor === lambdaSelect.value && lambdaSelect.value !== 'other';
      const v = keepExact ? lambdaState.exact
        : lambdaSelect.value === 'other' ? parseFloat(lambdaCustom.value) : ANODE_WAVELENGTHS[lambdaSelect.value];
      return { value: Number.isFinite(v) && v > 0 ? v : null, source: lambdaState.source };
    }

    function renderZSource() {
      zSourceEl.textContent = zState.source === 'file' ? `(from ${zState.label})`
        : zState.source === 'manual' ? '(set by hand)' : '';
    }
    function setZFromFile(z, label) {
      if (!(z > 0)) return false;
      zInput.value = z;
      zState = { source: 'file', label };
      zManuallyEdited = false;
      renderZSource();
      return true;
    }
    function resetZ() {
      zInput.value = '';
      zState = { source: null, label: null };
      zManuallyEdited = false;
      renderZSource();
    }
    function currentCellOrNull() {
      const v = {
        a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
      };
      return Object.values(v).every(Number.isFinite) ? v : null;
    }
    // Z from the field, otherwise estimated from V / (18 Å³ × non-H atoms).
    function currentZ(formula) {
      const typed = parseFloat(zInput.value);
      if (typed > 0) return { z: typed, estimated: false };
      const cell = currentCellOrNull();
      if (!cell || !formula || formula.length === 0) return { z: null, estimated: false };
      const det = estimateZDetails(cellVolume(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma), formula);
      return det ? { z: det.z, estimated: true } : { z: null, estimated: false };
    }
    // Exact SFAC/UNIT from a .ins wins until formula or Z are edited by hand.
    function useExactComposition() {
      return !!(lastInsComposition && lastInsComposition.length > 0 && !formulaManuallyEdited && !zManuallyEdited);
    }

    function fmtNum(x, digits) { return Number.isFinite(x) ? x.toFixed(digits) : '—'; }
    function fmtFpp(x) { return x >= 1 ? x.toFixed(2) : x >= 0.01 ? x.toFixed(3) : x.toPrecision(2); }

    // Z′ for the selected space-group candidate, or null.
    function selectedZPrime(z) {
      if (!z || !lastMatchResults || !lastMatchResults[selectedCandidateIndex]) return null;
      const entry = lastMatchResults[selectedCandidateIndex].entry;
      const info = zPrimeInfo(z, entry.hall);
      return info ? { ...info, entry } : null;
    }
    function zPrimeText(zp) {
      const val = Number.isInteger(zp.zPrime) ? zp.zPrime : (Math.abs(zp.zPrime * 4 - Math.round(zp.zPrime * 4)) < 1e-9 ? zp.zPrime.toString() : zp.zPrime.toFixed(3));
      return `Z′ = ${val} in ${formatHM(zp.entry.hm)} (${zp.mult} general positions)`;
    }

    function renderCrystalData() {
      const cell = currentCellOrNull();
      const formula = parseFormulaString(formulaInput.value);
      const zi = currentZ(formula);
      const lam = currentWavelength();
      zInput.placeholder = zi.estimated ? `≈ ${zi.z}` : '';
      if (!cell) {
        crystalDataResultEl.innerHTML = `<p><small>Enter or load the cell to see Z estimate, density, F(000) and μ.</small></p>`;
        zPrimeBoxEl.innerHTML = '';
        return;
      }
      const volume = cellVolume(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
      const parts = [];

      // Z estimate (always shown when a formula is there, also as a check of a given Z)
      const det = formula.length > 0 ? estimateZDetails(volume, formula) : null;
      if (det) {
        const off = Math.abs(det.zRaw - det.z) > 0.25;
        const given = parseFloat(zInput.value);
        const disagree = given > 0 && Math.abs(given - det.zRaw) / det.zRaw > 0.3;
        parts.push(`<p><small>Z estimate: V / 18 Å³ ≈ ${det.atomsPerCell.toFixed(1)} non-H atoms per cell; the formula has ${det.nonHPerUnit} →
          Z ≈ ${det.zRaw.toFixed(2)} → <strong>${det.z}</strong>${off ? ' <span style="color:var(--warning)">(far from an integer — check the formula, or solvent/Z′)</span>' : ''}.
          ${disagree ? `<span style="color:var(--warning)">The given Z = ${given} differs from this by more than 30 %.</span>` : ''}
          Rule of thumb for organic/organometallic compounds; less reliable for inorganic or very dense structures.</small></p>`);
      }

      // Composition for density/F(000)/μ
      let composition = null, compNote = '';
      if (useExactComposition()) {
        composition = lastInsComposition;
        compNote = 'exact SFAC/UNIT from the .ins';
      } else if (formula.length > 0 && zi.z) {
        composition = formula.map(f => ({ element: f.element, count: f.count * zi.z }));
        compNote = `formula × Z = ${zi.z}${zi.estimated ? ' (estimated)' : ''}`;
      }
      if (!composition) {
        parts.push(`<p><small>Enter the chemical formula (and Z, or let it be estimated) for density, F(000) and μ.</small></p>`);
        crystalDataResultEl.innerHTML = parts.join('');
        renderZPrimeBox(zi);
        return;
      }

      const size = sizeInputs.map(i => parseFloat(i.value));
      const cd = computeCrystalData({ volume, composition, z: zi.z, formula, wavelength: lam.value, size });
      if (cd && cd.unknown) {
        parts.push(`<p style="color:var(--error)"><small>Unknown element symbol(s): ${cd.unknown.join(', ')} — no crystal data calculated.</small></p>`);
        crystalDataResultEl.innerHTML = parts.join('');
        renderZPrimeBox(zi);
        return;
      }

      const rows = [];
      rows.push(['Cell volume', `${volume.toFixed(1)} Å³`]);
      rows.push(['Composition used', `${composition.map(c => c.element + (Number.isInteger(c.count) ? c.count : c.count.toFixed(2))).join(' ')} <small>(${compNote})</small>`]);
      if (cd.formulaWeight) rows.push(['Formula weight', `${cd.formulaWeight.toFixed(2)} g/mol`]);
      rows.push(['Density (calc.)', `${cd.density.toFixed(3)} g/cm³`]);
      rows.push(['F(000)', `${Math.round(cd.f000)}`]);
      const zp = selectedZPrime(zi.z);
      if (zp) rows.push(['Z′', zPrimeText(zp).replace(/^Z′ = /, '')]);
      if (cd.mu !== null) {
        rows.push([`μ (${cd.anode} Kα)`, `${cd.mu.toFixed(3)} mm⁻¹ <small>(NIST FFAST; SHELXL uses International Tables Vol. C and may differ by ~1 %)</small>`]);
      } else {
        rows.push(['μ', `<small>${lam.value ? `only for Mo, Cu and Ag Kα (λ = ${lam.value} Å has no tabulated values here)` : 'no wavelength'}</small>`]);
      }
      if (cd.absorption) {
        const a = cd.absorption;
        const advice = {
          'negligible': 'negligible — an absorption correction changes little',
          'moderate': 'moderate — a multi-scan correction (SADABS/SCALE3 ABSPACK) is usually enough',
          'strong': 'strong — multi-scan needed; consider a numerical (face-indexed) correction',
          'very strong': 'very strong — a numerical (face-indexed) correction is recommended, or a smaller crystal / shorter wavelength',
        }[a.level];
        rows.push(['μ·r (min / max)', `${a.muRMin.toFixed(2)} / ${a.muRMax.toFixed(2)} — <span class="absLevel ${a.level.replace(' ', '-')}">${advice}</span>`]);
        rows.push(['T<sub>min</sub> / T<sub>max</sub> (expected)', `${a.tMin.toFixed(3)} / ${a.tMax.toFixed(3)} <small>(exp(−μ·d<sub>max</sub>) / exp(−μ·d<sub>min</sub>) — the outer limits, as checkCIF's T<sub>min</sub>′ / T<sub>max</sub>; a real correction usually gives a somewhat higher T<sub>min</sub>)</small>`]);
      }
      const fpps = cd.perElement.filter(p => p.fpp !== null).sort((x, y) => y.fpp - x.fpp);
      if (fpps.length > 0) {
        rows.push([`f″ (${cd.anode} Kα)`, fpps.map(p => `${p.element} ${fmtFpp(p.fpp)}`).join(' · ')
          + ` <small>— anomalous signal for the absolute structure${fpps[0].fpp < 0.1 ? '; very weak here (light atoms only) — the Flack parameter will be imprecise' + (cd.anode !== 'Cu' ? ', Cu Kα would help' : '') : ''}</small>`]);
      }
      parts.unshift(`<table class="crystalDataTable">${rows.map(([k, v]) => `<tr><td>${k}</td><td>${v}</td></tr>`).join('')}</table>`);
      parts.push(`<p><small>μ/ρ and f″: NIST FFAST (Chantler et al., 2005), interpolated at the anode wavelength.</small></p>`);
      crystalDataResultEl.innerHTML = parts.join('');
      renderZPrimeBox(zi);
    }

    // Short Z′ line under the candidate table (Space Group tab).
    function renderZPrimeBox(zi) {
      const zp = selectedZPrime(zi.z);
      if (!zp) { zPrimeBoxEl.innerHTML = ''; return; }
      zPrimeBoxEl.innerHTML = `<p><strong>${zPrimeText(zp)}</strong> <small>for the selected row, with Z = ${zi.z}${zi.estimated ? ' (estimated — set Z on the File &amp; Cell tab)' : ''}.</small>
        ${zp.note ? `<br><span style="color:var(--warning)"><small>${zp.note}</small></span>` : ''}</p>`;
    }

    lambdaSelect.addEventListener('change', () => {
      lambdaState = { source: 'manual', label: null };
      renderLambdaSource();
      renderCrystalData();
      if (lastMatchResults) { renderTwinCheck(); generateIns(); } // f″ for the Bijvoet check depends on the anode
    });
    lambdaCustom.addEventListener('input', () => { lambdaState = { source: 'manual', label: null }; renderLambdaSource(); renderCrystalData(); });
    lambdaCustom.addEventListener('change', () => { if (lastMatchResults) { renderTwinCheck(); generateIns(); } });
    zInput.addEventListener('input', () => {
      zState = { source: zInput.value.trim() === '' ? null : 'manual', label: null };
      zManuallyEdited = zInput.value.trim() !== '';
      renderZSource();
      renderCrystalData();
    });
    zInput.addEventListener('change', applyFormulaChange); // Z feeds E² composition, .ins UNIT and the Bijvoet expectation
    sizeInputs.forEach(i => i.addEventListener('input', () => {
      sizeSourceEl.textContent = sizeInputs.some(x => x.value.trim() !== '') ? '(set by hand)' : '';
      renderCrystalData();
    }));
    sizeInputs.forEach(i => i.addEventListener('change', () => { if (lastMatchResults) generateIns(); })); // SIZE line in the .ins
    renderLambdaSource();
    renderCrystalData();
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
    let lastFpp = null; // f'' per element from the CIF atom-type loop (anomalous-check.js)
    function applyCellFromCIF(cell, text) {
      const fpp = extractDispersionFromCIF(text);
      if (fpp) lastFpp = fpp; // a reflection-only CIF without the loop keeps the earlier values
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
      if (cell.reportedHallSymbol !== null) extras.push(`Hall: ${cell.reportedHallSymbol}`);
      cellExtraEl.innerHTML = extras.length > 0 ? `<p>${extras.join(' · ')}</p>` : '';

      const formula = parseChemicalFormula(text);
      formulaInput.value = formula ? formula.map(f => f.element + f.count).join(' ') : '';
      formulaManuallyEdited = false; // fresh load
      if (!setZFromFile(cell.z, 'CIF')) resetZ();
      if (!setWavelengthFromFile(cell.wavelength, 'CIF')) resetWavelength();
      if (!setSizeFromFile(cell.size, 'CIF')) resetSize();
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
      if (!setZFromFile(insData.z, '.ins')) resetZ();
      if (!setWavelengthFromFile(insData.wavelength, '.ins')) resetWavelength();
      if (!setSizeFromFile(insData.size, 'SIZE in the .ins/.res')) resetSize();

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
      setWavelengthFromFile(p4pData.wavelength, '.p4p');

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
        noteLoaded('cell', file.name);
        return;
      }

      // A .ins carries cell + exact whole-cell composition (SFAC/UNIT), but
      // no reflections at all — handle it separately rather than forcing
      // it through the reflection-parsing path, which would just report
      // "no reflection data found".
      if (looksLikeSHELXIns(text)) {
        handleInsFile(text);
        noteLoaded('cell', file.name);
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
            noteLoaded('cell', file.name);
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
      noteLoaded('reflections', file.name, `(${reflections.length} refl.)`);
      determinedLaueClass = null;
      determinedLaueClassSource = null;
      laueResultEl.innerHTML = '';
      absencesResultEl.innerHTML = '';
      spaceGroupResultEl.innerHTML = '';
      eStatResultEl.innerHTML = '';
      const stats = summarizeHKL(reflections);

      sourceInfoEl.innerHTML = `<p><strong>Source:</strong> ${SOURCE_LABELS[source] || source}</p>`;

      // HKL export (Tab 7): only offered when the reflections came out of a
      // CIF — a raw SHELX .hkl is already in the right format. Written
      // untransformed; the embedded block is kept verbatim, a _refln loop is
      // re-written in SHELX format (3I4,2F8.2,I4).
      if (source === 'raw_hkl') {
        lastHklExport = null;
      } else {
        lastHklExport = {
          source,
          count: reflections.length,
          text: source === 'embedded_shelx_hkl_file' && result.rawHklText
            ? result.rawHklText.replace(/^\s*\n/, '').replace(/\s*$/, '\n')
            : writeSHELXHKL(reflections),
        };
      }
      renderHklExport();

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
        // Exception: a reflection-only file (e.g. a COD-style .hkl in CIF
        // format) usually repeats the cell but carries no formula/Z. If
        // the very same cell is already loaded together with a formula
        // (from the real CIF/.ins dropped before), keep that richer data
        // instead of wiping the formula.
        const sameCellAsLoaded = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
          .every(k => Math.abs(parseFloat(cellFields[k].value) - cell[k]) < 1e-3);
        const fileHasNoFormula = !parseChemicalFormula(text);
        if (fileHasNoFormula && sameCellAsLoaded && (formulaInput.value.trim() !== '' || lastInsComposition)) {
          cellNoteEl.innerHTML = `<p><em>Same cell as already loaded — keeping formula/Z from the earlier file.</em></p>`;
        } else {
          applyCellFromCIF(cell, text);
          noteLoaded('cell', file.name);
        }
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
        resetZ();
        resetWavelength();
        resetSize();
        lastFpp = null; // f″ of an earlier CIF belong to another structure
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
    async function handleFiles(fileList) {
      const list = Array.from(fileList || []).filter(Boolean);
      if (list.length === 0) return;
      // Files that carry reflections go first, cell/formula-only files
      // (.ins/.p4p/CIF without reflection loop) last, so the latter can
      // add their formula/composition on top instead of being overwritten.
      const texts = await Promise.all(list.map(f => f.text()));
      const rank = (t) => (looksLikeP4P(t) || looksLikeSHELXIns(t)) ? 1 : (parseReflectionInput(t).reflections.length > 0 ? 0 : 1);
      const ordered = list.map((f, i) => ({ f, r: rank(texts[i]) })).sort((a, b) => a.r - b.r).map(x => x.f);
      batchDepth++;
      try {
        for (const f of ordered) await handleFile(f);
      } finally {
        batchDepth--;
      }
      if (batchDepth === 0 && pendingRun) {
        pendingRun = false;
        await autoRunPipeline();
      }
    }

    async function autoRunPipeline() {
      // While several files are being read in one go, run the chain only
      // once at the end (otherwise it would run — and jump — per file).
      if (batchDepth > 0) { pendingRun = true; return; }
      if (!lastReflections) return;
      const values = {
        a: parseFloat(cellFields.a.value), b: parseFloat(cellFields.b.value), c: parseFloat(cellFields.c.value),
        alpha: parseFloat(cellFields.alpha.value), beta: parseFloat(cellFields.beta.value), gamma: parseFloat(cellFields.gamma.value),
      };
      if (!Object.values(values).every(Number.isFinite)) return;

      await runLaueDetermination();
      renderLaueTab();
      await runAbsencesScan();
      const matched = await runSpaceGroupMatch();
      if (matched) {
        // Wilson-plot/E² statistics only if a composition is already known
        // (exact .ins SFAC/UNIT, or a filled-in formula + Z) — otherwise
        // leave it for a manual click once the formula is entered.
        await runEStatistics(true);
      }
      // Jump to the candidate table only while the person is still on the
      // upload tab. Anyone who is already further along (e.g. adding the
      // CIF afterwards from another tab) stays exactly where they are.
      const activePanel = document.querySelector('.tabPanel.active');
      const compositionKnown = formulaInput.value.trim() !== '' || !!lastInsComposition;
      if (matched && typeof window.switchTab === 'function' && activePanel && activePanel.id === 'tab-file') {
        if (compositionKnown) {
          // Everything needed is in: show the candidates.
          awaitingJump = false;
          window.switchTab('tab-spacegroup');
        } else {
          // Reflections + cell are there, but formula/Z (usually from the
          // CIF) is not — likely another file is still coming. Stay on the
          // upload tab instead of jumping away from it.
          awaitingJump = true;
        }
        renderLoadedFiles();
      } else if (matched && awaitingJump) {
        awaitingJump = false;
        renderLoadedFiles();
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
        const mergingResult = determineLaueClass(lastReflections, classification.laueClasses, { monoclinicAxis });
        const { jumpDetectedAt, evaluations, redundancyInformative } = mergingResult;
        // A CIF space-group assignment is an independent symmetry source.
        // Its point group, completed by Friedel inversion, fixes the Laue
        // class; merged deposited Fo/Fc datasets often have no equivalent
        // repeats, so Rint cannot recover it. Prefer it when consistent with
        // the observed cell metric; otherwise retain the data-driven result.
        const cifLaueClass = laueClassFromCIFMetadata(lastCellExtras, classification.laueClasses, monoclinicAxis);
        const chosen = cifLaueClass || mergingResult.chosen;

        // With no redundant (duplicate) observations under ANY candidate,
        // Rint is undefined everywhere and this method has no signal at
        // all — "chosen" would otherwise silently default to the lowest
        // candidate, which looks like a confident answer but isn't one.
        // Don't treat that as determined; leave it unset so the
        // space-group matching step keeps warning the user honestly
        // instead of quietly trusting it.
        determinedLaueClass = cifLaueClass || (redundancyInformative ? chosen : null);
        determinedLaueClassSource = cifLaueClass ? 'cif' : (redundancyInformative ? 'merging' : null);

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

        const verdictBlock = cifLaueClass
          ? `<p><strong>Laue class from CIF space group: ${chosen}</strong> — ${lastCellExtras.reportedSpaceGroup || `space-group #${lastCellExtras.reportedSpaceGroupNumber}`} (${lastCellExtras.reportedHallSymbol || 'Hall symbol unavailable'}). The R<sub>int</sub> table is shown as an independent check; merged/refinement HKL may not contain redundancy to redetermine the class.</p>
             ${jumpNote}`
          : redundancyInformative
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

    const spaceGroupResultEl = document.getElementById('spaceGroupResult');

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
      renderCrystalData(); // Z′ depends on the selected candidate
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
      const lam = currentWavelength();
      const wavelength = lam.value;
      const zi = currentZ(formula);
      const z = zi.z;
      const zWasEstimated = zi.estimated;

      // If the selected candidate isn't the one whose own systematic
      // absences the as-entered cell/hkl natively satisfy (see
      // describeSetting, same logic), the cell needs transforming into the
      // selected setting's frame, and the .hkl file needs the
      // corresponding HKLF matrix so SHELXL re-indexes it to match.
      const sameNumberResults = lastMatchResults.filter(r => r.entry.number === top.number);
      const nativeEntry = sameNumberResults.reduce((best, r) => (r.score > best.score ? r : best), sameNumberResults[0]).entry;
      let settingTransform = null;
      if (nativeEntry.hall !== top.hall) {
        const matrix = deriveTransformMatrix(nativeEntry.hall, top.hall, values);
        if (matrix && !isIdentity3(matrix)) settingTransform = { matrix, nativeHM: nativeEntry.hm, isPurePermutation: isPurePermutation(matrix), obliqueNote: obliquenessNote(top, matrix, values) };
      }

      const spaceGroupTag = top.hm.trim().replace(/\s+/g, ' '); // H-M with its spaces (P 1 21/a 1), as in the Hall table
      const ins = generateInsFile({
        title: lastLoadedFileBaseName ? `${lastLoadedFileBaseName} in ${spaceGroupTag}` : spaceGroupTag,
        cell: values,
        wavelength,
        hallSymbol: top.hall,
        z,
        formula: formula.length > 0 ? formula : null,
        composition: useExactComposition() ? lastInsComposition : null, // exact whole-cell counts, if a .ins/.res was loaded — takes priority over formula×Z UNLESS the person has since typed their own formula or Z
        esd: lastCellEsd, // cell esds, if the source file provided them (CIF parens / .ins ZERR / .p4p CELLSD)
        settingTransform,
        twin: twinForIns(),
        size: sizeInputs.map(i => parseFloat(i.value)),
      });

      const hasComposition = useExactComposition() || formula.length > 0;
      const missingParts = [!wavelength && 'wavelength', !hasComposition && 'chemical formula'].filter(Boolean);
      const lambdaNote = wavelength && lam.source === 'assumed'
        ? `<p style="color:#b90"><small>λ = ${wavelength} Å (Mo Kα) is assumed — none of the loaded files gives a wavelength. Set it on the File &amp; Cell tab if your data are not Mo.</small></p>`
        : '';
      const missingNote = missingParts.length > 0
        ? `<p style="color:#b00"><small>Missing: ${missingParts.join(', ')} — placeholder(s) inserted, fill in manually before use.</small></p>`
        : '';
      const zNote = zWasEstimated
        ? `<p style="color:#b90"><small>Z = ${z} was estimated from cell volume ÷ (~18 Å³ × non-H atoms) — a starting guess, not a measured value. Sanity-check it (Z is usually 1, 2, 3, 4, 6, or 8) and enter it on the File &amp; Cell tab.</small></p>`
        : '';
      // If the SELECTED candidate is tied with another that disagrees on
      // centrosymmetry, absences alone couldn't break that tie on their own.
      // Only worth flagging for an unreviewed default pick or to confirm an
      // E²-statistics-backed choice — a deliberate manual selection doesn't
      // need to be second-guessed here.
      const tieNote = tiedOpposite && selectionSource === 'cif'
        ? `<p style="color:#069"><strong>${formatHM(top.hm)} is the space group the CIF reports.</strong> The absences cannot distinguish it from other top-scoring candidates
           (see the Space Group tab); the data agree with it as well as with any of them.</p>`
        : tiedOpposite && selectionSource !== 'manual'
        ? (selectionSource === 'auto-e-stat'
            ? `<p style="color:#080"><strong>✓ ${formatHM(top.hm)} was picked over the tied ${formatHM(tiedOpposite.entry.hm)} (${tiedOpposite.entry.hall}) based on the E² statistics (tab 5)</strong> —
               absences alone couldn't distinguish them (both ${(topScore * 100).toFixed(0)}%), but the intensity-distribution check can, and did.</p>`
            : `<p style="color:#b00"><strong>⚠ This is tied at ${(topScore * 100).toFixed(0)}% with ${formatHM(tiedOpposite.entry.hm)} (${tiedOpposite.entry.hall}), which is ${topIsCentro ? 'acentric' : 'centrosymmetric'}
               where this one is ${topIsCentro ? 'centrosymmetric' : 'acentric'}.</strong> Systematic absences alone cannot break a centrosymmetric/acentric tie — Friedel's law makes both look
               identical in a diffraction pattern regardless of which one the crystal actually is. This is just the top-scoring row, not a real determination between the two — run
               <strong>Compute E² statistics</strong> on the Space Group tab, then select whichever row it points to.</p>`)
        : '';

      // Twinning found on the Twinning tab but not (all) written to this .ins
      const findings = twinFindings() || [];
      const written = twinForIns();
      const notWritten = findings.filter(x => x.kind === 'law' ? !(written && written.matrix.flat().join(' ') !== '-1 0 0 0 -1 0 0 0 -1')
        : !(written && (written.n === -4 || written.matrix.flat().join(' ') === '-1 0 0 0 -1 0 0 0 -1'))).filter(x => x.kind === 'law' || !hallSymbolIsCentrosymmetric(top.hall));
      const twinInsNote = notWritten.length
        ? `<p style="color:#b60"><strong>⚠ Twinning found but not in this .ins:</strong> ${notWritten.map(x => x.text).join(' · ')} — tick it on the <a onclick="window.switchTab('tab-twin')" style="cursor:pointer">Twinning tab</a> to add TWIN/BASF.</p>`
        : '';

      const settingNote = settingTransform
        ? `<p style="color:#069"><small>Cell transformed from as-measured (${formatHM(settingTransform.nativeHM)}) into this setting; the .hkl file itself is untouched — SHELXL re-indexes it via the HKLF matrix.</small>${settingTransform.obliqueNote || ''}</p>`
        : '';

      insResultEl.innerHTML = `
        <h2>SHELX .ins (space group ${top.number}${top.qualifier ? ':' + top.qualifier : ''}, ${formatHM(top.hm)})</h2>
        <p><small>Kept in sync automatically with whichever candidate is selected on the Space Group tab — edit the cell/formula fields or the selection there and this updates on its own.</small></p>
        ${tieNote}
        ${twinInsNote}
        ${settingNote}
        ${missingNote}
        ${lambdaNote}
        ${zNote}
        <pre style="background:#f7f7f7; padding:0.8rem; border:1px solid #ddd; overflow-x:auto;">${ins}</pre>`;
      lastGeneratedIns = ins;
      insActionsEl.style.display = 'flex';
      insCopyNoteEl.textContent = '';
    }

    // Tab 7: "Save as .hkl" block, visible only while lastHklExport is set.
    const hklSaveSectionEl = document.getElementById('hklSaveSection');
    const hklSaveInfoEl = document.getElementById('hklSaveInfo');
    const hklSaveBtn = document.getElementById('hklSaveBtn');
    function renderHklExport() {
      if (!lastHklExport) { hklSaveSectionEl.style.display = 'none'; return; }
      const origin = lastHklExport.source === 'embedded_shelx_hkl_file'
        ? 'the HKL block embedded in the CIF (<code>_shelx_hkl_file</code>), saved verbatim'
        : 'the CIF reflection loop (<code>_refln</code>), written in SHELX format';
      hklSaveInfoEl.innerHTML = `<small>${lastHklExport.count} reflections from ${origin}. Indices and intensities are written exactly as read —
        no re-indexing or merging; a setting change is handled by the HKLF line in the .ins. Saved under the same base name as the .ins.</small>`;
      hklSaveSectionEl.style.display = 'block';
    }
    hklSaveBtn.addEventListener('click', () => {
      if (!lastHklExport) return;
      const blob = new Blob([lastHklExport.text], { type: 'text/plain' });
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = `${lastLoadedFileBaseName || 'structure'}.hkl`;
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    });

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
        const laueClassNote = determinedLaueClassSource === 'cif'
          ? `from CIF space-group symmetry`
          : determinedLaueClassSource === 'merging'
            ? `confirmed via merging statistics`
            : `assumed — highest-symmetry candidate for this cell; no source confirmed a Laue class`;

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
        // Candidates the absences cannot tell apart at all (identical predicted
        // conditions = same extinction symbol, e.g. I 4 2 2 / I 4 m m / I -4 m 2 /
        // I -4 2 m / I 4/m m m for "I - - -"). If the CIF names one of the
        // top-scoring rows, preselect it: the data don't contradict it, and
        // nothing in them could prefer another member of the tie.
        const topSignature = conditionSignature(results[0].entry.hall);
        const sameAbsences = results.filter(r => r.score === results[0].score && conditionSignature(r.entry.hall) === topSignature);
        const cifIndex = cifReportedCandidateIndex(results);
        if (cifIndex >= 0 && results[cifIndex].score === results[0].score) { // also when it is already the top row: E² must not overturn it
          selectedCandidateIndex = cifIndex;
          selectionSource = 'cif';
        }
        generateIns(); // pre-populate the .ins tab with the top candidate right away
        renderTwinCheck(); // the Bijvoet (inversion-twin) check needs the candidates

        const topTen = results.slice(0, 10);
        // The candidate list is ranked across ALL space-group numbers, so a
        // crowd of similarly-scoring pseudo-symmetric candidates from OTHER
        // numbers can push a number's own qualifier-null (ITA standard)
        // sibling below the visible cutoff — even when it's clearly the
        // right partner for the top hit (e.g. #61 Pbca scoring below several
        // unrelated #29/#54/#57/#60 candidates, while its own Pcab sibling
        // is the #1 match). Guarantee that sibling stays reachable: if the
        // top candidate isn't itself already standard and its standard
        // sibling didn't make the top 10 on score alone, append it.
        const topEntry = results[0].entry;
        // "Standard" here means whatever this tool actually knows how to
        // transform TO: qualifier null/"abc" (orthorhombic-style) or "b1"
        // (monoclinic b-unique choice 1) — NOT just null/"abc". Monoclinic
        // space groups never even HAVE a null-qualifier entry (checked:
        // #7/9/13/14/15 only ever use a1../b1../c1..), so the earlier
        // null/abc-only check would silently never find a monoclinic
        // sibling at all, regardless of its rank.
        const isStandardQualifier = q => q == null || q === 'abc' || q === 'b1';
        const topIsStandard = isStandardQualifier(topEntry.qualifier);
        let extraRow = null;
        if (!topIsStandard) {
          const siblingIndex = results.findIndex(r => r.entry.number === topEntry.number && isStandardQualifier(r.entry.qualifier));
          if (siblingIndex >= 10) extraRow = { index: siblingIndex, result: results[siblingIndex] };
        }

        const rowHtml = (r, i) => `
          <tr>
            <td><input type="radio" name="candidateSelect" value="${i}" ${i === selectedCandidateIndex ? 'checked' : ''}></td>
            <td>${r.entry.number}${r.entry.qualifier ? ':' + r.entry.qualifier : ''}</td>
            <td>${formatHM(r.entry.hm)}</td>
            <td>${r.entry.hall}</td>
            <td>${(r.score * 100).toFixed(1)}%</td>
            <td>${describeSetting(r.entry)}</td>
          </tr>`;
        const rows = topTen.map((r, i) => rowHtml(r, i)).join('')
          + (extraRow ? `<tr><td colspan="6"><small>↓ standard setting of the top match, shown despite its rank:</small></td></tr>${rowHtml(extraRow.result, extraRow.index)}` : '');

        const topScore = results[0].score;
        const topIsCentro = hallSymbolIsCentrosymmetric(results[0].entry.hall);
        const tiedOpposite = results.find(r =>
          r !== results[0] && r.score === topScore && hallSymbolIsCentrosymmetric(r.entry.hall) !== topIsCentro);
        const sameAbsencesHasBoth = sameAbsences.some(r => hallSymbolIsCentrosymmetric(r.entry.hall)) && sameAbsences.some(r => !hallSymbolIsCentrosymmetric(r.entry.hall));
        const tieWarning = tiedOpposite && !(sameAbsences.length > 2 && sameAbsencesHasBoth)
          ? `<p style="color:#b00"><strong>⚠ The top ${(topScore * 100).toFixed(0)}% match is tied between an acentric and a centrosymmetric option</strong>
             (${formatHM(results[0].entry.hm)} vs. ${formatHM(tiedOpposite.entry.hm)}) — systematic absences alone can never break that particular tie, by Friedel's law.
             The <strong>E² Statistics</strong> tab shows which one the intensity distribution favours — pick that one using the radio buttons — the .ins tab uses
             whichever row is selected, not just the top one.</p>`
          : '';

        // A check is a CLEAR contradiction if the scan reached a verdict (holds / violated)
        // and the candidate predicts the opposite. (Not "match < 20 %": a condition judged
        // "violated" at 67 % scores 0.33 for a candidate that predicts it, yet is clearly wrong.)
        const contradicts = c => c.verdict && c.expectHolds !== undefined &&
          ((c.verdict === 'holds' && !c.expectHolds) || (c.verdict === 'violated' && c.expectHolds));
        // Weight floor against noise: 8 reflections for zone checks. An axial row is
        // short by nature (b = 5.5 A: only 020/040/060 allowed, 7 positions in all), so
        // there the scan's own minimum of 3 allowed reflections (absences.js) counts —
        // otherwise a 2_1 along a short axis could never separate P 21/a from P 2/a.
        // A clear axial VIOLATION (strong forbidden reflections, see judgeCondition) needs
        // no allowed reference row: two significant forbidden reflections suffice.
        const isWeighty = c => (c.axial
          ? (c.allowedCount || 0) >= 3 || (c.verdict === 'violated' && (c.forbiddenCount || 0) >= 2)
          : c.weight >= 8);
        const decisiveAgainst = r => r.checks.filter(c => contradicts(c) && isWeighty(c)).filter(c => {
          const t = results[0].checks.find(x => x.label === c.label);
          return !t || !contradicts(t);
        });
        const describeContradiction = c => {
          const stricter = c.label.match(/^(\S+) period (\d+) \(stricter than the predicted (\d+)\)/);
          if (stricter) return `${stricter[1]}: the data obey period ${stricter[2]}, but this candidate only predicts period ${stricter[3]}`;
          const none = c.label.match(/^(\S+) period 1$/);
          if (none) return `${none[1]}: the data show absences along this axis (a screw-axis condition), but this candidate predicts none`
            + (c.forbiddenCount === 0 ? ` <small>(the forbidden reflections are not in the file at all — typical for deposited data with the absences removed, so this rests on the ${c.allowedCount} allowed ones)</small>` : '');
          return c.expectHolds ? `${c.label}: predicted by this candidate, but the data violate it`
                               : `${c.label}: the data obey it, but this candidate does not predict it`;
        };

        // Close runners-up that LOSE on a clear difference — otherwise their percentages
        // look like a tie (the shared zone checks involve thousands of reflections, a
        // decisive axial check only dozens).
        const closeRows = results.slice(1, 8).filter(r => topScore - r.score <= 0.03 && r.score !== topScore);
        const runnerUpNotes = closeRows
          .map(r => ({ r, bad: decisiveAgainst(r) })).filter(x => x.bad.length).slice(0, 5)
          .map(({ r, bad }) => `<li><strong>${formatHM(r.entry.hm)}</strong> (${(r.score * 100).toFixed(1)}%): ${bad.map(describeContradiction).join('; ')}</li>`);
        const runnerUpHtml = runnerUpNotes.length
          ? `<details class="explain" open><summary>Why close runners-up score lower</summary><p><small>The percentages above can look like a tie because the shared checks involve far more reflections than the decisive ones. These rows lose on a clear difference:</small></p><ul>${runnerUpNotes.join('')}</ul></details>`
          : '';

        // Unresolved: within 1.5 points of the top row and NO clear contradiction separates them
        // (an exact acentric/centrosymmetric tie is already covered by the E² warning above).
        const unresolved = results.slice(1, 6).filter(r =>
          topScore - r.score <= 0.015 &&
          !sameAbsences.includes(r) && // covered by the "identical absences" note below — not a data-quality question
          !(r.score === topScore && hallSymbolIsCentrosymmetric(r.entry.hall) !== topIsCentro) &&
          decisiveAgainst(r).length === 0);
        const cifRow = selectionSource === 'cif' ? results[selectedCandidateIndex] : null;
        const sameAbsencesHtml = sameAbsences.length > 1
          ? `<p style="color:#069"><strong>ℹ Identical systematic absences:</strong> ${sameAbsences.map(r => formatHM(r.entry.hm)).join(', ')}
              predict exactly the same absent reflections. No absence analysis can tell them apart (XPREP neither) — the equal percentages are not a data problem.
              ${sameAbsencesHasBoth ? 'The <strong>E² Statistics</strong> tab decides centrosymmetric vs. acentric; between the remaining ones' : 'Between them,'}
              only structure solution and refinement can decide (try each; check the result with ADDSYM in PLATON).
              ${cifRow ? `<br><strong>${formatHM(cifRow.entry.hm)} is preselected because the CIF reports it</strong> — the data agree with it equally well.` : ''}</p>`
          : (cifRow && selectedCandidateIndex > 0 ? `<p style="color:#069"><strong>${formatHM(cifRow.entry.hm)} is preselected because the CIF reports it</strong> — it scores as high as the top row.</p>` : '');
        const unresolvedHtml = unresolved.length
          ? `<p style="color:#b60"><strong>⚠ ${formatHM(results[0].entry.hm)} (${(topScore * 100).toFixed(1)}%) is not clearly separated from
              ${unresolved.map(r => `${formatHM(r.entry.hm)} (${(r.score * 100).toFixed(1)}%)`).join(', ')}</strong> — no check contradicts the runner-up
              <small>(within 1.5 points, and none of the differing conditions is clearly held or clearly violated).
              Typical causes: weak data, pseudo-symmetry (a weak pseudo-centering, a glide plane that is only partly obeyed), a condition this scan does not cover, or an
              enantiomorphic pair (P 4₁ / P 4₃ …) that absences can never tell apart. Look at the layer viewer (Laue tab) and the zone plots (Absences tab), and treat the top row
              as a proposal — the refinement (R1, Flack parameter, ADDSYM in PLATON) has the last word.</small></p>`
          : '';

        spaceGroupResultEl.innerHTML = `
          <h2>Space group candidates</h2>
          <p><small>Assumed Laue class: ${targetLaueClass} (${laueClassNote})</small></p>
          ${tieWarning}
          ${sameAbsencesHtml}
          ${unresolvedHtml}
          <table>
            <tr><th>Use</th><th>No.</th><th>H-M symbol</th><th>Hall symbol</th><th>Match</th><th>Setting</th></tr>
            ${rows}
          </table>
          ${runnerUpHtml}
          <p><small>Pick a different row any time before generating the .ins — it always uses whichever one is selected here, not just the top-scoring one.</small></p>
          <details class="explain"><summary>How the match score works</summary><p><small>Match % = share of calculated reflection conditions that agree with the observed scan, weighted by how many reflections
          actually informed each check — a solid signal from thousands of reflections outweighs a shaky one from a sparse axial row (a check
          where the data was inconclusive — too few reflections either way — doesn't count for or against at all). A less-than-100% match on an
          otherwise clear top candidate is often just the same small-sample noise the weighting already discounts, not a wrong space group.
          Acentric candidates are included alongside centrosymmetric ones — Friedel's law makes diffraction always look centrosymmetric
          regardless of the true crystal symmetry, so absences alone usually can't rule acentric options out; a close acentric/centrosymmetric
          tie here is exactly the situation E² statistics is for.</small></p></details>`;

        spaceGroupResultEl.querySelectorAll('input[name="candidateSelect"]').forEach(input => {
          input.addEventListener('change', (e) => {
            selectedCandidateIndex = parseInt(e.target.value, 10);
            selectionSource = 'manual';
            generateIns(); // rebuild the .ins for the newly selected candidate right away
            renderTwinCheck(); // its acentric/inversion-twin note depends on the selected space group
          });
        });
        resolve(true);
        }, 10);
      });
    }

    const eStatResultEl = document.getElementById('eStatResult');

    // Composition for Wilson/E-statistics: prefer an exact whole-cell
    // count from a loaded .ins (SFAC/UNIT), otherwise derive it from the
    // "chemical formula per formula unit" field × Z (same Z logic already
    // used for .ins generation above — reused rather than re-estimated).
    function getUnitCellComposition() {
      if (useExactComposition()) {
        return { composition: lastInsComposition, source: 'exact — SFAC/UNIT from the loaded .ins' };
      }
      const formula = parseFormulaString(formulaInput.value);
      if (formula.length === 0) return null;

      const { z, estimated: zWasEstimated } = currentZ(formula);
      if (!z) return null;

      return {
        composition: formula.map(f => ({ element: f.element, count: f.count * z })),
        source: `formula × Z=${z}${zWasEstimated ? ' (Z estimated from cell volume — check it)' : ''}`,
      };
    }


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
          eStatResultEl.innerHTML = `<p style="color:#b00">Need the unit cell's element composition first — either drop the matching
            SHELX .ins (exact SFAC/UNIT counts), or fill in the chemical formula field on the File &amp; Cell tab (also needs Z, e.g. from a loaded CIF) — the statistics then appear here automatically.</p>`;
          resolve(false);
          return;
        }

        eStatResultEl.innerHTML = `<p><em>Merging reflections and fitting the Wilson plot…</em></p>`;

        setTimeout(() => {
        const classification = classifyCellMetric(values);
        const monoclinicAxis = classification.system.includes('c-unique') ? 'c' : 'b';
        const targetLaueClass = determinedLaueClass || classification.laueClasses[0];
        const group = generateGroup(getLaueClassGenerators(targetLaueClass, monoclinicAxis));

        // Centering absences are zero by construction (lattice, not structure):
        // exclude them (observed lattice) before merging.
        const latticeInfo = observedLatticeFilter();
        const stats = computeEStatistics(latticeInfo.reflections, values, compInfo.composition, group, latticeInfo.allowedTest);

        if (stats.error) {
          eStatResultEl.innerHTML = `<p style="color:#b00">${stats.error}</p>`;
          resolve(false);
          return;
        }

        const ref = WILSON_REFERENCE;
        const verdictColor = stats.verdict === 'centrosymmetric' ? '#080' : '#06c';
        const nzColor = stats.nz.closest === 'centrosymmetric' ? '#080' : (stats.nz.closest === 'acentric' ? '#06c' : '#b60');
        // Plain-language reading of the N(z) result next to the <|E²−1|> verdict.
        let nzNote = '';
        if (stats.nz.closest === 'acentric, twinned?') {
          nzNote = `<p style="color:#b60"><small>⚠ The distribution is closest to a <strong>twinned acentric</strong> crystal. A true centrosymmetric structure would be closer to the green curve —
            if the &lt;|E²−1|&gt; verdict above says "centrosymmetric", consider twinning (e.g. a pseudo-merohedral twin of an acentric structure) before trusting it.</small></p>`;
        } else if ((stats.nz.closest === 'centrosymmetric') !== (stats.verdict === 'centrosymmetric')) {
          nzNote = `<p style="color:#b60"><small>⚠ N(z) and &lt;|E²−1|&gt; point in different directions — the data are not clean enough for either to be decisive (twinning, a dominant heavy atom or pseudosymmetry are the usual causes).</small></p>`;
        }

        // If the space-group candidates included a centrosymmetric/acentric
        // tie at the top score, this verdict IS the evidence needed to
        // break it — apply it directly instead of leaving the person to
        // notice and click a radio button themselves.
        let autoSelectNote = '';
        if (lastMatchResults && lastMatchResults.length > 0) {
          const topScore = lastMatchResults[0].score;
          const tiedGroup = lastMatchResults.filter(m => m.score === topScore);
          const matchesVerdict = m => hallSymbolIsCentrosymmetric(m.entry.hall) === (stats.verdict === 'centrosymmetric');
          // Keep the current pick (e.g. preselected from the CIF) if it already fits the verdict.
          const current = lastMatchResults[selectedCandidateIndex];
          const resolved = (current && tiedGroup.includes(current) && matchesVerdict(current)) ? current : tiedGroup.find(matchesVerdict);
          const wasActuallyTied = tiedGroup.some(m => !matchesVerdict(m));
          if (current && selectionSource === 'cif' && tiedGroup.includes(current) && !matchesVerdict(current) && wasActuallyTied) {
            // The CIF names a refined space group; E² is a heuristic that
            // pseudo-symmetry and heavy atoms mislead (nm75_26: P n a 21 with a
            // pseudo n-glide reads "centrosymmetric"). Keep the CIF's pick and say so.
            autoSelectNote = `<p style="color:#b60"><strong>E² points to a ${stats.verdict} structure, but the CIF reports ${formatHM(current.entry.hm)}</strong>
              (${hallSymbolIsCentrosymmetric(current.entry.hall) ? 'centrosymmetric' : 'acentric'}) — the CIF's choice is kept. A refined structure outranks the intensity statistics,
              which pseudo-symmetry, heavy atoms and twinning can mislead; select another row on the Space Group tab if you want to try it.</p>`;
          } else if (resolved && wasActuallyTied) {
            const newIndex = lastMatchResults.indexOf(resolved);
            if (!(selectionSource === 'cif' && resolved === current)) selectionSource = 'auto-e-stat'; // this pick is now backed by E² statistics either way
            if (newIndex !== selectedCandidateIndex) {
              selectedCandidateIndex = newIndex;
              generateIns(); // rebuild the .ins for the now-resolved candidate right away
              renderTwinCheck(); // its acentric/inversion-twin note depends on the selected space group
              const radio = document.querySelector(`input[name="candidateSelect"][value="${newIndex}"]`);
              if (radio) radio.checked = true;
              autoSelectNote = `<p style="color:#080"><strong>✓ Candidate selection updated on the Space Group tab:</strong> switched the space-group pick to
                ${formatHM(resolved.entry.hm)} (${resolved.entry.hall}) to match this ${stats.verdict} result — the SHELX .ins tab has been updated to match.</p>`;
            } else {
              autoSelectNote = `<p><small>This confirms the already-selected candidate (${formatHM(resolved.entry.hm)}${selectionSource === 'cif' ? ', preselected from the CIF' : ''}).</small></p>`;
            }
          }
        }

        // --- per-shell E² statistics (table + plot container + reading) ---
        let eShellHtml = '';
        const es = stats.eShells;
        if (es && !es.tooFew && es.shells.length >= 2) {
          const col = v => (v === 'centrosymmetric' ? '#080' : '#06c');
          const rows = es.shells.map(sh => `
            <tr${sh.few ? ' style="opacity:0.6;"' : ''}><td>${sh.dMax === null ? '∞' : sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)} Å</td><td>${sh.n}</td><td>${sh.meanE2raw.toFixed(2)}</td>
            <td><strong>${sh.meanAbs.toFixed(3)}</strong></td>
            <td style="color:${sh.few ? 'inherit' : col(sh.verdict)}">${sh.few ? 'few reflections' : sh.verdict}</td>
            <td>${(sh.fractionEgt3 * 100).toFixed(2)}%</td></tr>`).join('');
          const usedShells = es.shells.filter(sh => !sh.few);
          const first = usedShells[0], last = usedShells[usedShells.length - 1];
          const fewNote = es.nFew ? ` ${es.nFew} shell(s) with fewer than ${E_SHELL_MIN} reflections are shown (hollow circles) but not used for this comparison.` : '';
          const shellNote = es.consistent
            ? `<p style="color:#080"><small>✓ All ${usedShells.length} resolution shells used point to the same class (${first.verdict}) — the verdict does not depend on the resolution range.${fewNote}</small></p>`
            : `<p style="color:#b60"><small>⚠ The shells do not agree: ${first.verdict} at low resolution (${first.dMax === null ? '∞' : first.dMax.toFixed(2)} – ${first.dMin.toFixed(2)} Å), ${last.verdict} at high resolution
               (${last.dMax.toFixed(2)} – ${last.dMin.toFixed(2)} Å). A heavy atom that dominates the low-resolution data, or pseudosymmetry, can do that — treat the overall verdict with caution.
               If only the outermost shell differs, weak, noisy data are the more likely reason.${fewNote}</small></p>`;
          eShellHtml = `
          <h2>E² statistics by resolution shell</h2>
          <details class="explain"><summary>How to read this</summary><p><small>The same ⟨|E²−1|⟩ as above, worked out separately for resolution shells
            (same shell boundaries as the merging table on the Laue tab, so row n covers the same resolution range in both; shells with fewer than ${E_SHELL_MIN} reflections are not used for the comparison). Within a shell, E² is rescaled to ⟨E²⟩ = 1 so that the single overall Wilson fit does not leak into the comparison;
            "⟨E²⟩ (global fit)" shows what the overall fit predicted before that rescaling — values far from 1 mark shells the straight-line Wilson fit misses.
            Closest = the reference (0.736 acentric / 0.968 centrosymmetric) that the shell value is nearer to.</small></p></details>
          <table>
            <tr><th>d range</th><th>Reflections</th><th>⟨E²⟩ (global fit)</th><th>⟨|E²−1|⟩</th><th>Closest</th><th>|E| &gt; 3</th></tr>
            ${rows}
          </table>
          ${shellNote}
          <div id="eShellPlotContainer"></div>`;
        }

        eStatResultEl.innerHTML = `
          <h2>E² statistics (Wilson plot)</h2>
          <p><small>Composition: ${compInfo.source} · merged under Laue class ${targetLaueClass}${determinedLaueClassSource === 'merging' ? '' : determinedLaueClassSource === 'cif' ? ' (from the CIF space group)' : ' (assumed, not yet confirmed above)'}
          · ${stats.uniqueReflections} unique reflections (${stats.usableReflections} usable)<br>${latticeInfo.note}</small></p>
          <table>
            <tr><th></th><th>this data</th><th>ref. centrosymmetric</th><th>ref. acentric</th></tr>
            <tr><td>&lt;|E²−1|&gt;</td><td><strong>${stats.meanE2minus1.toFixed(3)}</strong></td><td>${ref.centrosymmetric.meanE2minus1}</td><td>${ref.acentric.meanE2minus1}</td></tr>
            <tr><td>&lt;|E|&gt;</td><td>${stats.meanE.toFixed(3)}</td><td>${ref.centrosymmetric.meanE}</td><td>${ref.acentric.meanE}</td></tr>
            <tr><td>fraction |E|&gt;3</td><td>${(stats.fractionEgt3 * 100).toFixed(3)}%</td><td>${(ref.centrosymmetric.fractionEgt3 * 100).toFixed(2)}%</td><td>${(ref.acentric.fractionEgt3 * 100).toFixed(2)}%</td></tr>
            <tr><td>&lt;E²&gt;</td><td>${stats.meanE2.toFixed(3)}</td><td colspan="2"><small>should be ≈1.000 if the Wilson scaling fit well</small></td></tr>
          </table>
          <p style="color:${verdictColor}"><strong>Best match: ${stats.verdict}</strong> <small>(margin ${(stats.margin * 100).toFixed(0)}% between the two reference distances — near 0% means right on the fence)</small></p>
          ${autoSelectNote}
          <p><small>Wilson fit: B = ${stats.B.toFixed(2)} Å² (overall isotropic falloff), scale k = ${stats.k.toFixed(3)}.</small></p>
          <details class="explain"><summary>How to read these statistics</summary><p><small>If &lt;E²&gt; above is far from 1.000, treat the verdict cautiously —
          common real causes are a dominant heavy atom (breaks the "randomly distributed atoms" Wilson assumption, especially at low resolution), twinning, or non-ideal data scaling; none of those are
          errors in this calculation, they're properties of the data itself. This is a different question from the absences scan above — Friedel's law makes every diffraction pattern look centrosymmetric
          regardless of the crystal's true symmetry, so this is the one of the two that can actually distinguish an acentric structure from a centrosymmetric one.</small></p></details>
          <div id="wilsonPlotContainer"></div>
          <h2>N(z) test (cumulative E² distribution)</h2>
          <details class="explain"><summary>How to read the N(z) test</summary><p><small>Shows the fraction of reflections with E² ≤ z against the ideal curves. Unlike the single number &lt;|E²−1|&gt; it uses the whole distribution,
          so it also reveals <strong>twinning</strong> (a twinned acentric crystal sits between the two ideal curves, with too few weak reflections for a centric structure).
          Closeness is the Kolmogorov distance D = largest vertical gap to each curve (smaller = better). With real data D rarely drops below ~0.02–0.05.</small></p></details>
          <table>
            <tr><th></th><th>centrosymmetric</th><th>acentric</th><th>acentric, perfect twin</th></tr>
            <tr><td>D</td><td>${stats.nz.dCentro.toFixed(3)}</td><td>${stats.nz.dAcentric.toFixed(3)}</td><td>${stats.nz.dTwin.toFixed(3)}</td></tr>
          </table>
          <p style="color:${nzColor}"><strong>N(z) is closest to: ${stats.nz.closest}</strong> <small>(${stats.nz.n} reflections${stats.nz.negClipped ? `, ${stats.nz.negClipped} with negative intensity counted as z = 0` : ''})</small></p>
          ${nzNote}
          <div id="nzPlotContainer"></div>
          ${eShellHtml}`;

        renderWilsonPlot(document.getElementById('wilsonPlotContainer'), stats.wilsonPoints, stats.B, stats.k);
        renderNZPlot(document.getElementById('nzPlotContainer'), stats.nz);
        if (eShellHtml) renderEShellPlot(document.getElementById('eShellPlotContainer'), es);
        resolve(true);
        }, 10);
      });
    }

    input.addEventListener('change', (e) => { handleFiles(e.target.files); input.value = ''; });

    // The whole window is a drop target, so more files can be added from
    // any tab without going back to tab 1. Several files at once are fine.
    let dragDepth = 0;
    document.addEventListener('dragenter', (e) => { e.preventDefault(); dragDepth++; dropZone.classList.add('dragover'); });
    document.addEventListener('dragover', (e) => { e.preventDefault(); });
    document.addEventListener('dragleave', () => { dragDepth = Math.max(0, dragDepth - 1); if (dragDepth === 0) dropZone.classList.remove('dragover'); });
    document.addEventListener('drop', (e) => {
      e.preventDefault();
      dragDepth = 0;
      dropZone.classList.remove('dragover');
      if (e.dataTransfer && e.dataTransfer.files.length) handleFiles(e.dataTransfer.files);
    });
