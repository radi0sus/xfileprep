// report-writer.js
// Tab 8: a report of the current state, as a self-contained HTML file (print
// it to PDF from the browser) and as a plain-text listing in the spirit of
// XPREP's .prp file (saved as <name>_report.txt).
//
// Nothing is calculated here: every number comes from what the other steps
// stored in `reportState` (app.js) or from the app's own state, so the report
// always shows exactly what the tabs show. The report is a snapshot: it holds
// numbers, tables and the plots as images, not the reflections.
//
// Structure: reportCollect() builds a list of sections made of simple blocks
// (key/value list, table, paragraph, list, plots, preformatted text). Two
// renderers turn the same list into HTML and into the text listing, so both
// always carry the same content.

const XFILEPREP_BUILD = '2026-11-26';

// ---- SHA-256 (synchronous; crypto.subtle is missing on some file:// pages) --
const _SHA256_K = new Uint32Array([
  0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
  0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
  0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
  0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
  0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
  0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
  0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
  0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]);
function sha256Hex(data) {
  const bytes = typeof data === 'string' ? new TextEncoder().encode(data) : data;
  const len = bytes.length;
  const total = Math.ceil((len + 9) / 64) * 64;
  const buf = new Uint8Array(total);
  buf.set(bytes);
  buf[len] = 0x80;
  const dv = new DataView(buf.buffer);
  dv.setUint32(total - 8, Math.floor(len / 0x20000000));
  dv.setUint32(total - 4, (len * 8) >>> 0);
  const h = new Uint32Array([0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]);
  const w = new Uint32Array(64);
  const rotr = (x, n) => (x >>> n) | (x << (32 - n));
  for (let off = 0; off < total; off += 64) {
    for (let i = 0; i < 16; i++) w[i] = dv.getUint32(off + 4 * i);
    for (let i = 16; i < 64; i++) {
      const a = w[i - 15], b = w[i - 2];
      const s0 = rotr(a, 7) ^ rotr(a, 18) ^ (a >>> 3);
      const s1 = rotr(b, 17) ^ rotr(b, 19) ^ (b >>> 10);
      w[i] = (w[i - 16] + s0 + w[i - 7] + s1) | 0;
    }
    let a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (let i = 0; i < 64; i++) {
      const S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
      const ch = (e & f) ^ (~e & g);
      const t1 = (hh + S1 + ch + _SHA256_K[i] + w[i]) | 0;
      const S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
      const maj = (a & b) ^ (a & c) ^ (b & c);
      const t2 = (S0 + maj) | 0;
      hh = g; g = f; f = e; e = (d + t1) | 0; d = c; c = b; b = a; a = (t1 + t2) | 0;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }
  return Array.from(h, x => x.toString(16).padStart(8, '0')).join('');
}

// ---- Small formatting helpers ------------------------------------------------
const rptEsc = s => String(s ?? '').replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const rptNum = (x, d = 2) => (x === null || x === undefined || !Number.isFinite(x) ? '—' : x.toFixed(d));
const rptPct = (x, d = 1) => (x === null || x === undefined || !Number.isFinite(x) ? '—' : (x * 100).toFixed(d) + ' %');

// HTML snippet (as the tabs show it) → plain Unicode text.
function rptHtmlToText(html) {
  if (html === null || html === undefined) return '';
  let s = String(html)
    .replace(/<span class="sgBar">(\d)<\/span>/g, '-$1')
    .replace(/<span class="tag tag-[a-z]+"[^>]*>([^<]*)<\/span>/g, '[$1]')
    .replace(/<br\s*\/?>/gi, ' ')
    .replace(/<\/(p|li|tr|div|h\d)>/gi, ' ');
  const div = document.createElement('div');
  div.innerHTML = s;
  return div.textContent.replace(/[  ]/g, ' ').replace(/\s+/g, ' ').trim();
}

// Plain Unicode → ASCII for the text listing (files keep plain ASCII, as the
// .ins does). Space-group symbols are written from the ASCII H-M symbol.
const _RPT_ASCII = {
  'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta', 'θ': 'theta', 'λ': 'lambda', 'μ': 'mu', 'ν': 'nu',
  'σ': 'sig', 'χ': 'chi', 'ρ': 'rho', 'π': 'pi', 'Δ': 'Delta', 'ε': 'eps', 'τ': 'tau',
  '≈': '~', '≤': '<=', '≥': '>=', '≠': '!=', '±': '+/-', '×': 'x', '·': '.', '−': '-', '–': '-', '—': '-',
  '⟨': '<', '⟩': '>', '→': '->', '←': '<-', '↓': 'v', '°': ' deg', 'Å': 'A', '′': "'", '″': "''",
  '²': '^2', '³': '^3', '¹': '^1', '⁻': '^-', '½': '1/2', '⅓': '1/3', '⅔': '2/3', '¼': '1/4', '¾': '3/4',
  '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5', '∞': 'inf', '⚠': '!', '✓': '*', '✔': '*', '▶': '*', '✕': 'x',
  'ℹ': 'i', '‘': "'", '’': "'", '“': '"', '”': '"', '…': '...', '∗': '*', 'ⓘ': '(i)', '●': '*', '◎': 'o', '√': 'sqrt',
};
function rptAscii(text) {
  return String(text ?? '')
    .replace(/(\d)̅/g, '-$1') // "4̅" (formatSGText) → "-4"
    .replace(/⁻([¹²³])/g, (m, d) => '^-' + { '¹': '1', '²': '2', '³': '3' }[d])
    .replace(/[^\x00-\x7f]/g, ch => (_RPT_ASCII[ch] !== undefined ? _RPT_ASCII[ch] : ch.normalize('NFD').replace(/[^\x00-\x7f]/g, '') || '?'));
}
const rptT = html => rptAscii(rptHtmlToText(html)); // HTML → ASCII text

// value(su) from a value and its esd, as entered when there is no esd
function rptWithSu(value, su) {
  if (!Number.isFinite(value)) return '—';
  if (!(su > 0)) return String(value);
  return fmtSU(value, su);
}

function rptFileSize(n) {
  if (!Number.isFinite(n)) return '—';
  return n < 1024 ? `${n} B` : n < 1048576 ? `${(n / 1024).toFixed(1)} kB` : `${(n / 1048576).toFixed(2)} MB`;
}
function rptDate(d) {
  if (!d) return '—';
  const p = x => String(x).padStart(2, '0');
  return `${d.getFullYear()}-${p(d.getMonth() + 1)}-${p(d.getDate())} ${p(d.getHours())}:${p(d.getMinutes())}`;
}

// Provenance tag: who made a decision / where a value came from.
const RPT_TAG = {
  auto: ['auto', 'determined automatically from the data'],
  user: ['by hand', 'set or chosen by the user'],
  file: ['file', 'taken from a loaded file'],
  cif: ['CIF', 'taken from the loaded CIF'],
  assumed: ['assumed', 'no source — a default was used'],
  estimated: ['estimated', 'estimated, not measured'],
};
const rptTag = kind => `<span class="tag tag-${kind}" title="${RPT_TAG[kind][1]}">${RPT_TAG[kind][0]}</span>`;

const RPT_SELECTION = {
  default: ['auto', 'top-ranked candidate of the absence analysis'],
  'auto-e-stat': ['auto', 'tie between a centrosymmetric and an acentric candidate broken by the E² statistics'],
  manual: ['user', 'selected by hand on the Space Group tab'],
  cif: ['cif', 'the space group the loaded CIF reports (the data agree with it)'],
};

// ---- Plots as images ------------------------------------------------------------
// The tabs draw on canvas in the page's colour scheme; the report always uses
// the light scheme (it is a document to print), at twice the resolution.
function rptPlot(title, draw) {
  const canvas = document.createElement('canvas');
  const savedDark = window.plotIsDarkMode;
  try {
    window.plotIsDarkMode = () => false;
    draw(canvas, 2);
    return { title, src: canvas.toDataURL('image/png'), width: parseFloat(canvas.style.width) / 2 };
  } catch (err) {
    return null;
  } finally {
    window.plotIsDarkMode = savedDark;
  }
}

// ---- Collecting the content -----------------------------------------------------
function reportHasData() {
  return !!lastReflections;
}

function reportBaseName() {
  const f = loadedFiles.reflections || loadedFiles.cell;
  const name = f ? f.name : (lastLoadedFileBaseName || 'structure');
  return name.replace(/\.[^./\\]+$/, '') || 'structure';
}

function reportCollect(opts) {
  const S = reportState;
  const sections = [];
  const findings = [];
  const add = (level, html, ref) => findings.push({ level, html, ref });

  const cell = currentCellOrNull();
  const lam = currentWavelength();
  const formula = parseFormulaString(formulaInput.value);
  const zi = currentZ(formula);
  const sgRes = S.sg && S.sg.results ? S.sg.results : null;
  const sel = sgRes && sgRes[selectedCandidateIndex] ? sgRes[selectedCandidateIndex] : null;
  const selEntry = sel ? sel.entry : null;
  const shell = S.shell && S.shell.stats && !S.shell.stats.error ? S.shell : null;
  const estat = S.estat && S.estat.stats ? S.estat : null;
  const laue = S.laue;
  const dash = '—';

  // ---- provenance of the inputs and decisions ----
  const cellSrc = loadedFiles.cell
    ? (cellEditedByHand ? `${rptTag('user')} edited after loading ${rptEsc(loadedFiles.cell.name)}` : `${rptTag('file')} ${rptEsc(loadedFiles.cell.name)}`)
    : (cell ? `${rptTag('user')} entered by hand` : dash);
  const lamSrc = lam.source === 'file' ? `${rptTag('file')} ${rptEsc(lambdaState.label)}`
    : lam.source === 'manual' ? rptTag('user') : rptTag('assumed');
  const exactComp = useExactComposition();
  const formulaSrc = exactComp ? `${rptTag('file')} exact SFAC/UNIT of the .ins`
    : formula.length === 0 ? dash
    : formulaManuallyEdited ? rptTag('user')
    : (loadedFiles.cell ? `${rptTag('file')} ${rptEsc(loadedFiles.cell.name)}` : rptTag('user'));
  const zSrc = !zi.z ? dash : zi.estimated ? `${rptTag('estimated')} V / (18 Å³ × non-H atoms)`
    : zState.source === 'file' ? `${rptTag('file')} ${rptEsc(zState.label)}` : rptTag('user');
  const laueSrc = !laue ? dash
    : laue.source === 'cif' ? `${rptTag('cif')} point group of the CIF space group`
    : laue.source === 'merging' ? `${rptTag('auto')} merging statistics (R<sub>int</sub>)`
    : `${rptTag('assumed')} highest-symmetry candidate of the cell metric — not confirmed`;
  const laueClass = laue ? (laue.determined || (sgRes ? S.sg.targetLaueClass : null)) : null;
  const selSrc = selEntry ? (() => { const [k, t] = RPT_SELECTION[selectionSource] || RPT_SELECTION.default; return `${rptTag(k)} ${t}`; })() : dash;

  // ---- findings: everything worth a look, most severe first ----
  if (!cell) add('bad', 'No complete unit cell — nothing could be analysed.', null);
  const tf = twinFindings() || [];
  const twinNotInIns = !!(S.ins && S.ins.notes && S.ins.notes.twin);
  for (const f of tf) add('warn', f.kind === 'nonmero' ? `${f.text} — needs an HKLF 5 file (PLATON TwinRotMat, TWINABS); the .ins cannot model it.`
    : `${f.kind === 'inversion' ? 'Inversion twin' : 'Twinning'}: ${f.text}${twinNotInIns && S.ins.notes.twin.includes(f.text) ? ' — <em>not written to the .ins</em> (tick it on tab 3 to add TWIN/BASF)' : ''}`, 'sec-twin');
  for (const f of absFindings()) add(f.text.includes('inverted') ? 'bad' : 'warn', `Absolute structure: ${f.text}`, 'sec-twin');
  for (const sec of S.twinSections || []) {
    if (sec.level === 'bad' && sec.id !== 'abs') add('bad', `${sec.short}: ${sec.verdict}`, 'sec-twin');
  }
  if (S.ins && S.ins.twin) add('info', `The .ins contains <code>TWIN ${S.ins.twin.matrix.flat().join(' ')}${S.ins.twin.n !== 2 ? ' ' + S.ins.twin.n : ''}</code> with <code>BASF ${S.ins.twin.basf.map(b => b.toFixed(2)).join(' ')}</code> (ticked on tab 3).`, 'sec-ins');
  if (S.ins && S.ins.tiedOpposite && selectionSource === 'default') {
    add('bad', `${formatSG(selEntry.hm)} is tied with ${formatSG(S.ins.tiedOpposite.entry.hm)} (${hallSymbolIsCentrosymmetric(S.ins.tiedOpposite.entry.hall) ? 'centrosymmetric' : 'acentric'}) — absences cannot decide; no E² statistics were available to break the tie.`, 'sec-sg');
  }
  if (S.sg && S.sg.unresolved && S.sg.unresolved.length) {
    add('warn', `${formatSG(sgRes[0].entry.hm)} is not clearly separated from ${S.sg.unresolved.map(r => `${formatSG(r.entry.hm)} (${(r.score * 100).toFixed(1)} %)`).join(', ')} — treat the top row as a proposal.`, 'sec-sg');
  }
  if (S.sg && S.sg.sameAbsences && S.sg.sameAbsences.length > 1) {
    add('info', `Identical systematic absences: ${S.sg.sameAbsences.map(r => formatSG(r.entry.hm)).join(', ')} — only solution and refinement can decide between them.`, 'sec-sg');
  }
  const pt = pseudoTranslationInfo();
  if (pt && pt.strongFound && pt.strongFound.length) add('warn', `Pseudo-translation: ${pseudoTranslationSummary(pt)}. E² looks more centrosymmetric and the twin fraction is overestimated.`, 'sec-abs');
  if (laue && !laue.source) add('warn', `Laue class not confirmed — the data have no redundant reflections; ${formatPG(laueClass)} is assumed from the cell metric.`, 'sec-laue');
  if (estat && selEntry) {
    const selCentro = hallSymbolIsCentrosymmetric(selEntry.hall);
    if ((estat.stats.verdict === 'centrosymmetric') !== selCentro) {
      add(selectionSource === 'cif' ? 'info' : 'warn', `E² statistics point to ${estat.stats.verdict === 'acentric' ? 'an' : 'a'} ${estat.stats.verdict} structure, the selected ${formatSG(selEntry.hm)} is ${selCentro ? 'centrosymmetric' : 'acentric'}.`, 'sec-estat');
    }
  }
  if (estat && estat.nzNote) add('warn', rptHtmlToText(estat.nzNote).replace(/^⚠\s*/, ''), 'sec-estat');
  if (estat && estat.stats.eShells && !estat.stats.eShells.tooFew && estat.stats.eShells.shells.length >= 2 && !estat.stats.eShells.consistent) {
    add('warn', 'The E² verdict changes with resolution (low vs. high resolution shells disagree).', 'sec-estat');
  }
  if (S.estat && S.estat.error && sgRes) add('info', `No E² statistics: ${rptEsc(S.estat.error)}.`, 'sec-estat');
  if (shell) {
    const ov = shell.stats.overall;
    if (ov.completeness !== null && ov.completeness < 0.95) add('warn', `Overall completeness ${rptPct(ov.completeness)} (to ${ov.dMin.toFixed(2)} Å).`, 'sec-data');
    if (shell.stats.cutoffIndex >= 0) {
      const sh = shell.stats.shells[shell.stats.cutoffIndex];
      add('info', `⟨I/σ⟩ drops below 2 at ${sh.dMax === null ? '∞' : sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)} Å.`, 'sec-data');
    }
  }
  const cr = S.crystal;
  if (cr && cr.cd.absorption && (cr.cd.absorption.level === 'strong' || cr.cd.absorption.level === 'very strong')) {
    add('warn', `Absorption is ${cr.cd.absorption.level}: μ·r up to ${cr.cd.absorption.muRMax.toFixed(2)} — ${cr.cd.absorption.level === 'very strong' ? 'a numerical (face-indexed) correction is recommended' : 'a multi-scan correction is needed'}.`, 'sec-cell');
  }
  if (cell && lam.source === 'assumed') add('warn', `Wavelength assumed (${lam.value} Å, Mo Kα) — none of the loaded files gives one.`, 'sec-cell');
  if (cell && formula.length === 0 && !exactComp) add('warn', 'No chemical formula — no E² statistics, density, F(000) or μ.', 'sec-cell');
  if (zi.z && zi.estimated) add('warn', `Z = ${zi.z} is estimated from the cell volume.`, 'sec-cell');
  if (cellEditedByHand && loadedFiles.cell) add('info', 'The cell constants were edited by hand after loading.', 'sec-cell');
  if (S.ins && S.ins.settingTransform) add('info', `The cell is transformed from the measured setting (${formatSG(S.ins.settingTransform.nativeHM)}); the .ins carries the matrix on the HKLF line.`, 'sec-sg');
  const order = { bad: 0, warn: 1, info: 2 };
  findings.sort((a, b) => order[a.level] - order[b.level]);

  // ---- Summary ----
  const zp = zi.z ? selectedZPrime(zi.z) : null;
  const result = [];
  result.push(['Space group', selEntry ? `<strong class="big">${formatSG(selEntry.hm)}</strong> &nbsp;No. ${selEntry.number}${selEntry.qualifier ? ':' + rptEsc(selEntry.qualifier) : ''} · Hall <code>${rptEsc(selEntry.hall)}</code> · ${hallSymbolIsCentrosymmetric(selEntry.hall) ? 'centrosymmetric' : 'acentric'}` : '<em>not determined</em>', selSrc]);
  if (S.ins && S.ins.settingTransform) result.push(['Setting', `transformed from ${formatSG(S.ins.settingTransform.nativeHM)} (as measured) — HKLF matrix <code>${rptEsc(formatHKLFMatrix(S.ins.settingTransform.matrix))}</code>`, `${rptTag('auto')}`]);
  result.push(['Laue class', laueClass ? `${formatPG(laueClass)} <small>(cell metric: ${rptEsc(laue.system)})</small>` : dash, laueSrc]);
  if (cell) result.push(['Cell', `a = ${rptWithSu(cell.a, lastCellEsd && lastCellEsd.a)}, b = ${rptWithSu(cell.b, lastCellEsd && lastCellEsd.b)}, c = ${rptWithSu(cell.c, lastCellEsd && lastCellEsd.c)} Å<br>α = ${rptWithSu(cell.alpha, lastCellEsd && lastCellEsd.alpha)}, β = ${rptWithSu(cell.beta, lastCellEsd && lastCellEsd.beta)}, γ = ${rptWithSu(cell.gamma, lastCellEsd && lastCellEsd.gamma)}° · V = ${cellVolume(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma).toFixed(1)} Å³`, cellSrc]);
  result.push(['Wavelength', lam.value ? `${lam.value} Å` : dash, lamSrc]);
  result.push(['Formula', exactComp ? `${rptEsc(formulaInput.value || '')} <small>(whole cell: ${lastInsComposition.map(c => c.element + c.count).join(' ')})</small>` : (formula.length ? rptEsc(formulaInput.value) : dash), formulaSrc]);
  result.push(['Z / Z′', zi.z ? `Z = ${zi.z}${zp ? ` · ${rptEsc(zPrimeText(zp).replace(/^Z′ = /, 'Z′ = ').replace(/ in .*$/, ''))} <small>(${zp.mult} general positions)</small>` : ''}` : dash, zSrc]);
  sections.push({
    id: 'sec-summary', title: 'Summary', blocks: [
      { t: 'kv3', head: ['', 'Result', 'Source / decision'], rows: result },
      { t: 'h', html: 'Points to check' },
      findings.length
        ? { t: 'findings', items: findings }
        : { t: 'p', level: 'ok', html: 'Nothing flagged by the checks.' },
    ],
  });

  // ---- 1 Input files ----
  const fileRows = [];
  const seen = new Set();
  for (const [role, f] of [['Reflections', loadedFiles.reflections], ['Cell / composition', loadedFiles.cell]]) {
    if (!f) continue;
    const m = f.meta || {};
    const key = f.name + '|' + (m.sha256 || '');
    const roles = seen.has(key) ? null : role;
    if (seen.has(key)) { fileRows.forEach(r => { if (r._key === key) r[0] += ' + ' + role.toLowerCase(); }); continue; }
    seen.add(key);
    const row = [roles, rptEsc(f.name), rptFileSize(m.size), rptDate(m.modified), m.sha256 ? `<code class="hash">${m.sha256}</code>` : dash];
    row._key = key;
    fileRows.push(row);
  }
  const hk = summarizeHKL(lastReflections);
  sections.push({
    id: 'sec-files', title: 'Input files', blocks: [
      { t: 'table', head: ['Role', 'File', 'Size', 'Modified', 'SHA-256'], rows: fileRows.map(r => r.slice()) },
      { t: 'kv', rows: [
        ['Reflection data', `${hk.count} reflections${hk.count ? ` · h ${hk.hRange.join(' … ')}, k ${hk.kRange.join(' … ')}, l ${hk.lRange.join(' … ')}` : ''}`],
      ] },
      { t: 'p', small: true, html: 'The SHA-256 checksum identifies the exact file (compare with <code>sha256sum</code> or <code>Get-FileHash</code>).' },
    ],
  });

  // ---- 2 Cell & crystal data ----
  const crystalRows = [];
  if (cell) {
    crystalRows.push(['Crystal system (metric)', rptEsc(laue ? laue.system : classifyCellMetric(cell).system)]);
    crystalRows.push(['Cell volume', `${cellVolume(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma).toFixed(1)} Å³`]);
  }
  if (cr) {
    const cd = cr.cd;
    crystalRows.push(['Composition used', `${cr.composition.map(c => c.element + (Number.isInteger(c.count) ? c.count : c.count.toFixed(2))).join(' ')} <small>(${rptEsc(cr.compNote)})</small>`]);
    if (cd.formulaWeight) crystalRows.push(['Formula weight', `${cd.formulaWeight.toFixed(2)} g/mol`]);
    crystalRows.push(['Density (calc.)', `${cd.density.toFixed(3)} g/cm³`]);
    crystalRows.push(['F(000)', `${Math.round(cd.f000)}`]);
    if (cr.det) crystalRows.push(['Z estimate', `V / 18 Å³ → Z ≈ ${cr.det.zRaw.toFixed(2)} → ${cr.det.z}`]);
    if (cd.mu !== null) crystalRows.push([`μ (${cd.anode} Kα)`, `${cd.mu.toFixed(3)} mm⁻¹ <small>(NIST FFAST)</small>`]);
    const sz = (cr.size || []).filter(x => x > 0);
    if (sz.length === 3) crystalRows.push(['Crystal size', `${sz.map(x => x.toFixed(2)).join(' × ')} mm`]);
    if (cd.absorption) {
      const a = cd.absorption;
      crystalRows.push(['μ·r (min / max)', `${a.muRMin.toFixed(2)} / ${a.muRMax.toFixed(2)} — ${rptEsc(a.level)}`]);
      crystalRows.push(['T<sub>min</sub> / T<sub>max</sub> (expected)', `${a.tMin.toFixed(3)} / ${a.tMax.toFixed(3)}`]);
    }
    const fpps = cd.perElement.filter(p => p.fpp !== null).sort((x, y) => y.fpp - x.fpp);
    if (fpps.length) crystalRows.push([`f″ (${rptEsc(cd.fppLine)})`, fpps.map(p => `${p.element} ${fmtFpp(p.fpp)}`).join(' · ')]);
  }
  sections.push({
    id: 'sec-cell', title: 'Cell and crystal data', blocks: crystalRows.length
      ? [{ t: 'kv', rows: crystalRows }]
      : [{ t: 'p', html: 'Needs the cell (and a formula for density, F(000) and μ).' }],
  });

  // ---- 3 Data quality ----
  const dataBlocks = [];
  if (shell) {
    const ov = shell.stats.overall;
    const thetaMax = lam.value && ov.dMin ? Math.asin(Math.min(1, lam.value / (2 * ov.dMin))) * 180 / Math.PI : null;
    dataBlocks.push({ t: 'kv', rows: [
      ['Merged under', `Laue class ${formatPG(shell.laue)}; ${rptEsc(shell.lattice)}`],
      ['Resolution', `${ov.dMax === null ? '∞' : ov.dMax.toFixed(2)} – ${ov.dMin.toFixed(3)} Å${thetaMax !== null ? ` · θ<sub>max</sub> = ${thetaMax.toFixed(2)}°` : ''} · sin θ/λ ≤ ${(1 / (2 * ov.dMin)).toFixed(3)} Å⁻¹`],
      ['Reflections', `${ov.nTotal} measured, ${ov.nUnique} unique · redundancy ${rptNum(ov.redundancy)}`],
      ['Completeness', rptPct(ov.completeness)],
      ['⟨I/σ⟩ · I > 2σ', `${rptNum(ov.meanISigma, 1)} · ${rptPct(ov.fractionGt2)}`],
      ['R<sub>int</sub> · R<sub>σ</sub>', `${rptPct(ov.rInt, 2)} · ${rptPct(ov.rSigma, 2)}`],
      ['R<sub>meas</sub> · R<sub>pim</sub> · CC<sub>½</sub>', `${rptPct(ov.rMeas, 2)} · ${rptPct(ov.rPim, 2)} · ${rptPct(ov.cc12)}`],
    ] });
    if (opts.shellTable) {
      const dL = sh => `${sh.dMax === null ? '∞' : sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)}`;
      const row = sh => [dL(sh), sh.nUnique, rptPct(sh.completeness), rptNum(sh.redundancy), rptNum(sh.meanISigma, 1), rptPct(sh.fractionGt2), rptPct(sh.cc12), rptPct(sh.rInt), rptPct(sh.rMeas), rptPct(sh.rPim), rptPct(sh.rSigma)];
      dataBlocks.push({
        t: 'table', caption: 'Resolution shells', num: true,
        head: ['d (Å)', 'Unique', 'Compl.', 'Red.', '⟨I/σ⟩', 'I > 2σ', 'CC<sub>½</sub>', 'R<sub>int</sub>', 'R<sub>meas</sub>', 'R<sub>pim</sub>', 'R<sub>σ</sub>'],
        rows: [...shell.stats.shells.map(row), ['all', ...row(shell.stats.overall).slice(1)]],
        mark: shell.stats.cutoffIndex >= 0 ? [shell.stats.cutoffIndex] : [], total: true,
        note: shell.stats.cutoffIndex >= 0 ? 'Marked row: ⟨I/σ⟩ drops below 2 there (orientation only — nothing was cut).' : null,
      });
    }
    if (opts.plots && shell.stats.shells.length >= 2) {
      const specs = shellCurveSpecs(shell.stats);
      dataBlocks.push({ t: 'plots', items: [
        rptPlot('Shell statistics (percent)', (cv, z) => drawShellCurve(cv, shell.stats, specs[0], z)),
        rptPlot('⟨I/σ⟩ per shell', (cv, z) => drawShellCurve(cv, shell.stats, specs[1], z)),
      ].filter(Boolean) });
    }
  } else {
    dataBlocks.push({ t: 'p', html: S.shell && S.shell.stats && S.shell.stats.error ? rptEsc(S.shell.stats.error) : 'Needs reflections and all six cell constants.' });
  }
  sections.push({ id: 'sec-data', title: 'Data quality', blocks: dataBlocks });

  // ---- 4 Laue symmetry ----
  const laueBlocks = [];
  if (laue) {
    laueBlocks.push({
      t: 'table', num: true, head: ['Laue class', 'Order', 'Unique refl.', 'Redundancy', 'R<sub>int</sub>'],
      rows: laue.evaluations.map(r => [`${formatPG(r.laueClass)}${r.laueClass === laue.chosen ? ' ✓' : ''}`, r.groupOrder, r.uniqueReflections, r.redundancy.toFixed(2), r.rInt !== null ? (r.rInt * 100).toFixed(2) + ' %' : '—']),
      mark: [laue.evaluations.findIndex(r => r.laueClass === laue.chosen)].filter(i => i >= 0),
    });
    laueBlocks.push({ t: 'p', html: laue.source === 'cif'
      ? `Laue class ${formatPG(laue.chosen)} from the space group of the CIF; the R<sub>int</sub> table is an independent check.`
      : laue.source === 'merging'
        ? `${formatPG(laue.chosen)}: highest-symmetry class without an R<sub>int</sub> jump (50 % tolerance).${laue.jumpDetectedAt ? ` R<sub>int</sub> jumps at ${formatPG(laue.jumpDetectedAt)}.` : ''}`
        : 'No redundant reflections under any candidate — merging statistics cannot choose; the highest-symmetry class of the metric is assumed.',
    level: laue.source ? null : 'warn' });
  } else laueBlocks.push({ t: 'p', html: 'Not determined (needs reflections and the cell).' });
  sections.push({ id: 'sec-laue', title: 'Laue symmetry', blocks: laueBlocks });

  // ---- 5 Twins & absolute structure ----
  const twinBlocks = [];
  if (S.twinSections) {
    twinBlocks.push({ t: 'table', head: ['Check', 'Result', 'Verdict'],
      rows: S.twinSections.map(s => [rptEsc(s.short), s.result || '—', `<span class="lv-${s.level}">${s.verdict}</span>`]) });
    const tc = lastTwinCheck;
    if (tc) twinBlocks.push({ t: 'p', small: true, html: `Laue class ${formatPG(tc.laue)} · lattice symmetry ${formatPG(tc.holohedry)}${tc.holohedryNote ? ` (${rptEsc(tc.holohedryNote)})` : ''}` });
    if (tc && tc.applicable && tc.laws.length) {
      twinBlocks.push({ t: 'table', caption: 'Merohedral twin laws (H test)', num: true,
        head: ['Twin law (SHELX)', 'Obliquity', 'Pairs', '⟨H⟩', 'α', 'Verdict'],
        rows: tc.laws.map(L => [`<code>${L.shelxMatrix.flat().join(' ')}</code>`, L.deviation ? `${L.deviation.angle.toFixed(2)}°` : '—', L.pairs,
          L.meanH === null ? '—' : L.meanH.toFixed(3), L.alpha === null ? '—' : `${L.alpha.toFixed(2)} ± ${L.alphaSe.toFixed(3)}`, rptEsc(L.verdict)]) });
      if (opts.plots) twinBlocks.push({ t: 'plots', items: [rptPlot('N(H) — H test', (cv, z) => drawNHPlot(cv, tc.laws, z))].filter(Boolean) });
    }
    const nm = lastNonMeroCheck;
    if (nm && nm.applicable && nm.laws.length) {
      twinBlocks.push({ t: 'table', caption: 'Non-merohedral twin laws (overlapping reflections)', num: true,
        head: ['2-fold about', 'Obliquity', 'Overlapping', 'Pairs', 'α (overlaps)', 'Verdict'],
        rows: nm.laws.map(L => [L.about === 'direct' ? `[${reducedDirection(L.axis).join(' ')}]` : `(${reducedDirection(L.plane).join(' ')})*`, `${L.delta.toFixed(2)}°`, L.overlaps, L.pairs,
          L.alpha === null ? '—' : `${L.alpha.toFixed(2)} ± ${L.alphaSe.toFixed(3)}`, rptEsc(L.verdict)]) });
    }
    const an = lastAnomalousCheck;
    if (an && an.enough) {
      const rows = [['Bijvoet pairs', `${an.pairs}${an.calibrated ? '' : ' <small>(σ not calibrated — merged data)</small>'}`], ['Anomalous signal χ²', `${an.chi2.toFixed(2)} <small>(1 = noise only)</small>`]];
      if (an.twin) {
        rows.push(['Expected χ² (no inversion twin)', an.twin.expectedChi2.toFixed(1)]);
        rows.push(['Inversion-twin fraction x', `${an.twin.x.toFixed(2)} ± ${an.twin.xSe.toFixed(2)} — ${rptEsc(an.twin.verdict)} <small>(model-free, no sign)</small>`]);
      }
      twinBlocks.push({ t: 'kv', caption: 'Inversion twin from the Bijvoet pairs', rows });
    }
    const ab = lastAbsStructure;
    if (ab && ab.applicable) {
      const row = (name, x, su, n) => [name, fmtSU(x, su), n ?? '—', rptEsc(absStructureVerdict(x, su).text)];
      const rows = [];
      if (ab.parsons) rows.push(row('Parsons quotients', ab.parsons.x, ab.parsons.su, ab.parsons.n));
      if (ab.hooft) rows.push(row('Hooft y', ab.hooft.y, ab.hooft.su, ab.hooft.n));
      if (ab.flack) rows.push(row('Flack, classical (model fixed)', ab.flack.x, ab.flack.su, ab.flack.n));
      if (ab.reported && ab.reported.flack) rows.push(['In the CIF', rptEsc(ab.reported.flack.raw), '—', '']);
      twinBlocks.push({ t: 'table', caption: 'Absolute structure (refined model from the CIF)', head: ['Method', 'x (su)', 'Pairs / refl.', 'Verdict'], rows });
    }
    if (pt && pt.strongFound && pt.strongFound.length) twinBlocks.push({ t: 'p', level: 'warn', html: 'A pseudo-translation is present: pseudo-symmetry raises the twin fraction of the H test — refine BASF and judge by the refinement.' });
    twinBlocks.push({ t: 'p', small: true, html: S.ins && S.ins.twin
      ? `Written to the .ins: <code>TWIN ${S.ins.twin.matrix.flat().join(' ')}${S.ins.twin.n !== 2 ? ' ' + S.ins.twin.n : ''}</code>, <code>BASF ${S.ins.twin.basf.map(b => b.toFixed(2)).join(' ')}</code>.`
      : 'Nothing from these checks was written to the .ins.' });
  } else twinBlocks.push({ t: 'p', html: 'Not run (needs reflections and the cell).' });
  sections.push({ id: 'sec-twin', title: 'Twins and absolute structure', blocks: twinBlocks });

  // ---- 6 Absences ----
  const absBlocks = [];
  const A = S.absences;
  if (A && A.scan) {
    const scan = A.scan;
    const vtxt = r => `${{ holds: 'satisfied', violated: 'violated', inconclusive: 'undetermined' }[r.verdict] || r.verdict}${r.weakTail ? ' (weak tail)' : ''}`;
    const vcell = r => `<span class="v-${r.verdict}">${vtxt(r)}</span>`;
    const counts = r => r.forbiddenCount > 0 && Number.isFinite(r.forbiddenPresentCount) ? `${r.forbiddenPresentCount} / ${r.forbiddenCount}` : (r.forbiddenCount === 0 ? '0 / 0' : '—');
    const CENT = { A: 'k+l = 2n', B: 'h+l = 2n', C: 'h+k = 2n', I: 'h+k+l = 2n', F: 'h+k, h+l, k+l = 2n', 'R (obverse)': '−h+k+l = 3n', 'R (reverse)': 'h−k−l = 3n' };
    absBlocks.push({ t: 'kv', rows: [['Observed lattice', rptEsc(observedCenteringName(scan.centering))]] });
    absBlocks.push({ t: 'table', caption: 'Lattice centring', head: ['Type', 'Condition', 'Forbidden > 3σ', 'Status'], num: false,
      rows: Object.entries(scan.centering).filter(([n]) => CENT[n]).map(([n, r]) => [n, CENT[n], `${counts(r)} <small>(${r.violationPercent.toFixed(1)} %)</small>`, vcell(r)]) });
    const axRows = [];
    for (const [axis, results] of Object.entries(scan.axial)) {
      for (const r of results) axRows.push([axis, `${{ h00: 'h', '0k0': 'k', '00l': 'l' }[axis] || '?'} = ${r.period}n`, r.reflectionCount, `${counts(r)} <small>(${r.violationPercent.toFixed(1)} %)</small>`, vcell(r)]);
    }
    absBlocks.push({ t: 'table', caption: 'Axial conditions (screw axes)', head: ['Row', 'Condition', 'N refl.', 'Forbidden > 3σ', 'Status'], rows: axRows });
    const znRows = [];
    for (const [zone, results] of Object.entries(scan.zonal)) {
      results.forEach((r, i) => {
        const pseudo = A.relevance && r.verdict === 'holds' && A.relevance.zonal[zone] && !A.relevance.zonal[zone][i];
        znRows.push([zone, r.name, r.reflectionCount, `${counts(r)} <small>(${r.violationPercent.toFixed(1)} %)</small>`, vcell(r) + (pseudo ? ` <small class="lv-warn">pseudo-symmetry — impossible in ${formatPG(A.laue)}</small>` : '')]);
      });
    }
    absBlocks.push({ t: 'table', caption: 'Zonal conditions (glide planes)', head: ['Zone', 'Condition', 'N refl.', 'Forbidden > 3σ', 'Status'], rows: znRows,
      note: 'Forbidden > 3σ: forbidden reflections still significant, and their share. “Weak tail”: satisfied although some forbidden reflections are significant — only the weak tail of a residual intensity (median I/σ < 3).' });
    if (pt) {
      const V = { 'pseudo-translation': 'pseudo-translation', 'weak hint': 'weak hint', absent: 'absent (sub-cell or unrecognised centring?)' };
      const shown = pt.table.filter(t => V[t.verdict]);
      absBlocks.push(shown.length
        ? { t: 'table', caption: 'Pseudo-translations (weak, but present)', head: ['Type', 't', 'Weak class', 'N weak / rest', '⟨I⟩ weak/rest', 'Verdict'],
            rows: shown.map(t => [ptLabel(t.label), t.vector, t.cond, `${t.stats.nWeak} / ${t.stats.nRest}`, `${(100 * t.stats.ratio).toFixed(0)} % / ${(100 * t.stats.medianRatio).toFixed(0)} %`, `<span class="${t.verdict === 'pseudo-translation' ? 'lv-warn' : ''}">${V[t.verdict]}</span>`]) }
        : { t: 'p', small: true, html: `Pseudo-translations: none — no reflection class is markedly weaker than the rest (${pt.nUsed} reflections without a zero index).` });
    }
    if (opts.zonePlots && cell) {
      absBlocks.push({ t: 'plots', items: Object.keys(ZONE_PLOT_CONFIG).map(z =>
        rptPlot(`${z} (${ZONE_PLOT_CONFIG[z].fixed})`, (cv, zoom) => renderZonePlot(cv, lastReflections, z, scan, cell, zoom, 300))).filter(Boolean), small: true });
      absBlocks.push({ t: 'p', small: true, html: 'Zone plots: green = significant reflection, grey = measured but weak (I/σ < 3), blank = not measured.' });
    }
  } else absBlocks.push({ t: 'p', html: 'Not run (needs reflections).' });
  sections.push({ id: 'sec-abs', title: 'Systematic absences', blocks: absBlocks });

  // ---- 7 E² statistics ----
  const eBlocks = [];
  if (estat) {
    const st = estat.stats, ref = WILSON_REFERENCE;
    eBlocks.push({ t: 'p', small: true, html: `Composition: ${rptEsc(estat.compInfo.source)} · merged under ${formatPG(estat.laue)} · ${st.uniqueReflections} unique reflections (${st.usableReflections} usable) · ${rptEsc(estat.lattice)}` });
    eBlocks.push({ t: 'table', num: true, head: ['', 'This data', 'Centrosymmetric', 'Acentric'], rows: [
      ['⟨|E²−1|⟩', `<strong>${st.meanE2minus1.toFixed(3)}</strong>`, ref.centrosymmetric.meanE2minus1, ref.acentric.meanE2minus1],
      ['⟨|E|⟩', st.meanE.toFixed(3), ref.centrosymmetric.meanE, ref.acentric.meanE],
      ['|E| > 3', (st.fractionEgt3 * 100).toFixed(3) + ' %', (ref.centrosymmetric.fractionEgt3 * 100).toFixed(2) + ' %', (ref.acentric.fractionEgt3 * 100).toFixed(2) + ' %'],
      ['⟨E²⟩', st.meanE2.toFixed(3), '1.000', '1.000'],
    ] });
    eBlocks.push({ t: 'kv', rows: [
      ['⟨|E²−1|⟩ verdict', `<strong>${st.verdict}</strong> <small>(margin ${(st.margin * 100).toFixed(0)} %)</small>`],
      ['N(z) closest to', `<strong>${rptEsc(st.nz.closest)}</strong> <small>(D: centro ${st.nz.dCentro.toFixed(3)} · acentric ${st.nz.dAcentric.toFixed(3)} · twinned ${st.nz.dTwin.toFixed(3)})</small>`],
      ['Wilson fit', `B = ${st.B.toFixed(2)} Å² · k = ${st.k.toFixed(3)}`],
    ] });
    if (estat.autoSelect) eBlocks.push({ t: 'p', small: true, html: estat.autoSelect.replace(/<\/?p[^>]*>/g, '') });
    const es = st.eShells;
    if (es && !es.tooFew && es.shells.length >= 2) {
      eBlocks.push({ t: 'table', caption: 'E² by resolution shell', num: true, head: ['d (Å)', 'Refl.', '⟨|E²−1|⟩', 'Closest'],
        rows: es.shells.map(sh => [`${sh.dMax === null ? '∞' : sh.dMax.toFixed(2)} – ${sh.dMin.toFixed(2)}`, sh.n, sh.meanAbs.toFixed(3), sh.few ? '(few reflections)' : sh.verdict]) });
    }
    if (opts.plots) {
      eBlocks.push({ t: 'plots', items: [
        rptPlot('Wilson plot', (cv, z) => drawWilsonPlot(cv, st.wilsonPoints, st.B, st.k, z)),
        rptPlot('N(z) test', (cv, z) => drawNZPlot(cv, st.nz, z)),
        es && !es.tooFew && es.shells.length >= 2 ? rptPlot('⟨|E²−1|⟩ per shell', (cv, z) => drawEShellPlot(cv, es, z)) : null,
      ].filter(Boolean) });
    }
  } else eBlocks.push({ t: 'p', html: S.estat && S.estat.error ? `Not computed: ${rptEsc(S.estat.error)}.` : 'Not computed (needs reflections, the cell and a composition).' });
  sections.push({ id: 'sec-estat', title: 'E² statistics', blocks: eBlocks });

  // ---- 8 Space group ----
  const sgBlocks = [];
  if (sgRes) {
    sgBlocks.push({ t: 'p', small: true, html: `Laue class used: ${formatPG(S.sg.targetLaueClass)} (${rptEsc(S.sg.laueClassNote)}).` });
    const shown = sgRes.slice(0, 10).map((r, i) => ({ r, i }));
    if (S.sg.extraRow) shown.push({ r: S.sg.extraRow.result, i: S.sg.extraRow.index });
    sgBlocks.push({ t: 'table', head: ['', 'No.', 'H-M symbol', 'Hall symbol', 'Match', 'Setting'],
      rows: shown.map(({ r, i }) => [i === selectedCandidateIndex ? '▶' : '', `${r.entry.number}${r.entry.qualifier ? ':' + rptEsc(r.entry.qualifier) : ''}`, formatSG(r.entry.hm), `<code>${rptEsc(r.entry.hall)}</code>`, `${(r.score * 100).toFixed(1)} %`, S.sg.settingOf(r.entry)]),
      mark: shown.map((x, k) => (x.i === selectedCandidateIndex ? k : -1)).filter(k => k >= 0),
      sgTable: true,
      note: `▶ = selected: ${selSrc}. Match = share of the reflection conditions that agree with the scan, weighted by the reflections behind each check.` });
    if (S.sg.tiedOpposite) sgBlocks.push({ t: 'p', level: selectionSource === 'default' ? 'bad' : 'info', html: `${formatSG(sgRes[0].entry.hm)} and ${formatSG(S.sg.tiedOpposite.entry.hm)} are tied (centrosymmetric vs. acentric) — absences cannot tell them apart; the E² statistics can.` });
    if (S.sg.sameAbsences && S.sg.sameAbsences.length > 1) sgBlocks.push({ t: 'p', level: 'info', html: `Identical systematic absences: ${S.sg.sameAbsences.map(r => formatSG(r.entry.hm)).join(', ')}.` });
    if (S.sg.unresolved && S.sg.unresolved.length) sgBlocks.push({ t: 'p', level: 'warn', html: `Not clearly separated from the top row: ${S.sg.unresolved.map(r => `${formatSG(r.entry.hm)} (${(r.score * 100).toFixed(1)} %)`).join(', ')}.` });
    if (S.sg.runnerUps && S.sg.runnerUps.length) sgBlocks.push({ t: 'list', caption: 'Why close runners-up score lower', items: S.sg.runnerUps.map(x => `<strong>${formatSG(x.r.entry.hm)}</strong> (${(x.r.score * 100).toFixed(1)} %): ${x.bad.join('; ')}`) });
    if (zp && zp.note) sgBlocks.push({ t: 'p', level: 'warn', html: rptEsc(zp.note) });
  } else sgBlocks.push({ t: 'p', html: 'No candidates (needs reflections and the cell).' });
  sections.push({ id: 'sec-sg', title: 'Space group', blocks: sgBlocks });

  // ---- Appendix ----
  if (opts.ins && S.ins && S.ins.text) {
    const insNotes = ['missing', 'lambda', 'z'].map(k => S.ins.notes[k]).filter(Boolean).map(rptHtmlToText);
    sections.push({ id: 'sec-ins', title: 'Appendix A — SHELX .ins', appendix: true, blocks: [
      ...insNotes.map(n => ({ t: 'p', level: 'warn', html: rptEsc(n) })),
      { t: 'pre', text: S.ins.text, summary: `${reportBaseName()}.ins for ${formatSGText(S.ins.entry.hm)}` },
    ] });
  }
  sections.push({ id: 'sec-ref', title: `Appendix ${opts.ins && S.ins && S.ins.text ? 'B' : 'A'} — Methods and references`, appendix: true, blocks: [
    { t: 'list', small: true, items: [
      'Symmetry from Hall symbols: S. R. Hall, <i>Acta Cryst.</i> A37, 517 (1981); International Tables Vol. B.',
      'Lattice symmetry: Y. Le Page, <i>J. Appl. Cryst.</i> 15, 255 (1982). H test: T. O. Yeates, <i>Acta Cryst.</i> A44, 142 (1988).',
      'CC<sub>½</sub> (σ–τ method): Assmann, Brehm &amp; Diederichs, <i>J. Appl. Cryst.</i> 49, 1021 (2016).',
      'Absolute structure: Flack &amp; Bernardinelli, <i>J. Appl. Cryst.</i> 33, 1143 (2000); Parsons, Flack &amp; Wagner, <i>Acta Cryst.</i> B69, 249 (2013); Hooft, Straver &amp; Spek, <i>J. Appl. Cryst.</i> 41, 96 (2008).',
      'f″: International Tables Vol. C, Table 4.2.6.8. μ/ρ: NIST FFAST (Chantler et al., 2005). Atomic weights: IUPAC.',
    ] },
  ] });

  return {
    title: reportBaseName(),
    generated: new Date(),
    notes: opts.notes || '',
    selected: selEntry,
    sections,
  };
}

// ---- HTML renderer --------------------------------------------------------------
const REPORT_CSS = `
:root { color-scheme: light; --ink:#1a1d23; --muted:#5f6673; --line:#d9dde3; --soft:#f4f5f7; --accent:#2f4fc7;
  --bad:#b3261e; --bad-bg:#fdecea; --warn:#8a5a00; --warn-bg:#fff5df; --info:#1f5fa8; --info-bg:#eaf2fc; --ok:#147a3d; --ok-bg:#e8f5ed; }
* { box-sizing: border-box; }
html { background: #e9ebef; }
body { margin: 0; color: var(--ink); background: #fff; font: 14px/1.5 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
  max-width: 980px; margin: 24px auto; padding: 40px 48px 56px; box-shadow: 0 1px 4px rgba(0,0,0,.08); }
header.rh { border-bottom: 2px solid var(--ink); padding-bottom: 14px; margin-bottom: 8px; }
.rh .kicker { font-size: 12px; letter-spacing: .08em; text-transform: uppercase; color: var(--muted); }
.rh h1 { font-size: 26px; margin: 2px 0 6px; font-weight: 650; }
.rh .meta { font-size: 12.5px; color: var(--muted); display: flex; flex-wrap: wrap; gap: 4px 18px; }
.rh .notes { margin-top: 10px; padding: 8px 12px; background: var(--soft); border-radius: 6px; white-space: pre-wrap; font-size: 13px; }
nav.toc { font-size: 12.5px; margin: 10px 0 4px; display: flex; flex-wrap: wrap; gap: 2px 14px; }
nav.toc a { color: var(--accent); text-decoration: none; }
nav.toc a:hover { text-decoration: underline; }
section { margin-top: 26px; }
section > h2 { font-size: 17px; margin: 0 0 10px; padding-bottom: 4px; border-bottom: 1px solid var(--line); font-weight: 650; }
section > h2 { display: flex; align-items: baseline; }
section > h2 a.top { margin-left: auto; font-size: 11.5px; font-weight: 500; color: var(--muted); text-decoration: none; white-space: nowrap; }
section > h2 a.top:hover { color: var(--accent); text-decoration: underline; }
section > h2 .no { color: var(--muted); font-weight: 500; margin-right: 6px; }
h3 { font-size: 13.5px; margin: 16px 0 6px; font-weight: 650; }
p { margin: 6px 0; }
small, .small { font-size: 12px; color: var(--muted); }
code, pre { font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace; font-size: 12px; }
code { background: var(--soft); padding: 0 3px; border-radius: 3px; }
code.hash { font-size: 10.5px; word-break: break-all; background: none; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: 6px 0 4px; font-size: 13px; break-inside: avoid; }
th, td { text-align: left; padding: 4px 8px; border-bottom: 1px solid var(--line); vertical-align: top; }
th { font-size: 12px; font-weight: 600; color: var(--muted); background: var(--soft); }
table.num td:not(:first-child), table.num th:not(:first-child) { text-align: right; font-variant-numeric: tabular-nums; }
tr.mark td { background: #fff7d6; font-weight: 600; }
tr.total td { border-top: 1.5px solid var(--ink); font-weight: 600; }
table.kv td:first-child { width: 30%; color: var(--muted); }
table.kv3 td:first-child { width: 17%; color: var(--muted); }
table.kv3 td:last-child { width: 34%; font-size: 12.5px; }
table.kv3 .big { font-size: 19px; }
.caption { font-size: 12.5px; font-weight: 600; margin: 14px 0 2px; }
.tnote { font-size: 11.5px; color: var(--muted); margin: 2px 0 8px; }
.tag { display: inline-block; font-size: 10.5px; font-weight: 600; letter-spacing: .02em; padding: 0 6px; border-radius: 9px; margin-right: 4px;
  border: 1px solid; line-height: 16px; vertical-align: 1px; white-space: nowrap; }
.tag-auto { color: var(--info); border-color: #b8cfee; background: var(--info-bg); }
.tag-user { color: #6a2fa8; border-color: #d6c3ee; background: #f4edfb; }
.tag-file, .tag-cif { color: var(--ok); border-color: #b9dcc6; background: var(--ok-bg); }
.tag-assumed, .tag-estimated { color: var(--warn); border-color: #ecd29a; background: var(--warn-bg); }
ul.findings { list-style: none; padding: 0; margin: 6px 0; }
ul.findings li { padding: 6px 10px 6px 12px; margin: 4px 0; border-left: 4px solid; border-radius: 0 5px 5px 0; }
ul.findings li .lvl { font-size: 10.5px; font-weight: 700; text-transform: uppercase; letter-spacing: .05em; margin-right: 8px; }
ul.findings li a { font-size: 11.5px; color: var(--muted); margin-left: 6px; text-decoration: none; white-space: nowrap; }
.l-bad { border-color: var(--bad); background: var(--bad-bg); } .l-bad .lvl { color: var(--bad); }
.l-warn { border-color: #d9a43a; background: var(--warn-bg); } .l-warn .lvl { color: var(--warn); }
.l-info { border-color: #8fb3e3; background: var(--info-bg); } .l-info .lvl { color: var(--info); }
p.l-ok, p.l-bad, p.l-warn, p.l-info { padding: 6px 10px; border-left: 4px solid; border-radius: 0 5px 5px 0; }
p.l-ok { border-color: var(--ok); background: var(--ok-bg); color: var(--ok); }
.lv-bad, .v-violated { color: var(--bad); } .lv-warn { color: var(--warn); } .lv-ok, .v-holds { color: var(--ok); } .v-inconclusive { color: var(--warn); }
.lv-none { color: var(--muted); }
ul.plain { margin: 4px 0; padding-left: 20px; } ul.plain li { margin: 2px 0; }
.plots { display: flex; flex-wrap: wrap; gap: 14px 18px; margin: 10px 0; }
.plots figure { margin: 0; break-inside: avoid; }
.plots figure.small img { width: 210px; }
.plots img { display: block; border: 1px solid var(--line); border-radius: 4px; cursor: zoom-in; max-width: 100%; height: auto; }
.plots figcaption { font-size: 11.5px; color: var(--muted); margin-top: 3px; }
details { margin: 6px 0; }
details > summary { cursor: pointer; font-size: 12.5px; color: var(--accent); }
pre { background: var(--soft); border: 1px solid var(--line); border-radius: 6px; padding: 10px 12px; white-space: pre; overflow-x: auto; line-height: 1.35; }
.sgSym { white-space: nowrap; font-style: normal; } .sgSym i { font-style: italic; }
.sgSym sub { font-size: .72em; vertical-align: -.25em; line-height: 0; } .sgBar { text-decoration: overline; }
footer { margin-top: 36px; padding-top: 8px; border-top: 1px solid var(--line); font-size: 11.5px; color: var(--muted); display: flex; justify-content: space-between; }
#lightbox { position: fixed; inset: 0; background: rgba(20,22,26,.78); display: none; align-items: center; justify-content: center; z-index: 10; cursor: zoom-out; padding: 24px; }
#lightbox.on { display: flex; }
#lightbox figure { margin: 0; background: #fff; padding: 12px; border-radius: 8px; max-width: 96vw; max-height: 94vh; overflow: auto; }
#lightbox img { display: block; max-width: calc(96vw - 24px); max-height: calc(94vh - 60px); }
#lightbox figcaption { font-size: 12.5px; color: var(--muted); margin-top: 6px; }
@media (max-width: 700px) { body { margin: 0; padding: 20px 16px 40px; box-shadow: none; } table { font-size: 12px; } th, td { padding: 3px 5px; } }
@page { size: A4; margin: 15mm 13mm 16mm; }
@media print {
  html { background: #fff; }
  body { margin: 0; padding: 0; max-width: none; box-shadow: none; font-size: 10.5pt; }
  nav.toc, details > summary, #lightbox, a.top { display: none !important; }
  section { break-inside: auto; }
  section > h2 { break-after: avoid; }
  section#sec-ins, section.appendix.first { break-before: page; }
  .blk { break-inside: avoid; } .blk.long { break-inside: auto; }
  .caption, h3 { break-after: avoid; }
  .plots { gap: 8px 4%; }
  .plots figure { width: 48%; } .plots figure img { width: 100% !important; }
  .plots figure.small { width: 23%; }
  table { font-size: 9pt; } th, td { padding: 2.5px 6px; }
  .plots img { cursor: default; }
  a { color: inherit; text-decoration: none; }
  ul.findings li, p.l-ok, p.l-bad, p.l-warn, p.l-info, tr.mark td, th, code, pre, .tag { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
}`;

// Opens collapsed sections for printing (closed <details> do not print) and
// enlarges plots on click.
const REPORT_JS = `
(function () {
  var ds = [];
  window.addEventListener('beforeprint', function () {
    ds = Array.prototype.filter.call(document.querySelectorAll('details'), function (d) { return !d.open; });
    ds.forEach(function (d) { d.open = true; });
  });
  window.addEventListener('afterprint', function () { ds.forEach(function (d) { d.open = false; }); ds = []; });
  var lb = document.getElementById('lightbox');
  document.addEventListener('click', function (e) {
    var img = e.target.closest && e.target.closest('.plots img');
    if (img) {
      lb.querySelector('img').src = img.src;
      lb.querySelector('figcaption').textContent = img.alt;
      lb.classList.add('on');
      return;
    }
    if (e.target.closest && e.target.closest('#lightbox')) lb.classList.remove('on');
  });
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') lb.classList.remove('on'); });
})();`;

// Caption + table (or plots) stay together on a printed page; long tables may break.
function renderReportBlockHTML(b) {
  const html = renderReportBlockInner(b);
  const long = b.rows && b.rows.length > 22;
  return html ? `<div class="blk${long ? ' long' : ''}">${html}</div>` : '';
}
function renderReportBlockInner(b) {
  const lvl = b.level ? ` class="l-${b.level}"` : '';
  switch (b.t) {
    case 'h': return `<h3>${b.html}</h3>`;
    case 'p': return b.small ? `<p class="small"${b.level ? ` style="padding:6px 10px"` : ''}>${b.html}</p>` : `<p${lvl}>${b.html}</p>`;
    case 'kv': return `${b.caption ? `<div class="caption">${b.caption}</div>` : ''}<table class="kv">${b.rows.map(([k, v]) => `<tr><td>${k}</td><td>${v}</td></tr>`).join('')}</table>`;
    case 'kv3': return `<table class="kv3"><tr>${b.head.map(h => `<th>${h}</th>`).join('')}</tr>${b.rows.map(r => `<tr>${r.map(c => `<td>${c}</td>`).join('')}</tr>`).join('')}</table>`;
    case 'table': {
      const mark = new Set(b.mark || []);
      const rows = b.rows.map((r, i) => {
        const cls = [mark.has(i) ? 'mark' : '', b.total && i === b.rows.length - 1 ? 'total' : ''].filter(Boolean).join(' ');
        return `<tr${cls ? ` class="${cls}"` : ''}>${r.map(c => `<td>${c ?? ''}</td>`).join('')}</tr>`;
      }).join('');
      return `${b.caption ? `<div class="caption">${b.caption}</div>` : ''}<table${b.num ? ' class="num"' : ''}><tr>${b.head.map(h => `<th>${h}</th>`).join('')}</tr>${rows}</table>${b.note ? `<p class="tnote">${b.note}</p>` : ''}`;
    }
    case 'findings': return `<ul class="findings">${b.items.map(f => `<li class="l-${f.level}"><span class="lvl">${{ bad: 'Problem', warn: 'Check', info: 'Note' }[f.level]}</span>${f.html}${f.ref ? `<a href="#${f.ref}">→ details</a>` : ''}</li>`).join('')}</ul>`;
    case 'list': return `${b.caption ? `<div class="caption">${b.caption}</div>` : ''}<ul class="plain${b.small ? ' small' : ''}">${b.items.map(i => `<li>${i}</li>`).join('')}</ul>`;
    case 'plots': return b.items.length ? `<div class="plots">${b.items.map(p => `<figure${b.small ? ' class="small"' : ''}><img src="${p.src}" alt="${rptEsc(rptHtmlToText(p.title))}" style="width:${b.small ? 210 : p.width}px" title="Click to enlarge"><figcaption>${p.title}</figcaption></figure>`).join('')}</div>` : '';
    case 'pre': return `<details open><summary>${rptEsc(b.summary || 'show')}</summary><pre>${rptEsc(b.text)}</pre></details>`;
    default: return '';
  }
}

function renderReportHTML(model) {
  let n = 0;
  const secs = model.sections.map(s => {
    const no = s.appendix || s.id === 'sec-summary' ? '' : `<span class="no">${++n}</span>`;
    const firstApp = s.appendix && model.sections.find(x => x.appendix) === s;
    return `<section id="${s.id}"${s.appendix ? ` class="appendix${firstApp ? ' first' : ''}"` : ''}><h2>${no}${s.title}<a class="top" href="#top" title="Back to the top">↑ top</a></h2>${s.blocks.map(renderReportBlockHTML).join('\n')}</section>`;
  }).join('\n');
  const toc = model.sections.map(s => `<a href="#${s.id}">${s.title.replace(/^Appendix [AB] — /, '')}</a>`).join('');
  const sgLine = model.selected ? ` · ${formatSG(model.selected.hm)}` : '';
  return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>${rptEsc(model.title)} — XFilePrep report</title>
<style>${REPORT_CSS}</style>
</head>
<body>
<header class="rh" id="top">
  <div class="kicker">XFilePrep · space-group determination report</div>
  <h1>${rptEsc(model.title)}${sgLine}</h1>
  <div class="meta"><span>Generated ${rptDate(model.generated)}</span><span>XFilePrep build ${XFILEPREP_BUILD}</span><span>Snapshot of the analysis at this time — reload the data in XFilePrep to change it.</span></div>
  ${model.notes ? `<div class="notes">${rptEsc(model.notes)}</div>` : ''}
</header>
<nav class="toc">${toc}</nav>
${secs}
<footer><span>${rptEsc(model.title)} — XFilePrep report</span><span>${rptDate(model.generated)}</span></footer>
<div id="lightbox"><figure><img alt=""><figcaption></figcaption></figure></div>
<script>${REPORT_JS}<\/script>
</body>
</html>`;
}

// ---- Text listing renderer -----------------------------------------------------------
const PRP_WIDTH = 79;

function prpWrap(text, indent = '', first = indent) {
  const words = text.split(/\s+/).filter(Boolean);
  const lines = [];
  let line = first, fresh = true;
  for (const w of words) {
    if (!fresh && (line + ' ' + w).length > PRP_WIDTH) { lines.push(line); line = indent + w; }
    else line = fresh ? line + w : line + ' ' + w;
    fresh = false;
  }
  if (!fresh) lines.push(line);
  return lines.join('\n');
}

function prpTable(head, rows, numeric, marks) {
  const cells = [head.map(rptT), ...rows.map(r => r.map(c => rptT(c === null || c === undefined ? '' : String(c))))];
  const w = head.map((_, j) => Math.max(...cells.map(r => (r[j] || '').length)));
  const fmt = (r, i) => {
    const mark = i > 0 && marks && marks.has(i - 1) ? '>' : ' ';
    return mark + ' ' + r.map((c, j) => (numeric && j > 0 && i > 0 ? (c || '').padStart(w[j]) : (c || '').padEnd(w[j]))).join('  ').replace(/\s+$/, '');
  };
  const out = [fmt(cells[0], 0), '  ' + w.map(x => '-'.repeat(x)).join('  ')];
  cells.slice(1).forEach((r, i) => out.push(fmt(r, i + 1)));
  return out.join('\n');
}

function renderReportBlockText(b) {
  const tagText = html => rptT(String(html).replace(/<span class="tag tag-[a-z]+"[^>]*>([^<]*)<\/span>/g, '[$1]'));
  switch (b.t) {
    case 'h': return `${rptT(b.html)}:`;
    case 'p': return prpWrap((b.level === 'bad' ? '!! ' : b.level === 'warn' ? '! ' : '') + rptT(b.html), '  ', '  ');
    case 'kv': {
      const rows = b.rows.map(([k, v]) => [rptT(k), tagText(v)]);
      const kw = Math.max(...rows.map(r => r[0].length));
      return (b.caption ? `  ${rptT(b.caption)}\n` : '') + rows.map(([k, v]) => prpWrap(v, ' '.repeat(kw + 5), `  ${k.padEnd(kw)} : `)).join('\n');
    }
    case 'kv3': {
      const rows = b.rows.map(r => [rptT(r[0]), tagText(r[1]), tagText(r[2])]);
      const kw = Math.max(...rows.map(r => r[0].length));
      return rows.map(([k, v, s]) => prpWrap(v, ' '.repeat(kw + 5), `  ${k.padEnd(kw)} : `) + (s && s !== '-' ? '\n' + prpWrap(`(${s})`, ' '.repeat(kw + 5), ' '.repeat(kw + 5)) : '')).join('\n');
    }
    case 'table': {
      // the setting column of the candidate table is long HTML: short form in the listing
      const rows = b.sgTable ? b.rows.map(r => [r[0] ? '*' : '', r[1], rptAscii(rptHtmlToText(r[2]).replace(/\s+/g, ' ')), r[3], r[4], rptT(String(r[5]).replace(/<br><code[^>]*>([^<]*)<\/code>/, ' (HKLF $1)')).replace(/\s*!.*$/, '')]) : b.rows;
      return (b.caption ? `  ${rptT(b.caption)}\n` : '') + prpTable(b.head, rows, b.num, new Set(b.mark || [])) + (b.note ? '\n' + prpWrap(rptT(b.note), '  ', '  ') : '');
    }
    case 'findings': return b.items.map(f => prpWrap(rptT(f.html), '            ', `  ${{ bad: '[PROBLEM]', warn: '[CHECK]  ', info: '[NOTE]   ' }[f.level]} `)).join('\n');
    case 'list': return (b.caption ? `  ${rptT(b.caption)}:\n` : '') + b.items.map(i => prpWrap(rptT(i), '    ', '  - ')).join('\n');
    case 'plots': return '';
    case 'pre': return b.text.replace(/\s+$/, '').split('\n').map(l => '  ' + l).join('\n');
    default: return '';
  }
}

function renderReportText(model) {
  const rule = ch => ch.repeat(PRP_WIDTH);
  const out = [];
  const title = `XFilePrep report: ${model.title}`;
  out.push(rule('='));
  out.push(rptAscii(title));
  out.push(`Generated ${rptDate(model.generated)} - XFilePrep build ${XFILEPREP_BUILD}`);
  if (model.selected) out.push(`Space group: ${model.selected.hm.replace(/\s+/g, ' ').trim()} (No. ${model.selected.number})`);
  out.push(rule('='));
  if (model.notes) { out.push(''); model.notes.split('\n').forEach(l => out.push(prpWrap(rptAscii(l), '  ', '  '))); }
  let n = 0;
  for (const s of model.sections) {
    const no = s.appendix || s.id === 'sec-summary' ? '' : `${++n}. `;
    out.push('', '', rptAscii(`${no}${s.title}`).toUpperCase(), rule('-'));
    for (const b of s.blocks) {
      const t = renderReportBlockText(b);
      if (t) out.push('', t);
    }
  }
  out.push('', rule('='), 'Plots are in the HTML report only.', '');
  return out.join('\n').replace(/[ \t]+\n/g, '\n');
}

// ---- Tab 8: preview, options, downloads -------------------------------------------
function reportOptions() {
  const cb = id => { const el = document.getElementById(id); return !!(el && el.checked); };
  const notesEl = document.getElementById('reportNotes');
  return { plots: cb('reportOptPlots'), zonePlots: cb('reportOptZones'), shellTable: cb('reportOptShells'), ins: cb('reportOptIns'), notes: notesEl ? notesEl.value.trim() : '' };
}

let _reportTimer = null;
let _reportLastHtml = null;

function reportIsVisible() {
  const panel = document.getElementById('tab-report');
  return !!(panel && panel.classList.contains('active'));
}

function reportScheduleRefresh() {
  if (!reportIsVisible()) return;
  clearTimeout(_reportTimer);
  _reportTimer = setTimeout(reportRefresh, 250);
}

function reportRefresh() {
  const frame = document.getElementById('reportPreview');
  const status = document.getElementById('reportStatus');
  const buttons = ['reportSaveHtmlBtn', 'reportPrintBtn', 'reportSavePrpBtn'].map(id => document.getElementById(id));
  if (!frame) return;
  if (!reportHasData()) {
    _reportLastHtml = null;
    buttons.forEach(b => { if (b) b.disabled = true; });
    if (status) status.innerHTML = '<small>Load reflection data first — the report shows what the other tabs found.</small>';
    frame.style.display = 'none';
    return;
  }
  let html;
  try {
    html = renderReportHTML(reportCollect(reportOptions()));
  } catch (err) {
    if (status) status.innerHTML = `<span style="color:#b00"><small>Could not build the report: ${rptEsc(err.message)}</small></span>`;
    return;
  }
  _reportLastHtml = html;
  buttons.forEach(b => { if (b) b.disabled = false; });
  frame.style.display = '';
  const y = frame.contentWindow ? frame.contentWindow.scrollY : 0;
  frame.onload = () => { try { frame.contentWindow.scrollTo(0, y); } catch (e) { /* ignore */ } };
  frame.srcdoc = html;
  if (status) status.innerHTML = `<small>Preview of the current state (${rptDate(new Date())}). It follows every change on the other tabs; the saved file is a snapshot.</small>`;
}

function reportDownload(text, filename, type) {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

(function initReportTab() {
  const panel = document.getElementById('tab-report');
  if (!panel) return;
  const btn = document.querySelector('.tabBtn[data-tab="tab-report"]');
  if (btn) btn.addEventListener('click', () => setTimeout(reportRefresh, 0));
  const origSwitch = window.switchTab;
  if (typeof origSwitch === 'function') {
    window.switchTab = id => { origSwitch(id); if (id === 'tab-report') setTimeout(reportRefresh, 0); };
  }
  ['reportOptPlots', 'reportOptZones', 'reportOptShells', 'reportOptIns'].forEach(id => {
    const el = document.getElementById(id);
    if (el) el.addEventListener('change', reportRefresh);
  });
  const notes = document.getElementById('reportNotes');
  if (notes) notes.addEventListener('input', () => { clearTimeout(_reportTimer); _reportTimer = setTimeout(reportRefresh, 600); });

  document.getElementById('reportSaveHtmlBtn').addEventListener('click', () => {
    if (!reportHasData()) return;
    const html = renderReportHTML(reportCollect(reportOptions()));
    reportDownload(html, `${reportBaseName()}_report.html`, 'text/html');
  });
  document.getElementById('reportSavePrpBtn').addEventListener('click', () => {
    if (!reportHasData()) return;
    const text = renderReportText(reportCollect(reportOptions()));
    reportDownload(text.replace(/\n/g, '\r\n'), `${reportBaseName()}_report.txt`, 'text/plain');
  });
  document.getElementById('reportPrintBtn').addEventListener('click', () => {
    const frame = document.getElementById('reportPreview');
    if (!frame || !_reportLastHtml) return;
    // the print dialog of the preview frame ("Save as PDF" there)
    try { frame.contentWindow.focus(); frame.contentWindow.print(); }
    catch (err) {
      const w = window.open('', '_blank');
      if (w) { w.document.write(_reportLastHtml); w.document.close(); w.focus(); w.print(); }
    }
  });
  reportRefresh();
})();
