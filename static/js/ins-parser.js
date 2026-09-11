// ins-parser.js
// Reads a SHELX .ins/.res file for exactly what this project can use:
// the unit cell, wavelength, Z, and — directly from SFAC + UNIT — the
// unit cell's element composition. UNIT already gives atom counts for
// the WHOLE cell (SHELX's own convention), so this sidesteps the
// "chemical formula per formula unit, times Z" indirection the CIF path
// needs; when a .ins is available it's the more direct, exact source.
//
// Deliberately narrow: this project only ever needs these four lines,
// not a full .ins/.res interpreter (no atom list, no refinement
// instructions).

// A .ins/.res reliably has CELL and SFAC/UNIT — reflection files (.hkl)
// never do, and CIFs use different tag syntax entirely for their OWN
// cell/composition (_cell_length_a, _chemical_formula_sum, ...). But a
// CIF can legitimately EMBED a full SHELX .res as a text block
// (_shelx_res_file) alongside its own _refln loop or _shelx_hkl_file —
// that embedded block also contains real CELL/SFAC/UNIT lines, so a
// plain substring/regex scan across the WHOLE document would wrongly
// flag such a CIF as a standalone .ins and divert it away from its own
// reflection data entirely (and its own _cell_length_a would never even
// be read). CIFs always go through their own path instead — this must
// never fire for anything looksLikeCIF already claims.
function looksLikeSHELXIns(text) {
  if (looksLikeCIF(text)) return false;
  return /^\s*CELL\s+\S/im.test(text) && /^\s*SFAC\s+\S/im.test(text) && /^\s*UNIT\s+\S/im.test(text);
}

// SHELX line continuation: a line ending in "=" continues on the next
// line. Rare in these four instructions but cheap to handle correctly.
function joinContinuations(text) {
  const lines = text.split(/\r?\n/);
  const joined = [];
  let buffer = '';
  for (const line of lines) {
    const trimmed = line.replace(/\s+$/, '');
    if (trimmed.endsWith('=')) {
      buffer += trimmed.slice(0, -1) + ' ';
    } else {
      joined.push(buffer + trimmed);
      buffer = '';
    }
  }
  if (buffer) joined.push(buffer);
  return joined;
}

function findInstructionLine(lines, keyword) {
  const re = new RegExp(`^\\s*${keyword}\\s+(.*)$`, 'i');
  for (const line of lines) {
    const m = line.match(re);
    if (m) return m[1].trim();
  }
  return null;
}

// Parses CELL, ZERR, SFAC and UNIT into { cell, wavelength, z, composition }.
// Returns null if CELL/SFAC/UNIT aren't all present and consistent —
// callers should fall back to manual entry / the CIF-formula path in
// that case, same as extractCellFromCIF's contract.
function parseInsFile(text) {
  const lines = joinContinuations(text);

  const cellLine = findInstructionLine(lines, 'CELL');
  const sfacLine = findInstructionLine(lines, 'SFAC');
  const unitLine = findInstructionLine(lines, 'UNIT');
  if (!cellLine || !sfacLine || !unitLine) return null;

  // CELL wavelength a b c alpha beta gamma
  const cellParts = cellLine.split(/\s+/).map(Number);
  if (cellParts.length < 7 || cellParts.some(v => !Number.isFinite(v))) return null;
  const [wavelength, a, b, c, alpha, beta, gamma] = cellParts;

  const zerrLine = findInstructionLine(lines, 'ZERR');
  let z = null;
  let esd = null;
  if (zerrLine) {
    const zerrParts = zerrLine.split(/\s+/).map(Number);
    z = Number.isFinite(zerrParts[0]) ? zerrParts[0] : null;
    // ZERR Z sd(a) sd(b) sd(c) sd(alpha) sd(beta) sd(gamma) — an all-zero
    // set is SHELX's own placeholder for "not determined", not a real esd.
    if (zerrParts.length >= 7 && zerrParts.slice(1, 7).every(Number.isFinite)) {
      const allZero = zerrParts.slice(1, 7).every(v => v === 0);
      if (!allZero) {
        esd = {
          a: zerrParts[1], b: zerrParts[2], c: zerrParts[3],
          alpha: zerrParts[4], beta: zerrParts[5], gamma: zerrParts[6],
        };
      }
    }
  }

  // SFAC lists element symbols (possibly with scattering-factor overrides
  // we don't need); UNIT lists one atom count per SFAC entry, same order.
  // Inline custom scattering-factor coefficients after a symbol (rare)
  // aren't handled — SFAC is expected to be a plain list of symbols.
  const elements = sfacLine.split(/\s+/).filter(Boolean);
  const counts = unitLine.split(/\s+/).map(Number);
  if (counts.length < elements.length || elements.some(e => !/^[A-Za-z]{1,2}$/.test(e))) return null;

  const composition = elements.map((symbol, i) => ({
    // normalize to plain element symbol (e.g. "PD" -> "Pd") — SHELX SFAC
    // is case-insensitive and often all-caps
    element: symbol.length > 1 ? symbol[0].toUpperCase() + symbol.slice(1).toLowerCase() : symbol.toUpperCase(),
    count: counts[i],
  })).filter(e => Number.isFinite(e.count) && e.count > 0);

  if (composition.length === 0) return null;

  return {
    cell: { a, b, c, alpha, beta, gamma },
    wavelength: Number.isFinite(wavelength) && wavelength > 0 ? wavelength : null,
    z: Number.isFinite(z) ? z : null,
    composition,
    esd,
  };
}
