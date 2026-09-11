// cif-reflections.js
// Extracts reflection data (h,k,l,F2,sigma) from CIF files.
//
// Two supported cases:
//   1) Embedded raw HKL text, stored verbatim under a CIF text field —
//      most commonly "_shelx_hkl_file", which is exactly what tools like
//      ShredCIF pull back out. We just locate that text field and hand it
//      straight to parseHKL().
//   2) A genuine "_refln" loop (proper CIF reflection table with its own
//      column tags), as used e.g. for IUCr deposition.
//
// Requires parseHKL() from hkl-parser.js to be loaded first (plain script,
// no ES modules — keeps everything file:// compatible).

// Extracts a CIF semicolon-delimited text field for a given tag name.
// Returns the raw text content, or null if the tag / text field isn't found.
function extractCIFTextField(cifText, tagName) {
  const lines = cifText.split(/\r?\n/);
  const tagRegex = new RegExp(`^\\s*${tagName}\\s*$`, 'i');

  for (let i = 0; i < lines.length; i++) {
    if (!tagRegex.test(lines[i])) continue;

    let j = i + 1;
    while (j < lines.length && lines[j].trim() === '') j++;

    // a CIF text field opens with a line starting with ';'
    if (j < lines.length && lines[j].startsWith(';')) {
      const contentLines = [];
      const firstLineRemainder = lines[j].slice(1);
      if (firstLineRemainder.trim() !== '') contentLines.push(firstLineRemainder);
      j++;

      while (j < lines.length && !lines[j].startsWith(';')) {
        contentLines.push(lines[j]);
        j++;
      }
      return contentLines.join('\n');
    }
  }
  return null;
}

// Known aliases per reflection quantity — different software uses slightly
// different tag names, so we accept a small set of common variants.
const REFLN_TAG_ALIASES = {
  h: ['_refln_index_h'],
  k: ['_refln_index_k'],
  l: ['_refln_index_l'],
  f2: ['_refln_f_squared_meas', '_refln_intensity_meas'],
  sigma: ['_refln_f_squared_sigma', '_refln_intensity_sigma'],
};

// Parses a genuine "_refln" loop_ block into the same {h,k,l,f2,sigma}
// shape used by parseHKL(). Returns null if no matching loop is found.
function parseReflnLoop(cifText) {
  const lines = cifText.split(/\r?\n/);

  for (let i = 0; i < lines.length; i++) {
    if (!/^\s*loop_\s*$/i.test(lines[i])) continue;

    const tags = [];
    let j = i + 1;
    while (j < lines.length && /^\s*_/.test(lines[j])) {
      tags.push(lines[j].trim().toLowerCase());
      j++;
    }

    const findCol = (aliases) => tags.findIndex(t => aliases.includes(t));
    const hCol = findCol(REFLN_TAG_ALIASES.h);
    const kCol = findCol(REFLN_TAG_ALIASES.k);
    const lCol = findCol(REFLN_TAG_ALIASES.l);
    if (hCol === -1 || kCol === -1 || lCol === -1) continue; // not the refln loop

    const f2Col = findCol(REFLN_TAG_ALIASES.f2);
    const sigmaCol = findCol(REFLN_TAG_ALIASES.sigma);
    const nCols = tags.length;

    // collect data rows until the loop ends (blank line, next loop_/tag, or new data_ block)
    const dataLines = [];
    while (
      j < lines.length &&
      lines[j].trim() !== '' &&
      !/^\s*loop_/i.test(lines[j]) &&
      !/^\s*_/.test(lines[j]) &&
      !/^\s*data_/i.test(lines[j])
    ) {
      dataLines.push(lines[j]);
      j++;
    }

    const tokens = dataLines.join(' ').trim().split(/\s+/).filter(t => t !== '');
    const reflections = [];
    const skipped = [];

    for (let r = 0; r + nCols <= tokens.length; r += nCols) {
      const row = tokens.slice(r, r + nCols);
      const h = parseInt(row[hCol], 10);
      const k = parseInt(row[kCol], 10);
      const l = parseInt(row[lCol], 10);
      const f2 = f2Col >= 0 ? parseFloat(row[f2Col]) : NaN;
      const sigma = sigmaCol >= 0 ? parseFloat(row[sigmaCol]) : NaN;

      if ([h, k, l, f2, sigma].every(Number.isFinite)) {
        reflections.push({ h, k, l, f2, sigma, batch: null });
      } else {
        skipped.push({ lineNumber: null, content: row.join(' ') });
      }
    }

    return { reflections, skipped };
  }
  return null;
}

// Tries embedded raw HKL first (most common for SHELX/Olex2 archives),
// then falls back to a genuine _refln loop. Returns null if neither is found.
function parseCIFReflections(cifText) {
  const embeddedHKL = extractCIFTextField(cifText, '_shelx_hkl_file');
  if (embeddedHKL !== null) {
    const result = parseHKL(embeddedHKL);
    return { ...result, source: 'embedded_shelx_hkl_file' };
  }

  const loopResult = parseReflnLoop(cifText);
  if (loopResult !== null) {
    return { ...loopResult, source: 'refln_loop' };
  }

  return null;
}

// Quick content sniff to decide whether input text is CIF or raw HKL,
// so the caller doesn't have to rely on the file extension alone.
function looksLikeCIF(text) {
  return /^\s*data_/im.test(text) || /_refln_index_h/i.test(text) || /_shelx_hkl_file/i.test(text);
}

// Single entry point: takes raw file text, returns {reflections, skipped, source}
// regardless of whether it was a .hkl or a .cif file.
function parseReflectionInput(text) {
  if (looksLikeCIF(text)) {
    const cifResult = parseCIFReflections(text);
    if (cifResult !== null) return cifResult;
    // sniffed as CIF but no recognizable reflection data inside — fall through
  }
  return { ...parseHKL(text), source: 'raw_hkl' };
}
