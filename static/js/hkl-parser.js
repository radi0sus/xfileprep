// hkl-parser.js
// Parses SHELX-format HKL reflection data: h, k, l, F2, sigma[, batch]
//
// Standard fixed-width column format (FORTRAN spec): (3I4,2F8.2,I4)
//   columns:  h[1-4]  k[5-8]  l[9-12]  F2[13-20]  sigma[21-28]  batch[29-32]
// Terminator: a line "   0   0   0    0.00    0.00" marks end of data
// (SHELX itself ignores anything after this line).
//
// Some non-Bruker exports drift slightly from the fixed column grid
// (e.g. extra whitespace, different field widths for very large F2
// values). So: try fixed-width first, fall back to whitespace-separated
// parsing if that doesn't produce valid numbers.

function parseHKL(text) {
  const lines = text.split(/\r?\n/);
  const reflections = [];
  const skipped = [];

  for (let i = 0; i < lines.length; i++) {
    const rawLine = lines[i];
    if (rawLine.trim() === '') continue;

    let h, k, l, f2, sigma, batch = null;

    // 1) fixed-width attempt (SHELX specification)
    if (rawLine.length >= 28) {
      h = parseInt(rawLine.slice(0, 4), 10);
      k = parseInt(rawLine.slice(4, 8), 10);
      l = parseInt(rawLine.slice(8, 12), 10);
      f2 = parseFloat(rawLine.slice(12, 20));
      sigma = parseFloat(rawLine.slice(20, 28));
      const batchStr = rawLine.length >= 32 ? rawLine.slice(28, 32).trim() : '';
      batch = batchStr === '' ? null : parseInt(batchStr, 10);
    }

    // 2) fallback: whitespace-separated, in case fixed-width wasn't clean
    const fixedOk = [h, k, l, f2, sigma].every(Number.isFinite);
    if (!fixedOk) {
      const parts = rawLine.trim().split(/\s+/).map(Number);
      if (parts.length >= 5 && parts.slice(0, 5).every(Number.isFinite)) {
        [h, k, l, f2, sigma] = parts;
        batch = parts.length >= 6 && Number.isFinite(parts[5]) ? parts[5] : null;
      }
    }

    if (![h, k, l, f2, sigma].every(Number.isFinite)) {
      skipped.push({ lineNumber: i + 1, content: rawLine });
      continue;
    }

    // SHELX terminator: 0 0 0 0.00 0.00 (with sigma = 0)
    if (h === 0 && k === 0 && l === 0 && f2 === 0 && sigma === 0) break;

    reflections.push({ h, k, l, f2, sigma, batch });
  }

  return { reflections, skipped };
}

// Quick summary stats right after loading — useful as a first sanity
// check in the UI (reflection count, hkl range, F2 range).
function summarizeHKL(reflections) {
  if (reflections.length === 0) return { count: 0 };

  const hs = reflections.map(r => r.h);
  const ks = reflections.map(r => r.k);
  const ls = reflections.map(r => r.l);
  const f2s = reflections.map(r => r.f2);

  return {
    count: reflections.length,
    hRange: [Math.min(...hs), Math.max(...hs)],
    kRange: [Math.min(...ks), Math.max(...ks)],
    lRange: [Math.min(...ls), Math.max(...ls)],
    f2Range: [Math.min(...f2s), Math.max(...f2s)],
    meanF2: f2s.reduce((a, b) => a + b, 0) / f2s.length,
  };
}
