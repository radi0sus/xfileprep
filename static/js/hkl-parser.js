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

  // Plain loop instead of Math.min(...array): spreading a six-figure array
  // into a call overflows the call stack (JavaScriptCore on iPad/iPhone
  // already fails at ~65k arguments, V8 at ~120k), and real data sets
  // easily hold 140 000 reflections.
  const lo = { h: Infinity, k: Infinity, l: Infinity, f2: Infinity };
  const hi = { h: -Infinity, k: -Infinity, l: -Infinity, f2: -Infinity };
  let sumF2 = 0;
  for (const r of reflections) {
    if (r.h < lo.h) lo.h = r.h; if (r.h > hi.h) hi.h = r.h;
    if (r.k < lo.k) lo.k = r.k; if (r.k > hi.k) hi.k = r.k;
    if (r.l < lo.l) lo.l = r.l; if (r.l > hi.l) hi.l = r.l;
    if (r.f2 < lo.f2) lo.f2 = r.f2; if (r.f2 > hi.f2) hi.f2 = r.f2;
    sumF2 += r.f2;
  }

  return {
    count: reflections.length,
    hRange: [lo.h, hi.h],
    kRange: [lo.k, hi.k],
    lRange: [lo.l, hi.l],
    f2Range: [lo.f2, hi.f2],
    meanF2: sumF2 / reflections.length,
  };
}
