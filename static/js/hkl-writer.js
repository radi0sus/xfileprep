// hkl-writer.js
// Writes reflections back out as a SHELX-format HKL file, FORTRAN (3I4,2F8.2,I4),
// terminated by the usual "   0   0   0    0.00    0.00" line.
//
// Deliberately NO transformation, merging or filtering: indices, F2 and sigma
// are written exactly as parsed. Any setting change is done on the .ins side
// via the HKLF matrix ("HKLF 4 1 ..."), so this file stays the raw data.
//
// Only meant for reflections that came out of a CIF (a _refln loop). A raw
// SHELX .hkl dropped by the user is already in that format and is not re-saved.

// Formats a float into a fixed width field. Starts at 2 decimals and drops
// decimals only if a very large value would otherwise overflow the field
// (never rescales, never changes the value beyond rounding).
function formatHKLField(value, width) {
  for (let decimals = 2; decimals >= 0; decimals--) {
    const s = value.toFixed(decimals);
    if (s.length <= width) return s.padStart(width);
  }
  return value.toExponential(0).padStart(width); // practically unreachable
}

function writeSHELXHKL(reflections) {
  const out = [];
  for (const r of reflections) {
    const batch = Number.isFinite(r.batch) ? r.batch : 1;
    out.push(
      String(r.h).padStart(4) + String(r.k).padStart(4) + String(r.l).padStart(4) +
      formatHKLField(r.f2, 8) + formatHKLField(r.sigma, 8) + String(batch).padStart(4)
    );
  }
  out.push('   0   0   0    0.00    0.00');
  return out.join('\n') + '\n';
}
