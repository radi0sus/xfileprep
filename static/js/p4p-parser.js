// p4p-parser.js
// Reads a Bruker/CrysAlisPro .p4p file (as exported by e.g. CrysAlisPro,
// SAINT) for the unit cell (+ esds), radiation wavelength, and chemical
// formula — the same handful of things this project already pulls from a
// CIF or a SHELX .ins, just under a different, line-keyword file format.
//
// .p4p line shapes this cares about (whitespace-separated, one instruction
// per line, keyword first):
//   CELL     a b c alpha beta gamma volume
//   CELLSD   sd(a) sd(b) sd(c) sd(alpha) sd(beta) sd(gamma) sd(volume)
//   CHEM     <element><count> <element><count> ...   (formula, one formula unit)
//   SOURCE   <element> avgWavelength Ka1 Ka2
//
// Deliberately narrow, same spirit as ins-parser.js: only what this project
// actually uses, not a full .p4p interpreter (ignores FACE/MAT/orientation
// matrix and everything else in the file).

// Requires looksLikeCIF (cif-reflections.js), looksLikeSHELXIns
// (ins-parser.js), and parseFormulaString (ins-writer.js) to be loaded first.
function looksLikeP4P(text) {
  if (looksLikeCIF(text)) return false;
  if (looksLikeSHELXIns(text)) return false;
  const hasCellLine = /^\s*CELL\s+(-?\d+\.?\d*\s+){6}-?\d+\.?\d*\s*$/im.test(text);
  const hasP4PMarker = /^\s*(TITLE|SOURCE|CHEM)\b/im.test(text);
  return hasCellLine && hasP4PMarker;
}

function findP4PLine(lines, keyword) {
  const re = new RegExp(`^\\s*${keyword}\\s+(.*)$`, 'i');
  for (const line of lines) {
    const m = line.match(re);
    if (m) return m[1].trim();
  }
  return null;
}

// A zero on the CELLSD line means "not actually determined" — several
// pipelines write an all-zero placeholder rather than omitting the line
// entirely (as seen in practice), so an all-zero set is treated the same
// as "no esds available" rather than as a genuine esd of 0.
function parseP4PFile(text) {
  const lines = text.split(/\r?\n/);

  const cellLine = findP4PLine(lines, 'CELL');
  if (!cellLine) return null;
  const cellParts = cellLine.split(/\s+/).map(Number);
  if (cellParts.length < 6 || cellParts.slice(0, 6).some(v => !Number.isFinite(v))) return null;
  const [a, b, c, alpha, beta, gamma, volume] = cellParts;

  let esd = null;
  const cellsdLine = findP4PLine(lines, 'CELLSD');
  if (cellsdLine) {
    const sdParts = cellsdLine.split(/\s+/).map(Number);
    if (sdParts.length >= 6 && sdParts.slice(0, 6).every(Number.isFinite)) {
      const allZero = sdParts.slice(0, 6).every(v => v === 0);
      if (!allZero) {
        esd = {
          a: sdParts[0], b: sdParts[1], c: sdParts[2],
          alpha: sdParts[3], beta: sdParts[4], gamma: sdParts[5],
        };
      }
    }
  }

  // SOURCE <element> avgWavelength Ka1 Ka2 — the first numeric token is
  // the (weighted average) wavelength, which is what the rest of this
  // project already uses elsewhere (CIF _diffrn_radiation_wavelength).
  let wavelength = null;
  const sourceLine = findP4PLine(lines, 'SOURCE');
  if (sourceLine) {
    const nums = sourceLine.split(/\s+/).map(Number).filter(Number.isFinite);
    if (nums.length > 0) wavelength = nums[0];
  }

  let formula = null;
  const chemLine = findP4PLine(lines, 'CHEM');
  if (chemLine) {
    const parsed = parseFormulaString(chemLine);
    if (parsed.length > 0) formula = parsed;
  }

  return {
    cell: { a, b, c, alpha, beta, gamma },
    esd,
    volume: Number.isFinite(volume) ? volume : null,
    wavelength,
    formula,
  };
}
