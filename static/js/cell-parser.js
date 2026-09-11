// cell-parser.js
// Extracts unit cell constants and related metadata from a CIF file.
//
// Unlike the reflection tags (which vary across software), the core cell
// tags are part of the standard CIF dictionary and are consistently named,
// so no alias handling is needed here.

// Finds the value for a single (non-loop) CIF tag, e.g. "_cell_length_a".
// Handles both "_tag value" on one line and "_tag" followed by the value
// on the next non-blank line (used for long values). Returns the raw
// string value, or null if the tag isn't found.
function extractCIFTagValue(cifText, tagName) {
  const lines = cifText.split(/\r?\n/);
  const tagRegex = new RegExp(`^\\s*${tagName}\\s*(.*)$`, 'i');

  for (let i = 0; i < lines.length; i++) {
    const match = lines[i].match(tagRegex);
    if (!match) continue;

    const sameLine = match[1].trim();
    if (sameLine !== '') return stripQuotes(sameLine);

    let j = i + 1;
    while (j < lines.length && lines[j].trim() === '') j++;
    if (j < lines.length && !/^\s*(_|loop_|data_)/i.test(lines[j])) {
      return stripQuotes(lines[j].trim());
    }
    return null;
  }
  return null;
}

function stripQuotes(value) {
  if ((value.startsWith("'") && value.endsWith("'")) ||
      (value.startsWith('"') && value.endsWith('"'))) {
    return value.slice(1, -1);
  }
  return value;
}

// CIF numeric values often carry a standard uncertainty in parentheses,
// e.g. "12.3456(7)" — the esd is in units of the last printed digit, so
// "6.1455(3)" means esd = 0.0003, while "1132.65(12)" means esd = 0.12.
// Returns { value, esd } — esd is null when the source has no parenthesis
// (or the value itself couldn't be parsed at all).
function parseCIFNumberWithEsd(value) {
  if (value === null || value === undefined) return { value: NaN, esd: null };
  const match = value.match(/^(-?\d+\.?\d*)\((\d+)\)$/);
  if (!match) return { value: parseFloat(value), esd: null };

  const numStr = match[1];
  const decimalIndex = numStr.indexOf('.');
  const decimalPlaces = decimalIndex === -1 ? 0 : numStr.length - decimalIndex - 1;
  const esd = parseInt(match[2], 10) * Math.pow(10, -decimalPlaces);
  return { value: parseFloat(numStr), esd };
}

// CIF numeric values often carry a standard uncertainty in parentheses,
// e.g. "12.3456(7)" — strip that before parsing as a float.
function parseCIFNumber(value) {
  return parseCIFNumberWithEsd(value).value;
}

// Extracts unit cell constants (plus a few useful extras) from CIF text.
// Returns null if the core six cell parameters aren't all present —
// callers should fall back to manual entry in that case.
function extractCellFromCIF(cifText) {
  const aInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_length_a'));
  const bInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_length_b'));
  const cInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_length_c'));
  const alphaInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_angle_alpha'));
  const betaInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_angle_beta'));
  const gammaInfo = parseCIFNumberWithEsd(extractCIFTagValue(cifText, '_cell_angle_gamma'));

  const a = aInfo.value, b = bInfo.value, c = cInfo.value;
  const alpha = alphaInfo.value, beta = betaInfo.value, gamma = gammaInfo.value;

  if (![a, b, c, alpha, beta, gamma].every(Number.isFinite)) return null;

  const volume = parseCIFNumber(extractCIFTagValue(cifText, '_cell_volume'));
  const z = parseCIFNumber(extractCIFTagValue(cifText, '_cell_formula_units_Z'));
  const wavelength = parseCIFNumber(extractCIFTagValue(cifText, '_diffrn_radiation_wavelength'));

  // reported space group is only for cross-checking later against our own
  // determination — two common tag names are in use depending on CIF version
  const reportedSpaceGroup =
    extractCIFTagValue(cifText, '_symmetry_space_group_name_H-M') ||
    extractCIFTagValue(cifText, '_space_group_name_H-M_alt');

  const esdValues = [aInfo.esd, bInfo.esd, cInfo.esd, alphaInfo.esd, betaInfo.esd, gammaInfo.esd];
  const esd = esdValues.some(v => v !== null)
    ? { a: aInfo.esd, b: bInfo.esd, c: cInfo.esd, alpha: alphaInfo.esd, beta: betaInfo.esd, gamma: gammaInfo.esd }
    : null;

  return {
    a, b, c, alpha, beta, gamma,
    volume: Number.isFinite(volume) ? volume : null,
    z: Number.isFinite(z) ? z : null,
    wavelength: Number.isFinite(wavelength) ? wavelength : null,
    reportedSpaceGroup: reportedSpaceGroup || null,
    esd,
  };
}
