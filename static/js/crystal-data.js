// crystal-data.js
// Crystal data that follow from cell + composition + wavelength: Z estimate,
// density, F(000), linear absorption coefficient μ, μ·r and rough
// transmission limits, and Z′ for a chosen space group. Pure calculations —
// no DOM access; app.js does the wiring.
//
// Needs: absorption-data.js (ABSORPTION_DATA, μ/ρ and f″ from NIST FFAST),
// ins-writer.js (cellVolume), hall-symbol-parser.js (generateSpaceGroupOperators).

// Standard anode wavelengths (Kα, weighted mean of α1/α2) — the three for
// which ABSORPTION_DATA has values.
const ANODE_WAVELENGTHS = { Mo: 0.71073, Cu: 1.54184, Ag: 0.56087 };

// Standard atomic weights (IUPAC, abridged to 5 significant digits; for
// elements without stable isotopes the mass number of the longest-lived
// isotope). Cross-checked against the atomic weights implied by the FFAST
// data (A = 2 r_e λ N_A f2 / (μ/ρ)_pe): all 92 agree to better than 0.3 %.
// Element number = position + 1.
const ELEMENT_SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn ' +
  'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd ' +
  'Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U').split(' ');
const ATOMIC_WEIGHTS = [1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990,
  24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996,
  54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904, 83.798, 85.468,
  87.62, 88.906, 91.224, 92.906, 95.95, 98, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, 121.76,
  127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, 145, 150.36, 151.96, 157.25,
  158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22,
  195.08, 196.97, 200.59, 204.38, 207.2, 208.98, 209, 210, 222, 223, 226, 227, 232.04, 231.04, 238.03];

// SHELX allows D for deuterium — same electrons and (per atom) the same
// absorption cross-section as H, but twice the mass.
function elementInfo(symbol) {
  if (symbol === 'D') return { symbol: 'D', z: 1, mass: 2.0141, abs: ABSORPTION_DATA.H };
  const i = ELEMENT_SYMBOLS.indexOf(symbol);
  if (i < 0) return null;
  return { symbol, z: i + 1, mass: ATOMIC_WEIGHTS[i], abs: ABSORPTION_DATA[symbol] || null };
}

// Which standard anode a wavelength corresponds to (±0.002 Å covers the
// usual roundings: 0.7107, 1.5418, 1.54178, 0.5609 …), or null.
function anodeForWavelength(lambda) {
  if (!Number.isFinite(lambda)) return null;
  for (const [name, l] of Object.entries(ANODE_WAVELENGTHS)) {
    if (Math.abs(lambda - l) <= 0.002) return name;
  }
  return null;
}

// ~18 Å³ per non-hydrogen atom (small-molecule rule of thumb). Returns the
// raw quotient as well as the rounded Z, so the UI can show how far off an
// integer it is.
function estimateZDetails(volume, formula) {
  const nonH = (formula || []).filter(f => f.element !== 'H' && f.element !== 'D').reduce((s, f) => s + f.count, 0);
  if (!(volume > 0) || nonH === 0) return null;
  const atomsPerCell = volume / 18;
  const zRaw = atomsPerCell / nonH;
  return { atomsPerCell, nonHPerUnit: nonH, zRaw, z: Math.max(1, Math.round(zRaw)) };
}

// Number of general positions of a space group (incl. centring), for Z′.
function spaceGroupMultiplicity(hallSymbol) {
  try { return generateSpaceGroupOperators(hallSymbol).length; } catch (e) { return null; }
}

// Z′ = Z / multiplicity, with a short note for the UI.
function zPrimeInfo(z, hallSymbol) {
  const mult = spaceGroupMultiplicity(hallSymbol);
  if (!z || !mult) return null;
  const zp = z / mult;
  let note = '';
  if (zp < 1 - 1e-9) note = 'Z′ < 1: the molecule (or formula unit) must sit on a special position (inversion centre, mirror, 2-fold axis …), or Z is wrong.';
  else if (zp > 1 + 1e-9) note = `Z′ > 1: ${Number.isInteger(zp) ? zp : zp.toFixed(2)} independent formula units in the asymmetric unit — possible, but also a hint at missed symmetry (check with ADDSYM).`;
  return { mult, zPrime: zp, note };
}

// Main calculation.
//   volume       cell volume (Å³)
//   composition  whole-cell element counts [{element, count}]
//   z, formula   (optional) Z and per-formula-unit formula, only for the formula weight
//   wavelength   Å (null → no μ)
//   size         [d1, d2, d3] crystal dimensions in mm, or null
function computeCrystalData({ volume, composition, z, formula, wavelength, size }) {
  if (!(volume > 0) || !composition || composition.length === 0) return null;
  const unknown = [];
  let massCell = 0, electrons = 0;
  const anode = anodeForWavelength(wavelength);
  // f″ as SHELXL uses it: International Tables Vol. C Table 4.2.6.8 at the
  // ten tabulated lines (dispersion-table.js); NIST FFAST only as fallback.
  const ita = typeof dispersionForWavelength === 'function' && wavelength ? dispersionForWavelength(wavelength) : null;
  let muSum = 0; // Σ n·A·(μ/ρ)
  let muComplete = !!anode;
  const perElement = [];
  for (const c of composition) {
    const info = elementInfo(c.element);
    if (!info) { unknown.push(c.element); continue; }
    massCell += c.count * info.mass;
    electrons += c.count * info.z;
    const ab = anode && info.abs ? info.abs[anode] : null;
    if (anode && !ab) muComplete = false;
    if (ab) muSum += c.count * info.mass * ab[0];
    const fppIta = ita && ita.fpp[c.element] !== undefined ? ita.fpp[c.element] : null;
    perElement.push({ element: c.element, count: c.count, muRho: ab ? ab[0] : null,
      fpp: fppIta !== null ? fppIta : (ab ? ab[1] : null), fppSource: fppIta !== null ? 'ITA' : (ab ? 'FFAST' : null) });
  }
  if (unknown.length > 0) return { unknown, volume };

  // ρ = M_cell / (N_A V);  1 / (N_A · 1e-24 cm³/Å³) = 1.66054
  const density = massCell * 1.66054 / volume;
  // μ = Σ n A (μ/ρ) / (N_A V)  in cm⁻¹ → mm⁻¹
  const mu = muComplete ? muSum / (0.602214076 * volume) / 10 : null;

  let fw = null;
  if (formula && formula.length > 0) {
    fw = 0;
    for (const f of formula) {
      const info = elementInfo(f.element);
      if (!info) { fw = null; break; }
      fw += f.count * info.mass;
    }
  }

  let absorption = null;
  const dims = (size || []).filter(d => Number.isFinite(d) && d > 0);
  if (mu !== null && dims.length === 3) {
    const dMin = Math.min(...dims), dMax = Math.max(...dims);
    const muRMax = mu * dMax / 2, muRMin = mu * dMin / 2;
    // Crude plate model, as used for "expected" Tmin/Tmax in checkCIF-style
    // checks: T ≈ exp(−μ d) along the longest / shortest dimension.
    const level = muRMax < 0.1 ? 'negligible' : muRMax < 1 ? 'moderate' : muRMax < 3 ? 'strong' : 'very strong';
    absorption = { dMin, dMax, muRMin, muRMax, tMin: Math.exp(-mu * dMax), tMax: Math.exp(-mu * dMin), level };
  }

  return { volume, massCell, formulaWeight: fw, density, f000: electrons, anode, fppLine: ita ? ita.line : (anode ? `${anode} Kα` : null), mu, perElement, absorption };
}
