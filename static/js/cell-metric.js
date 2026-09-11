// cell-metric.js
// Determines the crystal system / metric symmetry (holohedry) compatible
// with a given unit cell, within tolerance, and derives which Laue classes
// are worth testing during reflection merging (laue-merge.js, next step).
//
// Important: this is METRIC symmetry only. A cell that looks orthorhombic
// by its lengths/angles doesn't guarantee the crystal actually has that
// symmetry — the true Laue class is only established by comparing merging
// statistics (Rint) for each candidate. This module just narrows down
// which candidates are geometrically possible, so the merging step doesn't
// waste time testing symmetry the cell metric can't even support.
//
// There are exactly 11 Laue classes in total (the centrosymmetric point
// groups) — this is a fixed crystallographic fact, not something derived
// here.
const ALL_LAUE_CLASSES = ['-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m'];

function anglesClose(x, y, tolDeg) {
  return Math.abs(x - y) < tolDeg;
}
function is90(angle, tolDeg) {
  return anglesClose(angle, 90, tolDeg);
}
function is120(angle, tolDeg) {
  return anglesClose(angle, 120, tolDeg);
}
function lengthsClose(x, y, relTol) {
  return Math.abs(x - y) / ((x + y) / 2) < relTol;
}

// Classifies the cell metric and returns the crystal system plus the
// ordered list of Laue classes to try during merging (highest symmetry
// first, down to -1). lengthTol is a relative tolerance on cell edges
// (default 0.5%), angleTol an absolute tolerance in degrees (default
// 0.5°) — loosen these for cells refined from lower-quality data.
function classifyCellMetric(cell, lengthTol = 0.005, angleTol = 0.5) {
  const { a, b, c, alpha, beta, gamma } = cell;

  const ab = lengthsClose(a, b, lengthTol);
  const bc = lengthsClose(b, c, lengthTol);
  const ac = lengthsClose(a, c, lengthTol);
  const allLengthsEqual = ab && bc && ac;

  const alphaBetaGamma90 = is90(alpha, angleTol) && is90(beta, angleTol) && is90(gamma, angleTol);
  const alphaGamma90 = is90(alpha, angleTol) && is90(gamma, angleTol);
  const alphaBeta90 = is90(alpha, angleTol) && is90(beta, angleTol);
  const allAnglesEqual = anglesClose(alpha, beta, angleTol) && anglesClose(beta, gamma, angleTol);

  // Cubic: a = b = c, all angles 90
  if (allLengthsEqual && alphaBetaGamma90) {
    return { system: 'cubic', laueClasses: ['m-3m', 'm-3', '4/mmm', 'mmm', '2/m', '-1'] };
  }

  // Rhombohedral setting: a = b = c, alpha = beta = gamma != 90
  if (allLengthsEqual && allAnglesEqual && !is90(alpha, angleTol)) {
    return { system: 'trigonal (rhombohedral setting)', laueClasses: ['-3m', '-3', '-1'] };
  }

  // Hexagonal metric: a = b, alpha = beta = 90, gamma = 120
  // (covers true hexagonal AND trigonal in hexagonal-axes setting — the
  // metric alone can't distinguish them, only absences/statistics can)
  if (ab && alphaBeta90 && is120(gamma, angleTol)) {
    return { system: 'hexagonal metric (hexagonal or trigonal)', laueClasses: ['6/mmm', '6/m', '-3m', '-3', '-1'] };
  }

  // Tetragonal: a = b != c, all angles 90
  if (ab && alphaBetaGamma90) {
    return { system: 'tetragonal', laueClasses: ['4/mmm', '4/m', 'mmm', '2/m', '-1'] };
  }

  // Orthorhombic: a != b != c, all angles 90
  if (alphaBetaGamma90) {
    return { system: 'orthorhombic', laueClasses: ['mmm', '2/m', '-1'] };
  }

  // Monoclinic, b-unique (the standard IUCr setting): alpha = gamma = 90, beta != 90
  if (alphaGamma90 && !is90(beta, angleTol)) {
    return { system: 'monoclinic (b-unique)', laueClasses: ['2/m', '-1'] };
  }

  // Monoclinic, c-unique (used by some software): alpha = beta = 90, gamma != 90
  if (alphaBeta90 && !is90(gamma, angleTol)) {
    return { system: 'monoclinic (c-unique)', laueClasses: ['2/m', '-1'] };
  }

  // No higher-symmetry metric constraint recognized — triclinic is always valid
  return { system: 'triclinic', laueClasses: ['-1'] };
}
