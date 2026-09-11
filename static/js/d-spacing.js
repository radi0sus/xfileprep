// d-spacing.js
// Resolution (d-spacing / sin(theta)/lambda) from Miller indices and the
// unit cell — general triclinic formula, so it works unmodified for every
// crystal system rather than needing per-system special cases.
//
// No wavelength is needed anywhere here: sin(theta)/lambda = 1/(2d), and
// d itself only depends on h,k,l and the (reciprocal) cell metric — the
// wavelength cancels out of Bragg's law when expressed this way. That's
// exactly what Wilson/E-statistics need, since they characterize the
// crystal's scattering, not the particular experiment's radiation.

const DEG_TO_RAD = Math.PI / 180;

// Standard crystallographic reciprocal-cell formulas (e.g. Giacovazzo,
// "Fundamentals of Crystallography", eq. 2.20-2.23).
function reciprocalCellFromDirect(cell) {
  const alpha = cell.alpha * DEG_TO_RAD;
  const beta = cell.beta * DEG_TO_RAD;
  const gamma = cell.gamma * DEG_TO_RAD;
  const { a, b, c } = cell;

  const cosA = Math.cos(alpha), cosB = Math.cos(beta), cosG = Math.cos(gamma);
  const sinA = Math.sin(alpha), sinB = Math.sin(beta), sinG = Math.sin(gamma);

  const volume = a * b * c * Math.sqrt(
    Math.max(0, 1 - cosA * cosA - cosB * cosB - cosG * cosG + 2 * cosA * cosB * cosG)
  );

  const aStar = (b * c * sinA) / volume;
  const bStar = (c * a * sinB) / volume;
  const cStar = (a * b * sinG) / volume;

  const cosAStar = (cosB * cosG - cosA) / (sinB * sinG);
  const cosBStar = (cosG * cosA - cosB) / (sinG * sinA);
  const cosGStar = (cosA * cosB - cosG) / (sinA * sinB);

  return { aStar, bStar, cStar, cosAStar, cosBStar, cosGStar, volume };
}

// 1/d^2 for a given hkl, from the general triclinic reciprocal-metric
// formula — reduces correctly to the simpler monoclinic/orthorhombic/
// cubic formulas as a special case, no separate code path needed per
// crystal system.
function inverseDSquared(h, k, l, recip) {
  const { aStar, bStar, cStar, cosAStar, cosBStar, cosGStar } = recip;
  return (
    h * h * aStar * aStar +
    k * k * bStar * bStar +
    l * l * cStar * cStar +
    2 * k * l * bStar * cStar * cosAStar +
    2 * l * h * cStar * aStar * cosBStar +
    2 * h * k * aStar * bStar * cosGStar
  );
}

// sin(theta)/lambda in 1/Angstrom, i.e. 1/(2d).
function sinThetaOverLambda(h, k, l, recip) {
  return Math.sqrt(Math.max(0, inverseDSquared(h, k, l, recip))) / 2;
}
