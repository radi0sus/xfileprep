// absences.js
// Scans the observed reflection data for systematic absences: lattice
// centering (whole dataset), axial screw-axis conditions (h00/0k0/00l),
// and zonal glide-plane conditions (0kl/h0l/hk0).
//
// (Zone/axial visualization for this data lives in absences-viz.js,
// which reuses ZONE_SELECTORS, AXIAL_CLASSES, and isWeak from here.)
//
// This module only LOOKS at the data — it doesn't assume any particular
// space group or use a lookup table of known conditions. For each
// candidate condition, we partition reflections by residue class and
// check whether the "forbidden" residues are statistically absent while
// the "allowed" residue is normally observed. This is exactly the
// evidence XPREP/PLATON use to build an extinction symbol, and it's what
// later gets matched against space-group candidates (via the symmetry
// operators from the Hall-symbol parser, not a second absences table).
//
// Deliberately works on unmerged reflections — there's already enough
// redundancy in typical single-crystal data for this statistical test,
// and it avoids a dependency on the merging step.

// A reflection counts as "weak/absent" if I/sigma is below this factor.
// Judging by the whole residue class (not single reflections) matters
// because real data always has some noisy weak reflections even where
// no absence applies.
const ABSENCE_SIGNIFICANCE_THRESHOLD = 3; // I / sigma(I)

function isWeak(f2, sigma) {
  if (sigma <= 0) return true;
  return f2 / sigma < ABSENCE_SIGNIFICANCE_THRESHOLD;
}

function mod(n, m) {
  return ((n % m) + m) % m;
}

// Partitions reflections by valueOf(r) mod modulus and computes basic
// statistics per residue class.
function residueBuckets(reflections, valueOf, modulus) {
  const buckets = Array.from({ length: modulus }, () => []);
  for (const r of reflections) buckets[mod(valueOf(r), modulus)].push(r);

  return buckets.map((group, residue) => {
    const count = group.length;
    const weakCount = group.filter(r => isWeak(r.f2, r.sigma)).length;
    return {
      residue,
      count,
      fractionWeak: count > 0 ? weakCount / count : null,
      meanF2: count > 0 ? group.reduce((a, r) => a + r.f2, 0) / count : null,
    };
  });
}

// Judges whether a condition holds: residue 0 is the "allowed" class by
// convention (index ≡ 0 mod n is always the non-forbidden residue for
// every condition used here), every other residue is "forbidden" and
// should be statistically absent.
// Judges whether a condition holds and reports the underlying percentages,
// so ambiguous cases aren't just a flat "inconclusive" label — the person
// can see the actual numbers XPREP-style and decide for themselves.
// Residue 0 is the "allowed" class by convention; every other residue is
// "forbidden" and should be statistically absent.
function judgeCondition(buckets, { minCount = 5, violationPercentThreshold = 15, presentPercentThreshold = 50 } = {}) {
  const allowed = buckets[0];
  const forbidden = buckets.slice(1);

  let forbiddenCount = 0;
  let forbiddenPresentCount = 0; // "present" = NOT weak, i.e. a real violation of the condition
  let forbiddenF2Sum = 0;
  for (const b of forbidden) {
    if (b.count === 0) continue; // never even generated — itself supports the condition, contributes 0 to violation%
    const presentCount = Math.round((1 - b.fractionWeak) * b.count);
    forbiddenCount += b.count;
    forbiddenPresentCount += presentCount;
    forbiddenF2Sum += b.meanF2 * b.count;
  }

  const allowedPresentPercent = allowed.count > 0 ? (1 - allowed.fractionWeak) * 100 : null;
  const violationPercent = forbiddenCount > 0 ? (forbiddenPresentCount / forbiddenCount) * 100 : 0;
  // mean F2 of the forbidden class as % of the allowed class's mean F2 —
  // near 0% is the clean case, high % suggests the class isn't actually absent
  const meanRatioPercent = forbiddenCount > 0 && allowed.meanF2 > 0
    ? ((forbiddenF2Sum / forbiddenCount) / allowed.meanF2) * 100
    : 0;

  let verdict;
  let reason = null;
  if (allowed.count < minCount) {
    verdict = 'inconclusive';
    reason = `too few reference reflections in the allowed class (n=${allowed.count}, need ≥${minCount})`;
  } else if (allowedPresentPercent < (100 - presentPercentThreshold)) {
    verdict = 'inconclusive';
    reason = `even the allowed class is mostly weak here (only ${allowedPresentPercent.toFixed(0)}% present) — nothing to compare the forbidden class against`;
  } else if (violationPercent > violationPercentThreshold) {
    verdict = 'violated';
  } else {
    verdict = 'holds';
  }

  return { verdict, reason, violationPercent, meanRatioPercent, forbiddenCount, allowedCount: allowed.count };
}

// --- Lattice centering (whole dataset) ---------------------------------

const CENTERING_CANDIDATES = {
  'P (primitive)': [],
  'A': [{ valueOf: r => r.k + r.l, modulus: 2 }],
  'B': [{ valueOf: r => r.h + r.l, modulus: 2 }],
  'C': [{ valueOf: r => r.h + r.k, modulus: 2 }],
  'I': [{ valueOf: r => r.h + r.k + r.l, modulus: 2 }],
  'F': [
    { valueOf: r => r.h + r.k, modulus: 2 },
    { valueOf: r => r.h + r.l, modulus: 2 },
    { valueOf: r => r.k + r.l, modulus: 2 },
  ],
  'R (obverse)': [{ valueOf: r => -r.h + r.k + r.l, modulus: 3 }],
  'R (reverse)': [{ valueOf: r => r.h - r.k - r.l, modulus: 3 }],
};

function scanLatticeCentering(reflections) {
  const results = {};
  for (const [name, conditions] of Object.entries(CENTERING_CANDIDATES)) {
    if (conditions.length === 0) {
      results[name] = { verdict: 'holds', violationPercent: 0, note: 'no restriction — always compatible', forbiddenCount: 0, allowedCount: reflections.length };
      continue;
    }
    const judged = conditions.map(c => judgeCondition(residueBuckets(reflections, c.valueOf, c.modulus)));
    // for a combined centering (F needs all three sub-conditions), report the
    // worst-case violation percentage — that's what would actually break the rule
    const worst = judged.reduce((a, b) => (b.violationPercent > a.violationPercent ? b : a));
    const overall = judged.some(j => j.verdict === 'violated') ? 'violated'
      : judged.every(j => j.verdict === 'holds') ? 'holds'
      : 'inconclusive';
    results[name] = { verdict: overall, violationPercent: worst.violationPercent, forbiddenCount: worst.forbiddenCount, allowedCount: worst.allowedCount };
  }
  return results;
}

// --- Axial screw-axis conditions (h00, 0k0, 00l) -----------------------

const AXIAL_CLASSES = {
  h00: { select: r => r.k === 0 && r.l === 0, indexOf: r => r.h },
  '0k0': { select: r => r.h === 0 && r.l === 0, indexOf: r => r.k },
  '00l': { select: r => r.h === 0 && r.k === 0, indexOf: r => r.l },
};

function scanAxialCondition(reflections, axisName, periods = [2, 3, 4, 6]) {
  const cls = AXIAL_CLASSES[axisName];
  const relevant = reflections.filter(cls.select);
  return periods.map(n => ({
    period: n,
    ...judgeCondition(residueBuckets(relevant, cls.indexOf, n)),
    reflectionCount: relevant.length,
  }));
}

// --- Zonal glide-plane conditions (0kl, h0l, hk0) -----------------------
// d-glide (quarter-translation) isn't included — it needs an additional
// same-parity check beyond a single mod-2/mod-4 condition and is rare
// enough to leave as a future extension.

const ZONE_SELECTORS = {
  '0kl': r => r.h === 0,
  h0l: r => r.k === 0,
  hk0: r => r.l === 0,
};

const ZONAL_CONDITIONS = {
  '0kl': [
    { name: 'b-glide (k=2n)', valueOf: r => r.k },
    { name: 'c-glide (l=2n)', valueOf: r => r.l },
    { name: 'n-glide (k+l=2n)', valueOf: r => r.k + r.l },
  ],
  h0l: [
    { name: 'a-glide (h=2n)', valueOf: r => r.h },
    { name: 'c-glide (l=2n)', valueOf: r => r.l },
    { name: 'n-glide (h+l=2n)', valueOf: r => r.h + r.l },
  ],
  hk0: [
    { name: 'a-glide (h=2n)', valueOf: r => r.h },
    { name: 'b-glide (k=2n)', valueOf: r => r.k },
    { name: 'n-glide (h+k=2n)', valueOf: r => r.h + r.k },
  ],
};

function scanZonalConditions(reflections, zoneName) {
  const relevant = reflections.filter(ZONE_SELECTORS[zoneName]);
  return ZONAL_CONDITIONS[zoneName].map(cond => ({
    name: cond.name,
    ...judgeCondition(residueBuckets(relevant, cond.valueOf, 2)),
    reflectionCount: relevant.length,
  }));
}

// --- Top-level scan ------------------------------------------------------

function scanAllAbsences(reflections) {
  return {
    centering: scanLatticeCentering(reflections),
    axial: {
      h00: scanAxialCondition(reflections, 'h00'),
      '0k0': scanAxialCondition(reflections, '0k0'),
      '00l': scanAxialCondition(reflections, '00l'),
    },
    zonal: {
      '0kl': scanZonalConditions(reflections, '0kl'),
      h0l: scanZonalConditions(reflections, 'h0l'),
      hk0: scanZonalConditions(reflections, 'hk0'),
    },
  };
}
