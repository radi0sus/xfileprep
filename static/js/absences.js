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

// Pure "fraction of reflections with I/sigma >= 3" is fragile on very strong
// data: sigmas are tiny compared with the intensity scale, so even a residual
// signal of ~0.1% of an allowed reflection (Renninger/lambda/2 leakage, tails of
// neighbouring strong reflections) clears 3 sigma for a noticeable share of the
// forbidden class. A class whose MEDIAN I/sigma is below the significance
// threshold is therefore treated as absent, unless more than
// MAX_TAIL_PERCENT of its reflections are individually significant.
const MAX_TAIL_PERCENT = 50;

// The median test above compares the forbidden class with an ABSOLUTE
// threshold (I/σ < 3). That is not enough on weakly diffracting data: when
// even the allowed class has a median I/σ of only ~4, a forbidden class that
// is really present (pseudo-centering, e.g. weak h+l odd reflections of a
// primitive cell with a pseudo-B translation) also falls below 3. So the
// weak-tail rule only applies if the forbidden class is weak RELATIVE to the
// allowed one: its mean intensity must stay below this percentage of the
// allowed class's mean. Truly absent classes sit at ≲1 % (Renninger and
// neighbour leakage); a class at tens of percent is simply weak data.
const MAX_TAIL_MEAN_RATIO_PERCENT = 10;

function median(values) {
  if (values.length === 0) return null;
  const a = [...values].sort((x, y) => x - y);
  const n = a.length;
  return n % 2 ? a[(n - 1) / 2] : (a[n / 2 - 1] + a[n / 2]) / 2;
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
      iSigmas: group.map(r => (r.sigma > 0 ? r.f2 / r.sigma : 0)),
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

  const medianISigma = median(forbidden.flatMap(b => b.iSigmas));

  let verdict;
  let reason = null;
  let weakTail = false; // "holds" only thanks to the weak-residual-tail rule
  if (allowed.count < minCount && allowed.count >= 1 && forbiddenCount >= 2
      && violationPercent > violationPercentThreshold && meanRatioPercent >= MAX_TAIL_MEAN_RATIO_PERCENT
      && (medianISigma === null || medianISigma >= ABSENCE_SIGNIFICANCE_THRESHOLD)) {
    // Too few allowed reflections to CONFIRM a condition — but a violation needs
    // no large reference: forbidden reflections that are clearly significant and
    // comparable in strength to the allowed ones refute it. (SrCu2(BO3)2, COD
    // 1501632, I -4 2 m: 00l has only 004/008 for l = 4n, while 002 and 006 are
    // strong — the 4_1 of I 4_1 2 2 is clearly absent, not "undetermined".)
    verdict = 'violated';
    reason = `only ${allowed.count} allowed reflection(s), but ${violationPercent.toFixed(0)}% of the ${forbiddenCount} forbidden ones exceed ${ABSENCE_SIGNIFICANCE_THRESHOLD}σ `
      + `with ${meanRatioPercent.toFixed(0)}% of the allowed mean intensity — clearly present`;
  } else if (allowed.count < minCount) {
    verdict = 'inconclusive';
    reason = `too few reference reflections in the allowed class (n=${allowed.count}, need ≥${minCount})`;
  } else if (allowedPresentPercent < (100 - presentPercentThreshold)) {
    verdict = 'inconclusive';
    reason = `even the allowed class is mostly weak here (only ${allowedPresentPercent.toFixed(0)}% present) — nothing to compare the forbidden class against`;
  } else if (violationPercent > violationPercentThreshold) {
    if (medianISigma !== null && medianISigma < ABSENCE_SIGNIFICANCE_THRESHOLD && violationPercent <= MAX_TAIL_PERCENT
        && meanRatioPercent < MAX_TAIL_MEAN_RATIO_PERCENT) {
      // typical forbidden reflection is insignificant; the significant ones are
      // the upper tail of residual intensity, not a real violation
      verdict = 'holds';
      weakTail = true;
      reason = `${violationPercent.toFixed(0)}% of single forbidden reflections exceed ${ABSENCE_SIGNIFICANCE_THRESHOLD}σ, but the class median is only I/σ = ${medianISigma.toFixed(1)} `
        + `(forbidden/allowed mean intensity ${meanRatioPercent.toFixed(2)}%) — a weak residual tail, not a violation`;
    } else {
      verdict = 'violated';
      reason = `${violationPercent.toFixed(0)}% of the forbidden reflections exceed ${ABSENCE_SIGNIFICANCE_THRESHOLD}σ and their mean intensity is ${meanRatioPercent.toFixed(1)}% of the allowed class — present, not absent`;
    }
  } else {
    verdict = 'holds';
  }

  return { verdict, reason, weakTail, violationPercent, meanRatioPercent, medianISigma, forbiddenCount, forbiddenPresentCount, allowedCount: allowed.count, allowedPresentPercent };
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
    results[name] = { verdict: overall, violationPercent: worst.violationPercent, forbiddenCount: worst.forbiddenCount,
      forbiddenPresentCount: worst.forbiddenPresentCount, allowedCount: worst.allowedCount,
      reason: worst.reason, weakTail: overall === 'holds' && judged.some(j => j.weakTail) };
  }
  return results;
}

// The lattice the data actually show. Several centerings can "hold" at once
// — an F lattice satisfies A, B and C as well (F = A + B + C), so taking the
// first held entry in table order picked "A" for every F-centred dataset and
// left no candidates in m-3m (COD 2023426, F m -3 m). Most specific first:
// F, then I, R, then a single face (two faces holding imply the third = F).
const CENTERING_SPECIFICITY = ['F', 'I', 'R (obverse)', 'R (reverse)', 'A', 'B', 'C'];
function observedCenteringName(centering) {
  const holds = name => centering[name] && centering[name].verdict === 'holds';
  const faces = ['A', 'B', 'C'].filter(holds);
  if (faces.length >= 2 && centering.F && centering.F.verdict !== 'violated') return 'F';
  return CENTERING_SPECIFICITY.find(holds) || 'P (primitive)';
}

// --- Axial screw-axis conditions (h00, 0k0, 00l) -----------------------

const AXIAL_CLASSES = {
  h00: { select: r => r.k === 0 && r.l === 0, indexOf: r => r.h },
  '0k0': { select: r => r.h === 0 && r.l === 0, indexOf: r => r.k },
  '00l': { select: r => r.h === 0 && r.k === 0, indexOf: r => r.l },
};

// Axial rows hold only a handful of reflections by nature: along a short
// axis (b = 5.5 A at d_min = 0.77 A: 0k0 runs to k = 7) there are just
// three allowed ones, 020/040/060. The zonal default of 5 reference
// reflections made a 2_1 along such an axis permanently "undetermined", so
// P 2/a, P a and P 21/a tied (COD 1543504, actually P 21/a). XPREP judges
// screw axes from as few reflections; 3 is the minimum used here — the
// allowed ones still have to be mostly significant (presentPercentThreshold).
const AXIAL_MIN_ALLOWED_COUNT = 3;

function scanAxialCondition(reflections, axisName, periods = [2, 3, 4, 6]) {
  const cls = AXIAL_CLASSES[axisName];
  const relevant = reflections.filter(cls.select);
  return periods.map(n => {
    // `given[p]`: the same period-n test, but only among the reflections that
    // period p (a divisor of n) already allows — "of the k = 2n reflections, are
    // the k = 4n+2 ones absent?". That is the question a candidate predicting
    // period p must answer for a stricter period n. Asked over the whole row,
    // the absences the candidate itself predicts (k odd) dilute the violation
    // percentage and cost the RIGHT candidate credit (P 21/n lost to P n on
    // 0k0 k = 6n: 42 % violation overall, but 100 % among the even k).
    const given = {};
    for (const p of periods) {
      if (p >= n || n % p !== 0) continue;
      const subset = relevant.filter(r => mod(cls.indexOf(r), p) === 0);
      given[p] = judgeCondition(residueBuckets(subset, r => cls.indexOf(r) / p, n / p), { minCount: AXIAL_MIN_ALLOWED_COUNT });
    }
    return {
      period: n,
      ...judgeCondition(residueBuckets(relevant, cls.indexOf, n), { minCount: AXIAL_MIN_ALLOWED_COUNT }),
      reflectionCount: relevant.length,
      given,
    };
  });
}

// --- Zonal glide-plane conditions (0kl, h0l, hk0 and the diagonal plane hhl) --
// Each condition is  a·p + b·q ≡ 0 (mod m)  on the two free indices (p, q) of
// the zone, stored as `test: [a, b, m]` so the theoretical side
// (reflection-conditions.js) can compare a candidate against exactly the same
// tests. Free indices: 0kl → (k, l); h0l → (h, l); hk0 → (h, k); hhl → (h, l).
//   mod 2 = a/b/c/n glides, mod 4 on a SUM = d-glide (quarter translation).
// hhl is the plane through the diagonal (k = h): the mirror/glide planes ⊥ [1-10]
// of tetragonal, ⊥ [210] of hexagonal/trigonal and the {110} planes of cubic
// groups. Under 4-fold (or 3/6-fold) symmetry its partner h-hl is equivalent, so
// one plane is enough; in other crystal systems nothing special happens there
// (only the centering shows up), which affects all candidates alike.

const ZONE_SELECTORS = {
  '0kl': r => r.h === 0,
  h0l: r => r.k === 0,
  hk0: r => r.l === 0,
  hhl: r => r.h === r.k,
};

const ZONAL_CONDITIONS = {
  '0kl': [
    { name: 'b-glide (k=2n)', valueOf: r => r.k, test: [1, 0, 2] },
    { name: 'c-glide (l=2n)', valueOf: r => r.l, test: [0, 1, 2] },
    { name: 'n-glide (k+l=2n)', valueOf: r => r.k + r.l, test: [1, 1, 2] },
    { name: 'd-glide (k+l=4n)', valueOf: r => r.k + r.l, modulus: 4, test: [1, 1, 4] },
  ],
  h0l: [
    { name: 'a-glide (h=2n)', valueOf: r => r.h, test: [1, 0, 2] },
    { name: 'c-glide (l=2n)', valueOf: r => r.l, test: [0, 1, 2] },
    { name: 'n-glide (h+l=2n)', valueOf: r => r.h + r.l, test: [1, 1, 2] },
    { name: 'd-glide (h+l=4n)', valueOf: r => r.h + r.l, modulus: 4, test: [1, 1, 4] },
  ],
  hk0: [
    { name: 'a-glide (h=2n)', valueOf: r => r.h, test: [1, 0, 2] },
    { name: 'b-glide (k=2n)', valueOf: r => r.k, test: [0, 1, 2] },
    { name: 'n-glide (h+k=2n)', valueOf: r => r.h + r.k, test: [1, 1, 2] },
    { name: 'd-glide (h+k=4n)', valueOf: r => r.h + r.k, modulus: 4, test: [1, 1, 4] },
  ],
  hhl: [
    { name: 'c-glide (l=2n)', valueOf: r => r.l, test: [0, 1, 2] },
    { name: 'n-glide (h+l=2n)', valueOf: r => r.h + r.l, test: [1, 1, 2] },
    { name: 'd-glide (2h+l=4n)', valueOf: r => 2 * r.h + r.l, modulus: 4, test: [2, 1, 4] },
    { name: 'h=2n', valueOf: r => r.h, test: [1, 0, 2] },
  ],
};

function scanZonalConditions(reflections, zoneName) {
  const relevant = reflections.filter(ZONE_SELECTORS[zoneName]);
  return ZONAL_CONDITIONS[zoneName].map(cond => ({
    name: cond.name,
    ...judgeCondition(residueBuckets(relevant, cond.valueOf, cond.modulus || 2)),
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
      hhl: scanZonalConditions(reflections, 'hhl'),
    },
  };
}


// --- Shared with absences-viz.js -----------------------------------------
// True if the reflection (h,k,l) is forbidden by a condition that the scan
// judged to HOLD. The plots use this so a cell is only drawn "present" if it
// is not systematically forbidden — keeping picture and table consistent
// (a single-cell I/sigma test cannot tell a 0.1% leakage signal from a real
// reflection; the class-level verdict can).
// The period that actually defines an axis: among the periods judged "holds",
// the one with the lowest violation% (ties within 1 point go to the larger
// period, e.g. a real 4_1 where periods 2 and 4 both look clean). Using EVERY
// holding period would be wrong: a larger period such as 4 or 6 can be
// accepted spuriously by the median rule when about half of its "forbidden"
// class (e.g. k = 2 mod 4 on a 2_1 axis) is genuinely strong, which would
// wrongly grey out those allowed reflections in the plots. Returns 1 if none holds.
function effectiveAxialPeriod(results) {
  const holding = results.filter(r => r.verdict === 'holds');
  if (holding.length === 0) return 1;
  let best = holding[0];
  for (const r of holding.slice(1)) {
    if (r.violationPercent < best.violationPercent - 1 ||
        (Math.abs(r.violationPercent - best.violationPercent) <= 1 && r.period > best.period)) best = r;
  }
  return best.period;
}

function isForbiddenByHeldConditions(scan, h, k, l) {
  const r = { h, k, l };

  for (const [name, conds] of Object.entries(CENTERING_CANDIDATES)) {
    if (conds.length === 0 || scan.centering[name].verdict !== 'holds') continue;
    if (conds.some(c => mod(c.valueOf(r), c.modulus) !== 0)) return true;
  }

  for (const [axis, cls] of Object.entries(AXIAL_CLASSES)) {
    if (!cls.select(r)) continue;
    const idx = cls.indexOf(r);
    const period = effectiveAxialPeriod(scan.axial[axis]);
    if (period > 1 && mod(idx, period) !== 0) return true;
  }

  for (const [zone, conds] of Object.entries(ZONAL_CONDITIONS)) {
    if (!ZONE_SELECTORS[zone](r)) continue;
    conds.forEach((c, i) => {
      if (scan.zonal[zone] && scan.zonal[zone][i] && scan.zonal[zone][i].verdict === 'holds' && mod(c.valueOf(r), c.modulus || 2) !== 0) r._forbidden = true;
    });
    if (r._forbidden) return true;
  }
  return false;
}
