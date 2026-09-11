// candidate-matcher.js
// Narrows the 530 Hall-symbol table entries down to space groups
// consistent with what was actually measured: the determined Laue class
// (laue-merge.js), the observed lattice centering, and the observed
// axial/zonal absences (absences.js) — scored against the CALCULATED
// conditions for each candidate (reflection-conditions.js).
//
// This is the step that finally connects every earlier module: cell
// metric narrows the Laue class candidates, merging statistics pick the
// real one, the absence scan characterizes the data, and this module
// asks "which of the 530 known settings could have produced exactly
// this pattern?".

// Extracts the rotation-only part of a space group's operator set (i.e.
// its point group, ignoring translations and lattice centering) so it can
// be compared against a Laue class's generator-closure matrix set from
// laue-merge.js.
function rotationOnlySet(operators) {
  const seen = new Map();
  for (const op of operators) {
    const key = op.R.map(r => r.join(',')).join(';');
    if (!seen.has(key)) seen.set(key, op.R);
  }
  return Array.from(seen.values());
}

// Whether a Hall symbol's own point group contains an inversion center —
// i.e. whether the CRYSTAL (not just its diffraction pattern) is truly
// centrosymmetric. Used to warn when a ".ins from top candidate" pick is
// actually an unresolved tie between an acentric and a centrosymmetric
// option — absences alone cannot break that tie (see the E² statistics
// module), so silently taking "whichever sorted first" risks generating
// a .ins for the wrong one exactly when it matters most.
function hallSymbolIsCentrosymmetric(hallSymbol) {
  const rotations = rotationOnlySet(generateSpaceGroupOperators(hallSymbol));
  return rotations.some(R => R.every((row, i) => row.every((v, j) => v === (i === j ? -1 : 0))));
}

function matrixSetKey(matrices) {
  return matrices.map(m => m.map(r => r.join(',')).join(';')).sort().join('|');
}

// Combines a rotation set with its own point-inversion images (R and -R
// for every R already present) — exactly what Friedel's law does to a
// crystal's true point group in a diffraction pattern: intensities can
// never distinguish a reflection from its inverse, so the OBSERVED
// (diffraction/Laue) symmetry is always the smallest centrosymmetric
// supergroup of the crystal's real point group, whatever that real point
// group is. A Laue class IS one such centrosymmetric supergroup by
// definition, so a candidate's OWN rotation group only needs its Friedel
// completion to equal the target Laue class — not to equal it outright.
// Skipping this step (as an earlier version of this function did, via
// plain equality) silently excluded every acentric space group (P21,
// Pca21, I41cd, ...) from ever being proposed as a candidate for ANY
// Laue class, no matter how well the data matched them.
function friedelCompletion(rotationMatrices) {
  const seen = new Map();
  for (const R of rotationMatrices) {
    for (const M of [R, R.map(row => row.map(v => -v))]) {
      const key = M.map(r => r.join(',')).join(';');
      if (!seen.has(key)) seen.set(key, M);
    }
  }
  return Array.from(seen.values());
}

// Determines whether a candidate's point group (derived from its Hall
// symbol), completed for Friedel's law, matches the given Laue class
// exactly (same set of rotation matrices — not just the same order).
function candidateMatchesLaueClass(hallSymbol, laueClassName, monoclinicAxis) {
  const operators = generateSpaceGroupOperators(hallSymbol);
  const candidateRotations = rotationOnlySet(operators);
  const completedRotations = friedelCompletion(candidateRotations);

  const targetGenerators = getLaueClassGenerators(laueClassName, monoclinicAxis);
  const targetRotations = generateGroup(targetGenerators);

  return matrixSetKey(completedRotations) === matrixSetKey(targetRotations);
}

// Scores one candidate's calculated conditions against the observed scan
// from absences.js. Returns a fraction 0..1 (weighted share of checks
// that agreed) plus the details, so the UI can show its reasoning rather
// than just a score.
//
// Checks are weighted by how many reflections actually informed them
// (forbidden+allowed count), not counted flatly one-vote-each. Without
// this, a rock-solid signal from thousands of reflections (e.g. a h0l
// zone glide check) counts exactly as much as a shaky one from a
// handful of reflections on a sparse axial row — real axial classes in
// particular are often a tiny fraction of the dataset (by definition,
// only one line through reciprocal space) and are notoriously prone to
// a few spurious "forbidden but present" reflections from multiple/Umweg
// diffraction, an experimental artifact rather than a real symmetry
// violation. Flat counting let such a shaky mismatch cancel out a
// strong, well-supported one and produce an artificial tie.
function scoreCandidate(hallSymbol, observedScan) {
  const calc = deriveAllConditions(hallSymbol);
  const checks = [];
  const weightOf = r => (r.forbiddenCount || 0) + (r.allowedCount || 0);

  // How much a single check supports the candidate, as a continuous 0..1
  // fraction rather than an all-or-nothing snap at judgeCondition's 15%
  // violationPercentThreshold. That threshold is meant for the human-
  // readable "holds"/"violated" label shown elsewhere in the UI — a
  // borderline case (say 17% violation, just past the 15% cut) is still
  // overwhelmingly clean evidence, and treating it as a full-weight
  // disagreement (as a binary snap would) lets one noisy real-world
  // reflection or two, right at the margin, cost a candidate the ENTIRE
  // weight of an otherwise strongly-supporting, thousands-of-reflections
  // check. Graduated credit means a candidate whose predictions are
  // mostly-but-not-perfectly clean is scored as mostly-but-not-perfectly
  // right, not as flatly wrong.
  const matchFraction = (expectHolds, violationPercent) =>
    expectHolds ? (100 - violationPercent) / 100 : violationPercent / 100;

  const observedCenteringForScore = Object.entries(observedScan.centering)
    .find(([name, r]) => name !== 'P (primitive)' && r.verdict === 'holds')?.[0] || 'P (primitive)';
  const centeringResult = observedScan.centering[observedCenteringForScore];
  checks.push({
    label: 'centering',
    match: calc.centering === observedCenteringForScore ? 1 : 0,
    weight: weightOf(centeringResult),
  });

  for (const axis of ['h00', '0k0', '00l']) {
    const results = observedScan.axial[axis];
    const calcPeriod = calc.axial[axis];
    // 'inconclusive' means the data itself couldn't confirm OR refute
    // anything at this specific period (too few reflections, or even the
    // allowed class was too weak to trust) — that's an absence of
    // evidence, not evidence of absence, so an inconclusive period is
    // never used as the reference for a check's credit.
    const informative = results.filter(r => r.verdict !== 'inconclusive');
    if (informative.length === 0) continue;

    let reference;
    if (calcPeriod === 1) {
      // No restriction expected (a trivial candidate on this axis) — best
      // supported by whichever TESTED period looks cleanest (lowest
      // violation%): if even that best-case period still shows real
      // violation, that's exactly what "no restriction" predicts. Using
      // the single cleanest period (rather than e.g. averaging all of
      // them) mirrors the existing "trust the cleanest evidence" logic
      // used below for restrictive candidates, applied symmetrically.
      reference = informative.reduce((a, b) => (b.violationPercent < a.violationPercent ? b : a));
    } else {
      // A specific restriction is expected — judge it on ITS OWN
      // evidence (that exact period's result), not on how some other
      // period happens to look.
      reference = informative.find(r => r.period === calcPeriod);
      if (!reference) continue; // that period wasn't informative for this axis — no evidence for this candidate's specific claim
    }

    checks.push({
      label: `${axis} period ${calcPeriod}`,
      match: matchFraction(calcPeriod !== 1, reference.violationPercent),
      weight: reference.reflectionCount, // every period partitions the SAME reflection set on this axis, so the weight is period-independent
    });
  }

  for (const zone of ['0kl', 'h0l', 'hk0']) {
    const zoneNames = { '0kl': ['b-glide (k=2n)', 'c-glide (l=2n)', 'n-glide (k+l=2n)'],
      h0l: ['a-glide (h=2n)', 'c-glide (l=2n)', 'n-glide (h+l=2n)'],
      hk0: ['a-glide (h=2n)', 'b-glide (k=2n)', 'n-glide (h+k=2n)'] }[zone];
    const c = calc.zonal[zone];
    const expectations = [
      { name: zoneNames[0], expectHolds: c.periodP === 2 },
      { name: zoneNames[1], expectHolds: c.periodQ === 2 },
      { name: zoneNames[2], expectHolds: c.periodSum === 2 },
    ];
    for (const exp of expectations) {
      const observed = observedScan.zonal[zone].find(r => r.name === exp.name);
      if (!observed || observed.verdict === 'inconclusive') continue; // no evidence either way
      checks.push({
        label: `${zone} ${exp.name}`,
        match: matchFraction(exp.expectHolds, observed.violationPercent),
        weight: weightOf(observed),
      });
    }
  }

  const totalWeight = checks.reduce((a, c) => a + c.weight, 0);
  const agreeWeight = checks.reduce((a, c) => a + c.match * c.weight, 0);
  return { score: totalWeight > 0 ? agreeWeight / totalWeight : null, checks };
}

// Two Hall-table entries are indistinguishable by systematic absences
// alone exactly when their calculated conditions are identical — this is
// the correct, physically-grounded reason to treat them as one candidate
// (not "same space-group number", which is too broad; see below).
function conditionSignature(hallSymbol) {
  return JSON.stringify(deriveAllConditions(hallSymbol));
}

// Collapses settings that produce the exact same calculated absence
// pattern (most commonly an origin-choice-1/2 pair for the same physical
// symmetry — shifting the origin never changes which reflections are
// systematically absent, so those two entries are never distinguishable
// this way and should count as one candidate, not two).
//
// Important: this does NOT collapse everything sharing a space-group
// NUMBER. The monoclinic "cell choice" alternatives (b1/b2/b3 — e.g. for
// No. 14: c-glide vs n-glide vs a-glide) all use the exact same rotation
// part (same unique-axis 2-fold matrix) and differ only in the glide's
// translation component — meaning they are directly, natively testable
// against the SAME given cell/hkl with no index transformation at all,
// and their calculated conditions genuinely differ (that's the whole
// point of "cell choice" as a concept: it's naming which of three
// physically distinct, independently possible glide directions the
// crystal actually has, in the cell as measured). Excluding two of the
// three a priori — as an earlier version of this function did, by
// keeping only the first Hall-table entry per number — silently threw
// away, for example, P21/n whenever the data actually supported an
// n-glide over the table's first-listed c-glide entry for No. 14.
function deduplicateByConditionSignature(entries) {
  const seen = new Set();
  const result = [];
  for (const entry of entries) {
    const key = entry.number + '|' + conditionSignature(entry.hall);
    if (seen.has(key)) continue;
    seen.add(key);
    result.push(entry);
  }
  return result;
}

// Top-level: filters the Hall table by Laue class and centering, collapses
// settings that are indistinguishable by absences alone (see
// deduplicateByConditionSignature), scores the remainder against the
// observed absences, returns them ranked best first. Pass
// skipConditionDedup: true to see every raw matching entry, including
// origin-choice duplicates that always score identically — mainly useful
// for debugging the table itself.
function matchSpaceGroupCandidates(hallTable, laueClassName, observedScan, options = {}) {
  const monoclinicAxis = options.monoclinicAxis || 'b';

  const laueMatches = hallTable.filter(entry =>
    candidateMatchesLaueClass(entry.hall, laueClassName, monoclinicAxis));

  const observedCentering = Object.entries(observedScan.centering)
    .find(([name, r]) => name !== 'P (primitive)' && r.verdict === 'holds')?.[0] || 'P (primitive)';
  const centeringMatches = laueMatches.filter(entry => {
    const { latticeTranslations } = parseHallSymbol(entry.hall);
    return latticeTranslationsToName(latticeTranslations) === observedCentering;
  });

  const candidateEntries = options.skipConditionDedup
    ? centeringMatches
    : deduplicateByConditionSignature(centeringMatches);

  const scored = candidateEntries.map(entry => ({
    entry,
    ...scoreCandidate(entry.hall, observedScan),
  }));

  scored.sort((a, b) => (b.score ?? -1) - (a.score ?? -1));
  return scored;
}
