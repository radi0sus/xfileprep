// pseudo-translation.js
// Pseudo-centring / pseudo-translation: most of the scattering matter repeats
// with a translation t that is not a lattice vector of the true cell, so the
// reflections with h·t non-integral are systematically WEAK, but present.
// Example: pf12c_26_o (P -1): h+k+l odd about 3x weaker than h+k+l even,
// a pseudo-I cell.
//
// Why it matters:
//   - with weak data the weak class looks absent: wrong (centred or smaller)
//     cell, wrong Z, an averaged structure;
//   - two intensity populations widen the E² distribution: ⟨|E²−1|⟩ goes up
//     and an acentric structure can look centrosymmetric;
//   - pseudo-symmetry makes twin-related intensities similar: the H test
//     overestimates α;
//   - solution and refinement: the substructure comes first, strong
//     correlations, odd ADPs; the weak class must be measured well.
//
// Test: for every candidate translation, the mean intensity of the class it
// makes weak over that of the rest. Reflections with a zero index are left
// out, so glide planes and screw axes cannot fake a low ratio; reflections
// that the true lattice forbids are left out as well. Several translations
// are found one after the other: once the strongest one is found, the
// others are tested within its strong class only, so the classes it already
// explains (A, B, C when I is pseudo) drop out.

const _ptMod = (x, n) => ((x % n) + n) % n;

const PSEUDO_TRANSLATION_CLASSES = [
  { key: 'I', vector: '(½, ½, ½)', label: 'pseudo-I', cond: 'h+k+l odd', weak: r => _ptMod(r.h + r.k + r.l, 2) !== 0 },
  { key: 'A', vector: '(0, ½, ½)', label: 'pseudo-A', cond: 'k+l odd', weak: r => _ptMod(r.k + r.l, 2) !== 0 },
  { key: 'B', vector: '(½, 0, ½)', label: 'pseudo-B', cond: 'h+l odd', weak: r => _ptMod(r.h + r.l, 2) !== 0 },
  { key: 'C', vector: '(½, ½, 0)', label: 'pseudo-C', cond: 'h+k odd', weak: r => _ptMod(r.h + r.k, 2) !== 0 },
  { key: 'Robv', vector: '(⅔, ⅓, ⅓)', label: 'pseudo-R (obverse)', cond: '−h+k+l ≠ 3n', weak: r => _ptMod(-r.h + r.k + r.l, 3) !== 0 },
  { key: 'Rrev', vector: '(⅓, ⅔, ⅓)', label: 'pseudo-R (reverse)', cond: 'h−k+l ≠ 3n', weak: r => _ptMod(r.h - r.k + r.l, 3) !== 0 },
  { key: 'a2', vector: '(½, 0, 0)', label: 'a/2 sub-cell', cond: 'h odd', weak: r => _ptMod(r.h, 2) !== 0 },
  { key: 'b2', vector: '(0, ½, 0)', label: 'b/2 sub-cell', cond: 'k odd', weak: r => _ptMod(r.k, 2) !== 0 },
  { key: 'c2', vector: '(0, 0, ½)', label: 'c/2 sub-cell', cond: 'l odd', weak: r => _ptMod(r.l, 2) !== 0 },
  { key: 'a3', vector: '(⅓, 0, 0)', label: 'a/3 sub-cell', cond: 'h ≠ 3n', weak: r => _ptMod(r.h, 3) !== 0 },
  { key: 'b3', vector: '(0, ⅓, 0)', label: 'b/3 sub-cell', cond: 'k ≠ 3n', weak: r => _ptMod(r.k, 3) !== 0 },
  { key: 'c3', vector: '(0, 0, ⅓)', label: 'c/3 sub-cell', cond: 'l ≠ 3n', weak: r => _ptMod(r.l, 3) !== 0 },
];

const PT_MIN_COUNT = 30;        // reflections needed in each class
const PT_STRONG_RATIO = 0.5;    // ⟨I⟩weak/⟨I⟩rest at or below: pseudo-translation
const PT_HINT_RATIO = 0.75;     // up to here: weak hint
const PT_ABSENT_RATIO = 0.03;   // below, with hardly any strong ones: a real absence

function _ptMedian(a) {
  if (!a.length) return NaN;
  const v = Float64Array.from(a).sort();
  const m = v.length >> 1;
  return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
}

function _ptClassStats(refl, weakTest) {
  let nW = 0, nS = 0, sumW = 0, sumS = 0, sigW = 0, sigS = 0, hash = 0;
  const iW = [], iS = [];
  for (let i = 0; i < refl.length; i++) {
    const r = refl[i];
    const strong = r.sigma > 0 && r.f2 / r.sigma >= 3;
    if (weakTest(r)) { nW++; sumW += r.f2; iW.push(r.f2); if (strong) sigW++; hash = (hash * 31 + i) % 1000000007; }
    else { nS++; sumS += r.f2; iS.push(r.f2); if (strong) sigS++; }
  }
  if (nW < PT_MIN_COUNT || nS < PT_MIN_COUNT) return null;
  const meanW = sumW / nW, meanS = sumS / nS;
  if (!(meanS > 0)) return null;
  // median ratio: robust against a few very strong low-angle reflections
  const medS = _ptMedian(iS);
  const medianRatio = medS > 0 ? Math.max(0, _ptMedian(iW)) / medS : NaN;
  return { nWeak: nW, nRest: nS, ratio: meanW / meanS, medianRatio, strongWeak: sigW / nW, strongRest: sigS / nS, signature: `${nW}:${hash}` };
}

// the larger of mean and median ratio decides (conservative: a few very
// strong reflections lower the mean, many weak ones the median)
function _ptRatio(st) {
  return Number.isFinite(st.medianRatio) ? Math.max(st.ratio, st.medianRatio) : st.ratio;
}

function _ptVerdict(st) {
  if (!st) return 'too few';
  const q = _ptRatio(st);
  if (q < PT_ABSENT_RATIO && st.strongWeak < 0.05) return 'absent';
  if (q <= PT_STRONG_RATIO) return 'pseudo-translation';
  if (q <= PT_HINT_RATIO) return 'weak hint';
  return 'none';
}

// reflections   raw (or merged) reflections {h, k, l, f2, sigma}
// latticeTest   r => true if the TRUE lattice allows r (null: primitive)
// Returns { table: [...all classes on the full set], found: [...chain],
//           weakTest: r => true for reflections off the pseudo-lattice | null,
//           nUsed }
function detectPseudoTranslations(reflections, latticeTest = null) {
  const base = reflections.filter(r => r.h !== 0 && r.k !== 0 && r.l !== 0
    && Number.isFinite(r.f2) && (!latticeTest || latticeTest(r)));
  // in a centred lattice some classes coincide (I: h+k odd = l odd): the
  // later one is folded into the first as an alias
  const table = [];
  for (const c of PSEUDO_TRANSLATION_CLASSES) {
    const st = _ptClassStats(base, c.weak);
    const twin = st && table.find(t => t.stats && t.stats.signature === st.signature);
    if (twin) { twin.aliases.push(c); continue; }
    table.push({ ...c, stats: st, verdict: _ptVerdict(st), aliases: [] });
  }

  // greedy chain: strongest translation first, the rest tested inside its strong class
  const found = [];
  let pool = base;
  const used = new Set();
  for (let step = 0; step < 3; step++) {
    let best = null;
    for (const c of PSEUDO_TRANSLATION_CLASSES) {
      if (used.has(c.key)) continue;
      const st = _ptClassStats(pool, c.weak);
      const v = _ptVerdict(st);
      if (v !== 'pseudo-translation' && v !== 'weak hint') continue;
      if (!best || _ptRatio(st) < _ptRatio(best.stats)) best = { ...c, stats: st, verdict: v };
    }
    if (!best) break;
    const row = table.find(t => t.key === best.key || t.aliases.some(a => a.key === best.key));
    best.aliases = row ? [row, ...row.aliases].filter(a => a.key !== best.key) : [];
    found.push(best);
    used.add(best.key);
    best.aliases.forEach(a => used.add(a.key));
    pool = pool.filter(r => !best.weak(r));
  }
  const strongFound = found.filter(f => f.verdict === 'pseudo-translation');
  const weakTest = strongFound.length ? (r => strongFound.some(f => f.weak(r))) : null;
  return { table, found, strongFound, weakTest, nUsed: base.length };
}

function pseudoTranslationRatio(st) { return _ptRatio(st); }

if (typeof module !== 'undefined') module.exports = { detectPseudoTranslations, pseudoTranslationRatio, PSEUDO_TRANSLATION_CLASSES };
