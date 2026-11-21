# XFilePrep

A browser-based tool for space-group determination from single-crystal diffraction data, with SHELX `.ins` file generation.

Runs entirely client-side — open `index.html` directly in a browser, no server or build step required.

## What it does

1. **File & Cell** — Load reflection data (SHELX `.hkl`, CIF, `.p4p`) or a `.ins` file; enter or edit unit cell parameters and the chemical formula.
2. **Laue Symmetry** — Intensity-weighted reciprocal-lattice planes (0kl, h0l, hk0) drawn with the true reciprocal metric; circle size follows the mean intensity per lattice point (averaged over exact repeats only, never over symmetry equivalents). A **layer viewer** pages through the reciprocal layers hk*L*, h*K*l and *H*kl with a slider (or ▶ Play); red rings mark gaps — positions the lattice allows that were neither measured nor covered by a symmetry equivalent. An optional checkbox **fill in symmetry equivalents** (off by default) completes the zone plots and the layer viewer for merged data: positions that were not measured themselves get the mean intensity of their measured equivalents under the Laue class and are drawn in orange; measured positions, the intensity scale and the gaps are unchanged. It shows the assumed symmetry, so it is meant for orientation, not for judging the Laue class. Below that, a table of **merging statistics per resolution shell** (shells of equal d*³ width): unique reflections, completeness, redundancy, ⟨I/σ⟩, I > 2σ, CC½, Rint, Rmeas and Rpim, plus curves of the main quantities — for orientation only, nothing is cut or written to the .ins/.hkl.
3. **Absences** — Scans systematic absences (centering, axial, glide/zone conditions) and visualizes them as zone plots (true reciprocal metric) and axial strips. Click a zone plot or axial strip (here or in the Laue tab), the Wilson plot or the N(z) plot, to enlarge it over a blurred page; Esc or a click outside closes it.
4. **E² Statistics** — Wilson plot, ⟨|E²−1|⟩ against the centrosymmetric/acentric references and an N(z) test (centric, acentric and twinned-acentric curves). The same ⟨|E²−1|⟩ is also shown per resolution shell, using the same shell boundaries as the merging table on the Laue tab (with a warning if the shells disagree, e.g. because of a heavy atom or pseudosymmetry). Centering-absent reflections of the observed lattice are excluded first. Needs a composition (formula + Z, or a `.ins`).
5. **Space Group** — Matches candidate space groups against the observed absences; a centrosymmetric/acentric tie is resolved with the E² verdict from the previous tab. Below the candidate table, close runners-up that lose on a clear difference are explained, and a warning appears when the top candidate is not clearly separated from the next one.
6. **SHELX .ins and .hkl** — Generates a ready-to-use `.ins` file for the selected space group, kept in sync with the current cell, formula, and candidate selection. If the selected candidate isn't the one the as-measured cell/hkl natively fit, the cell is transformed into that setting and the corresponding matrix is written on the `HKLF` line (the `.hkl` file itself is never touched). If the reflections came from a CIF (a `_refln` loop or an embedded `_shelx_hkl_file` block) rather than from a raw SHELX `.hkl`, this tab also offers to save them as a SHELX `.hkl` — untransformed, same base name as the `.ins`.

## Project structure

```
index.html
static/
  css/style.css       — styling, incl. automatic dark mode (prefers-color-scheme)
  js/
    hkl-parser.js, cif-reflections.js, cell-parser.js, cell-metric.js
                        — parsing reflection/cell data (HKL, CIF, .ins, .p4p)
    laue-merge.js       — point-group generation, Laue class determination
    absences.js, absences-viz.js
                        — systematic absence scan + plots
    recip-geometry.js, laue-viz.js
                        — reciprocal-lattice plane geometry shared by the plots; Laue intensity overview
    hall-symbols-data.js, hall-symbol-parser.js, reflection-conditions.js
                        — space-group symmetry data and derived reflection conditions
    candidate-matcher.js
                        — matches candidate space groups against observed data
    setting-transform.js
                        — derives+verifies a cell/HKL transformation matrix between two settings of the same space-group number, from their actual symmetry operators (not a trusted name lookup)
    ins-parser.js, p4p-parser.js, scattering-factors.js, d-spacing.js
    e-statistics.js     — Wilson plot / E² statistics / N(z) test / lattice filter
    shell-stats.js      — merging statistics per resolution shell (completeness, ⟨I/σ⟩, R values, CC½) and shell curves
    layer-viewer.js     — reciprocal-layer slider with gap markers
    ins-writer.js       — SHELX .ins generation
    hkl-writer.js       — SHELX .hkl export for reflections that came from a CIF
    plot-theme.js       — light/dark palette for canvas plots
    plot-zoom.js        — click-to-enlarge overlay for zone plots, axial strips and the Wilson plot
    tabs.js, app.js      — UI wiring
tests/                  — regression tests, run with node from the project root (see Notes)
```

## Tests

Run from the project root with Node (no installation needed):

| Command | What it checks | Time |
|---|---|---|
| `node tests/blind-run.js [folder]` | **Blind run** over real datasets (default `tests/data`): every CIF (reflections embedded, or in `name.hkl` / `name.hkl.txt` / `name.fcf` next to it) is solved from cell, reflections and formula only — Laue class from Rint, absences, candidates, E² tie-break — and only then compared with the space group the CIF reports. Result per file: `OK`, `TIE` (identical absences, e.g. I 4 2 2 / I -4 2 m — no absence analysis can decide), `RANK n`, `LAUE`, `MISS`, `ERROR`; exit code 1 unless all are OK/TIE. | seconds |
| `node tests/short-axis-screw-regression.js` | Screw axes along short axes, stricter-period check, violations with few allowed reflections (synthetic data). | < 1 s |
| `node tests/hall-operator-count.js` | All 530 Hall table entries: operator count, Laue-group order per IT number, origin choice 1 ↔ 2 = pure origin shift. | < 1 s |
| `node tests/shelx-symm-expansion.js` | All 530 entries: the emitted SYMM + LATT cards, expanded as SHELXL does, give exactly the Hall group, with the minimal number of cards. | < 1 s |
| `node tests/absence-signature-check.js` | All 530 entries: the true systematic absences (from the operators, up to symmetry equivalence) vs. the condition signature the app compares candidates by — settings with different absences must get different signatures, identical ones the same. Spot-checked against ITA A Table 1.6.4.8 (oC, mmm): all nine rows reproduced. | 2 s |
| `node tests/sginfo-compare.js /path/to/sginfo` | Optional: all 530 entries against **SgInfo** (Grosse-Kunstleve), an independent space-group library — IT number, operator set, LATT, SYMM count, and table completeness. Build SgInfo once with `cc -O -w -o sginfo sginfo.c sgclib.c sgio.c sgfind.c sghkl.c sgsi.c -lm` in its source folder. | seconds |
| `node tests/monoclinic-setting-regression.js` | Every setting transformation (all monoclinic pairs, orthorhombic permutations) against the systematic absences. | ~1 min |

To grow the blind test set, drop further CIF(+HKL) pairs into `tests/data` — any structure whose CIF names the space group works, ideally with raw (unmerged) reflections. Current set: COD 1543504 (P 1 2₁/a 1, merged, short b axis), COD 1501632 (I -4 2 m, merged, identical extinction symbol), COD 1561021 (P c a b = non-standard Pbca, 112 000 raw reflections), COD 2023426 (F m -3 m, merged, JANA2006 CIF), COD 2010746 (F m -3 m, only 116 merged reflections), jg6 (I 4₁ c d, 74 000 raw reflections).

## Notes

- CIF space-group metadata (Hall symbol, or H-M name + IT number) are used as an independent source for the Laue class; if they contradict the cell metric, the app falls back to the data analysis. Rint stays visible as a separate control table.
- Regression test for the CIF Laue class: `node tests/cif-laue-regression.js` (run from the project root; the other tests are named below).
- No build tools, bundlers, or external runtime dependencies — vanilla JavaScript throughout.
- Symmetry and systematic-absence logic is derived from Hall symbols and symmetry operators at runtime, not hardcoded per space group.
- Some ambiguities (e.g. centrosymmetric vs. acentric candidates with identical absences) are inherent to X-ray diffraction (Friedel's law) and cannot be resolved from absences alone — E² statistics is provided to help with those cases.
- Setting transformation covers: all monoclinic settings (cell choices 1/2/3, A/C/I centering, a-, b- and c-unique, incl. the -b1.. entries of Nos. 9/15) and orthorhombic-style axis permutation. Matrices are in SHELX/HKLF convention (a' = M a, h' = M h); the symmetry operators are checked with the coordinate matrix (M⁻¹)ᵀ — up to v2 they were conjugated with M itself, which is only correct for pure axis permutations, so every monoclinic shear (e.g. P2₁/n → P2₁/c, I2/a → C2/c) got a wrong HKLF matrix. For monoclinic groups all matrices with entries −1/0/1 are searched and the one giving the least oblique cell (shortest axes, angle closest to 90°, obtuse) is used; if the selected setting is still clearly more oblique than the measured one, the Space Group and .ins tabs say so (IUCr: least oblique cell choice is recommended). Regression test: `node tests/monoclinic-setting-regression.js` (checks every matrix independently against the systematic absences, ~1 min). Not covered: origin choice (e.g. Fd-3m Origin 1/2) and the hexagonal↔rhombohedral axes transformation for R-centered trigonal groups (a fractional, not integer, transform) — these are reported as "not found" rather than guessed at.
- Match scores are always computed against the raw, as-measured HKL indices — never against a transformed re-indexing. So a candidate's standard-setting sibling (e.g. Pbca vs. the top-scoring Pcab) can legitimately score far lower even though it describes the identical physical structure: the score reflects how well the *unmodified* indices happen to fit that candidate's own axis convention, not the structure's correctness. This is why the Space Group tab surfaces a top candidate's standard sibling even when its score is low — the low score is expected there, not a red flag.
- Completeness and the shell boundaries only count reflections as possible that the lattice centering and the systematic absences found by the absence scan (screw axes, glide planes) allow, so data already cleaned of absences (deposited CIF/.hkl files) are not rated incomplete. For CIF data, the Laue class used by all statistics (E², shells, layer viewer) comes from the CIF space group.
- Pseudo-centering (e.g. weak but present h+l-odd reflections of a primitive cell) is no longer read as a lattice absence: the "weak tail" rule of the absence scan only applies if the forbidden class is also weak relative to the allowed one (mean intensity < 10 %). Regression test: `node tests/pseudo-centering-regression.js`.
- `summarizeHKL` uses loops instead of `Math.min(...array)`: large files (more than about 65 000 reflections on iPad/iPhone, about 120 000 elsewhere) used to overflow the call stack.
- Candidate scoring: a candidate that predicts only the weaker axial condition (00l: l = 2n) is penalised if the data obey the stronger one (l = 4n, a 4₁ axis). Before, I 4 c m and I 4₁ c d scored the same. The Match column shows one decimal. Regression test: `node tests/stricter-period-regression.js`.
- The absence scan and the candidate scoring know the **d-glide** (quarter translation: a sum modulo 4, e.g. 0kl: k+l = 4n as in F d d 2) and the **diagonal plane hhl** (glide/mirror planes ⊥ [1-10] in tetragonal, ⊥ [210] in hexagonal/trigonal and {110} in cubic groups; e.g. 2h+l = 4n in I 4₁ c d and I a -3 d). The theoretical conditions are derived from the operators as the set of allowed index pairs (mod 12); a candidate "expects" a condition to hold if all its allowed pairs satisfy it, so a- and b-glide together, or a d-glide, also imply the n-condition. Neutral checks (a d/hhl condition that is violated while the candidate predicts nothing there) do not count, so the absolute match values do not drop. Regression test: `node tests/diagonal-conditions-regression.js`. Not covered: rhombohedral axes, the second diagonal plane h-hl (equivalent under 3-/4-/6-fold symmetry), a dedicated hhl plot.
- Space Group tab: below the candidate table, a block explains why close runners-up lose on a **clear** difference (a condition the scan found clearly held or violated, where the candidate predicts the opposite). If a candidate is within 1.5 points of the top and **no** such difference separates them, a warning ("not clearly separated") appears; weak data, pseudosymmetry, a condition the scan does not cover and enantiomorphic pairs (P 4₁/P 4₃) are listed as possible causes.
- Shell table (Laue tab): the CC½ column uses the σ–τ method (Assmann, Brehm & Diederichs 2016, as in XDSCC12 and CrystFEL's partialator; see the XDSwiki page "CC1/2"): deterministic, no random split, independent of the σ values; orange below 30 %. Two curves below the table show completeness, I > 2σ, CC½ and Rmeas (percent) against (sin θ/λ)², and ⟨I/σ⟩ on a log axis with a line at 2. For merged data (redundancy 1), CC½ and the R values are "—". Regression test against the worked example of the wiki page (0.9458): `node tests/cc12-regression.js`.
- Screw axes along a **short axis**: an axial row needs only 3 allowed reflections (zones: 5). With b = 5.5 Å at 0.77 Å resolution 0k0 holds just 020/040/060, so the 2₁ was always "undetermined" and P 2/a, P a and P 2₁/a tied (COD 1543504, really P 2₁/a). The "why runners-up score lower" block also counts such an axial check as decisive, and says when the forbidden reflections are missing from the file altogether (deposited data with absences removed). Regression test: `node tests/short-axis-screw-regression.js` (incl. a P 2/a control that must not gain a 2₁). Together with this, the "stricter period" check (I 4 c m vs. I 4₁ c d, see above) is now asked only among the reflections the predicted period allows (e.g. "of the k = 2n, are the k = 4n+2 absent?"); over the whole row the candidate's own absences diluted it, and P 2₁/n lost to P n / P 2/n.
- **Few allowed reflections, strong forbidden ones:** a condition with fewer allowed reflections than the minimum can't be confirmed, but it is now judged *violated* when at least two forbidden reflections are significant and reach ≥ 10 % of the allowed mean intensity (COD 1501632, I -4 2 m: 00l holds only 004/008 for l = 4n while 002/006 are strong, so I 4₁ 2 2 drops out instead of tying).
- **Identical extinction symbols:** candidates whose predicted absences are exactly the same (e.g. I 4 2 2, I 4 m m, I -4 m 2, I -4 2 m, I 4/m m m for "I - - -") are reported as indistinguishable by absences instead of "not clearly separated". If a loaded CIF reports one of the top-scoring candidates, that one is preselected (and kept by the E² step if it fits its verdict); otherwise the order stays as before. Regression cases in `tests/short-axis-screw-regression.js`.
- **F-centred data:** an F lattice satisfies the A, B and C conditions as well; the observed centring is now the most specific one that holds (F > I > R > A/B/C). Before, F data were taken as A-centred and e.g. m-3m had no candidates at all (COD 2023426).
- **Hall parser, N = 1:** "-1n", "-1ab", "-1bw" at the third position of origin-choice-1 symbols (Nos. 59, 125, 129, 130, 133, 137, 138, 141, 142) were read as -3 along [111]; they are the inversion with a translation. These settings had cubic rotation groups and never matched their Laue class.
- **Origin choices:** origin choice 1 vs 2 is found as a pure origin shift (shift search in 24ths — F d -3 m, F d d d, I 41/a m d differ by 1/8); shown as "matches as measured (other origin choice)", no HKLF matrix.
- **SYMM cards:** one per rotation, no bare centring vectors (I 41 c d: 7 instead of 15; F d -3 m:2: 23 instead of 95). The groups were correct before, only redundant.
- **Condition signature:** it was `JSON.stringify` of the derived conditions, in which the zone sets (`Set`) became `{}` — all zonal conditions of the full model, hhl in particular, were ignored. 50 pairs with different absences (P 4₂ 2 2 / P 4₂ m c, F 4 3 2 / F -4 3 c, P 6₃ 2 2 / P 6₃ m c, ...) looked identical, 6 identical ones (F d d d:1/:2, R 3 2 / R 3 m, ...) looked different. Affected the "identical absences" note, CIF preselection and setting deduplication — not the scores themselves. The signature is now exactly what the scoring tests.
