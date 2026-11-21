# XFilePrep

A browser-based tool for space-group determination from single-crystal diffraction data, with SHELX `.ins` file generation.

Runs entirely client-side — open `index.html` directly in a browser, no server or build step required.

## What it does

1. **File & Cell** — Load reflection data (SHELX `.hkl`, CIF, `.p4p`) or a `.ins` file; enter or edit unit cell parameters and the chemical formula.
2. **Laue Symmetry** — Intensity-weighted reciprocal-lattice planes (0kl, h0l, hk0) drawn with the true reciprocal metric; circle size follows the mean intensity per lattice point (averaged over exact repeats only, never over symmetry equivalents). A **layer viewer** pages through the reciprocal layers hk*L*, h*K*l and *H*kl with a slider (or ▶ Play); red rings mark gaps — positions the lattice allows that were neither measured nor covered by a symmetry equivalent. Below that, a table of **merging statistics per resolution shell** (shells of equal d*³ width): unique reflections, completeness, redundancy, ⟨I/σ⟩, I > 2σ, CC½, Rint, Rmeas and Rpim, plus curves of the main quantities — for orientation only, nothing is cut or written to the .ins/.hkl.
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
