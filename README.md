# XFilePrep

A browser-based tool for space-group determination from single-crystal diffraction data, with SHELX `.ins` file generation.

Runs entirely client-side — open `index.html` directly in a browser, no server or build step required.

## What it does

1. **File & Cell** — Load reflection data (SHELX `.hkl`, CIF, `.p4p`) or a `.ins` file; enter or edit unit cell parameters and the chemical formula.
2. **Absences** — Scans systematic absences (centering, axial, glide/zone conditions) and visualizes them as zone and axial-strip plots.
3. **Space Group** — Determines the Laue class from the reflection data, matches candidate space groups against the observed absences, and (optionally) runs E² statistics to help resolve centrosymmetric/acentric ties that absences alone cannot distinguish. Also shows, per candidate, whether its cell/hkl match what was actually measured or would need a setting transformation (see `setting-transform.js` below).
4. **SHELX .ins** — Generates a ready-to-use `.ins` file for the selected space group, kept in sync with the current cell, formula, and candidate selection. If the selected candidate isn't the one the as-measured cell/hkl natively fit, the cell is transformed into that setting and the corresponding matrix is written on the `HKLF` line (the `.hkl` file itself is never touched).

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
    hall-symbols-data.js, hall-symbol-parser.js, reflection-conditions.js
                        — space-group symmetry data and derived reflection conditions
    candidate-matcher.js
                        — matches candidate space groups against observed data
    setting-transform.js
                        — derives+verifies a cell/HKL transformation matrix between two settings of the same space-group number, from their actual symmetry operators (not a trusted name lookup)
    ins-parser.js, p4p-parser.js, scattering-factors.js, d-spacing.js
    e-statistics.js     — Wilson plot / E² statistics
    ins-writer.js       — SHELX .ins generation
    plot-theme.js       — light/dark palette for canvas plots
    tabs.js, app.js      — UI wiring
```

## Notes

- No build tools, bundlers, or external runtime dependencies — vanilla JavaScript throughout.
- Symmetry and systematic-absence logic is derived from Hall symbols and symmetry operators at runtime, not hardcoded per space group.
- Some ambiguities (e.g. centrosymmetric vs. acentric candidates with identical absences) are inherent to X-ray diffraction (Friedel's law) and cannot be resolved from absences alone — E² statistics is provided to help with those cases.
- Setting transformation currently covers: monoclinic cell choice (b-unique only — a/c-unique isn't implemented) and orthorhombic-style axis permutation. It does NOT cover origin choice (e.g. Fd-3m Origin 1/2) or the hexagonal↔rhombohedral axes transformation for R-centered trigonal groups (a fractional, not integer, transform) — these are reported as "not found" rather than guessed at.
- Match scores are always computed against the raw, as-measured HKL indices — never against a transformed re-indexing. So a candidate's standard-setting sibling (e.g. Pbca vs. the top-scoring Pcab) can legitimately score far lower even though it describes the identical physical structure: the score reflects how well the *unmodified* indices happen to fit that candidate's own axis convention, not the structure's correctness. This is why the Space Group tab surfaces a top candidate's standard sibling even when its score is low — the low score is expected there, not a red flag.
