# XFilePrep

A browser-based tool for space-group determination from single-crystal diffraction data, with SHELX `.ins` file generation.

Runs entirely client-side — open `index.html` directly in a browser, no server or build step required.

## What it does

1. **File & Cell** — Load reflection data (SHELX `.hkl`, CIF, `.p4p`) or a `.ins` file; enter or edit unit cell parameters and the chemical formula.
2. **Laue Symmetry** — Intensity-weighted reciprocal-lattice planes (0kl, h0l, hk0) drawn with the true reciprocal metric; circle size follows the mean intensity per lattice point (averaged over exact repeats only, never over symmetry equivalents). Below the plots, a table of merging statistics per resolution shell (equal d*³ width): unique reflections, completeness, redundancy, ⟨I/σ⟩, I > 2σ, Rint, Rmeas and Rpim — for orientation only, nothing is cut or written to the .ins/.hkl. Außerdem ein Layer-Viewer: per Schieberegler (oder ▶ Play) durch die reziproken Ebenen hk*L*, h*K*l, *H*kl blättern, mit roten Ringen für Lücken (erlaubt, aber weder gemessen noch über Symmetrie abgedeckt).
3. **Absences** — Scans systematic absences (centering, axial, glide/zone conditions) and visualizes them as zone plots (true reciprocal metric) and axial strips. Click a zone plot or axial strip (here or in the Laue tab), the Wilson plot or the N(z) plot, to enlarge it over a blurred page; Esc or a click outside closes it.
4. **E² Statistics** — Wilson plot, ⟨|E²−1|⟩ against the centrosymmetric/acentric references and an N(z) test (centric, acentric and twinned-acentric curves). The same ⟨|E²−1|⟩ is also shown per resolution shell, using the same shell boundaries as the merging table on the Laue tab (with a warning if the shells disagree, e.g. because of a heavy atom or pseudosymmetry). Centering-absent reflections of the observed lattice are excluded first. Needs a composition (formula + Z, or a `.ins`).
5. **Space Group** — Matches candidate space groups against the observed absences; a centrosymmetric/acentric tie is resolved with the E² verdict from the previous tab.
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
    ins-writer.js       — SHELX .ins generation
    hkl-writer.js       — SHELX .hkl export for reflections that came from a CIF
    plot-theme.js       — light/dark palette for canvas plots
    plot-zoom.js        — click-to-enlarge overlay for zone plots, axial strips and the Wilson plot
    tabs.js, app.js      — UI wiring
```

## Notes

- CIF-Raumgruppenmetadaten (Hall-Symbol bzw. H-M + IT-Nummer) werden als unabhängige Quelle für die Laue-Klasse verwendet; widersprechen sie der Zellmetrik, fällt die App auf die Datenanalyse zurück. Rint bleibt als separate Kontrolltabelle sichtbar.
- Reproduktion des Regressionstests: `node tests/cif-laue-regression.js`.
- No build tools, bundlers, or external runtime dependencies — vanilla JavaScript throughout.
- Symmetry and systematic-absence logic is derived from Hall symbols and symmetry operators at runtime, not hardcoded per space group.
- Some ambiguities (e.g. centrosymmetric vs. acentric candidates with identical absences) are inherent to X-ray diffraction (Friedel's law) and cannot be resolved from absences alone — E² statistics is provided to help with those cases.
- Setting transformation currently covers: monoclinic cell choice (b-unique only — a/c-unique isn't implemented) and orthorhombic-style axis permutation. It does NOT cover origin choice (e.g. Fd-3m Origin 1/2) or the hexagonal↔rhombohedral axes transformation for R-centered trigonal groups (a fractional, not integer, transform) — these are reported as "not found" rather than guessed at.
- Match scores are always computed against the raw, as-measured HKL indices — never against a transformed re-indexing. So a candidate's standard-setting sibling (e.g. Pbca vs. the top-scoring Pcab) can legitimately score far lower even though it describes the identical physical structure: the score reflects how well the *unmodified* indices happen to fit that candidate's own axis convention, not the structure's correctness. This is why the Space Group tab surfaces a top candidate's standard sibling even when its score is low — the low score is expected there, not a red flag.
- Die Vollständigkeit (Completeness) und die Schalengrenzen zählen nur Reflexe als möglich, die Gitterzentrierung und die im Absence-Scan als systematisch verboten erkannten Auslöschungen (Schraubenachsen, Gleitspiegelebenen) zulassen. So werden bereits von Auslöschungen bereinigte Daten (CIF/.hkl aus Datenbanken) nicht als unvollständig bewertet. Bei CIF-Daten wird die Laue-Klasse in allen Statistiken (E², Schalen, Layer-Viewer) aus der CIF-Raumgruppe übernommen.
- Pseudo-Zentrierung (z. B. schwache, aber vorhandene h+l-ungerade-Reflexe einer primitiven Zelle) wird nicht mehr als Gitterauslöschung gewertet: Die Regel „schwacher Rest" im Absence-Scan gilt nur, wenn die verbotene Klasse auch relativ zur erlaubten schwach ist (mittlere Intensität < 10 %). Regressionstest: `node tests/pseudo-centering-regression.js`.
- `summarizeHKL` nutzt Schleifen statt `Math.min(...array)`: große Dateien (> ca. 65 000 Reflexe auf iPad/iPhone, > ca. 120 000 sonst) liefen vorher in einen Stack-Überlauf.
- Kandidaten-Bewertung: Ein Kandidat, der nur die schwächere Achsenbedingung vorhersagt (00l: l = 2n), wird jetzt abgewertet, wenn die Daten die stärkere erfüllen (l = 4n, 4₁-Achse). Vorher bekamen I 4 c m und I 4₁ c d dieselbe Punktzahl. Unter der Kandidatentabelle erklärt ein Block, woran dicht folgende Kandidaten verlieren; die Match-Spalte zeigt eine Nachkommastelle. Regressionstest: `node tests/stricter-period-regression.js`. Noch nicht gescannt: Diagonalbedingungen (hhl, Gleitspiegelebenen auf den Diagonalen, d-Gleitung).
- Absence-Scan und Kandidaten-Bewertung kennen jetzt die **d-Gleitung** (Viertel-Translation: Summe mod 4, z. B. 0kl: k+l = 4n, wie in F d d 2) und die **Diagonalebene hhl** (Gleit-/Spiegelebenen ⊥ [1-10] tetragonal, ⊥ [210] hexagonal/trigonal, {110} kubisch; z. B. 2h+l = 4n bei I 4₁ c d und I a -3 d). Die theoretischen Bedingungen werden aus den Operatoren als Menge erlaubter Indexpaare (mod 12) abgeleitet; ein Kandidat „erwartet" eine Bedingung, wenn alle erlaubten Paare sie erfüllen. Dadurch gilt auch: a- und b-Gleitung zusammen bzw. eine d-Gleitung implizieren die n-Bedingung. Neutrale Checks (d-/hhl-Bedingung verletzt, Kandidat sagt nichts voraus) zählen nicht mit, damit die absoluten Match-Werte nicht sinken. Regressionstest: `node tests/diagonal-conditions-regression.js`. Nicht abgedeckt: rhomboedrische Achsen, die zweite Diagonalebene h-hl (bei 3-/4-/6-zähliger Symmetrie äquivalent), ein eigener Plot für hhl.
