# XFilePrep

A browser-based tool for space-group determination from single-crystal X-ray diffraction data, with SHELX `.ins` / `.hkl` generation — in the spirit of XPREP.

Runs entirely client-side: open `index.html` in a browser. No server, no build step, no external dependencies.

## Workflow

1. **File & Cell** — Load reflections (SHELX `.hkl`, CIF with a `_refln` loop or an embedded `_shelx_hkl_file`), a `.ins`/`.res` or a `.p4p`. Cell, wavelength, formula, Z and crystal size are taken from the files where available and can be edited. A crystal-data block shows the Z estimate, formula weight, density, F(000), Z′, μ, the f″ of each element and, with a crystal size, μ·r and the expected T<sub>min</sub>/T<sub>max</sub>.
2. **Laue Symmetry** — R<sub>int</sub> for each Laue class the metric allows. Intensity-weighted reciprocal-lattice planes and a layer viewer (hk*L*, h*K*l, *H*kl) that marks gaps in the data, optionally completed from symmetry equivalents. Merging statistics per resolution shell: completeness, redundancy, ⟨I/σ⟩, CC½, R<sub>int</sub>, R<sub>meas</sub>, R<sub>pim</sub>.
3. **Twins & Abs. Structure** — Its own tab, because twinning and an inverted model are easy to miss; findings are flagged on every tab. An overview table gives one line per check; the details sit in collapsible sections that open by themselves when something is wrong.
   - Lattice symmetry beyond the Laue class (Le Page search for 2-fold axes, independent of the cell choice) and the resulting (pseudo-)merohedral twin laws, tested with the H test, with an estimated twin fraction and an N(H) plot.
   - Non-merohedral twins: lattice 2-folds with 1–6° obliquity, tested on the overlapping reflections only. Also 2-folds of twin index 2 (matrix with halves, e.g. l′ = l + h/2 in a monoclinic cell with a long axis), which overlap exactly for every second reflection; these can be refined with TWIN in HKLF 4 mode, and the app gives the TWIN line. The overlap tolerance scales with the shortest reciprocal axis.
   - Inversion (racemic) twinning from the Bijvoet differences, compared with the anomalous signal expected from the composition — a model-free, sign-less analogue of the Flack parameter.
   - Absolute structure from a final CIF (refined model in `_shelx_res_file`, unmerged data in `_shelx_hkl_file`): Parsons quotients, Hooft y (with P2/P3) and a classical Flack x with the model held fixed, next to the value SHELXL wrote into the CIF. Supports LATT/SYMM, short and long SFAC, DISP, free variables, riding H, anisotropic U, EXTI, HKLF 4 matrices, OMIT/SHEL and merohedral twins (TWIN with BASF, also combined with inversion as in TWIN … -4: the model intensities include the twin law and x is the inversion share); not HKLF 5 or SWAT.
   - Merohedral and inversion twin laws can be written to the `.ins` as `TWIN`/`BASF` (opt-in).
4. **Absences** — Lattice centring, screw axes and glide planes (including *d*-glides and diagonal planes), with zone plots and axial strips. Conditions that no space group of the Laue class can have are flagged as pseudo-symmetry. A table of **pseudo-translations** (reflection classes that are weak but present) warns of pseudo-centring and sub-cells.
5. **E² Statistics** — Wilson plot, ⟨|E²−1|⟩ overall and per resolution shell, N(z) test. Needs a composition.
6. **Space Group** — Candidates ranked against the observed absences, with an explanation of why close runners-up lose, a warning when the top candidates cannot be separated, and a space-group diagram (triclinic, monoclinic, orthorhombic) in the style of International Tables. A centrosymmetric/acentric tie is resolved by E² — or by the space group a loaded CIF reports.
7. **SHELX .ins and .hkl** — A ready-to-use `.ins` for the selected candidate. If that candidate needs another setting than the measured one, the cell is transformed and the matrix written on the `HKLF` line; the `.hkl` itself stays untouched. Reflections from a CIF can be saved as a SHELX `.hkl`.
8. **Report** — Everything the other tabs found in one document, with a live preview. It lists the input files (size, date, SHA-256), the result with the source of every decision (automatic, by hand, from a file, assumed or estimated), a list of points to check (twins, ties, unresolved candidates, pseudo-translations, an E² verdict that contradicts the selected group, low completeness, strong absorption, an assumed wavelength or estimated Z …), and the details of every step with the plots, the zone plots and the space-group diagram (each can be switched off). Saved as a self-contained `.html` (prints cleanly to PDF from the browser; plots enlarge on click, collapsed parts open for printing) or as a plain-ASCII text listing (`…_report.txt`) in the spirit of XPREP's .prp file. A free-text field for sample, operator and remarks goes into the header.

## Tablet

The same page works on a tablet; there is no separate version and no switch.

- **Layout** follows the device: when the primary pointer is a finger (`pointer: coarse`, e.g. an iPad without trackpad), controls, tabs and fields are finger-sized, the tab bar stays at the top and scrolls sideways, wide tables scroll in their own box, and files are opened with a *Choose files…* button.
- **Behaviour** follows each input: a tap on a reflection, a space-group-diagram symbol or a dotted-underlined term shows its information in a bar at the bottom of the screen instead of a hover tooltip; switches with an explanation get an ⓘ. Swipe sideways over the layer viewer or the Bijvoet zones to change the layer; double-tap the 3D view to reset it. With a mouse or trackpad (also on an iPad) the hover tooltips work as on the desktop.

Safari on iPad cannot run `index.html` straight from the Files app; open it from a web address (e.g. GitHub Pages). *Add to Home Screen* in Safari gives it an icon.

## Methods

- **Symmetry** is derived at runtime from Hall symbols and their operators (all 530 settings of International Tables Vol. B), not hard-coded per space group. Reflection conditions follow from the operators.
- **Settings:** transformation matrices between settings are derived and verified from the symmetry operators. All monoclinic cell choices and axis settings are covered; the least oblique cell is preferred. Origin choices 1/2 are recognised as pure origin shifts.
- **Scoring:** a condition counts as violated relative to how many allowed reflections are significant, so weak data do not penalise the true space group. Candidates with identical predicted absences are reported as indistinguishable rather than ranked.
- **Pseudo-symmetry vs. absences:** reflections are only treated as forbidden if a space group of the current Laue class can actually forbid them.
- **Twin fraction:** α = ½ − ⟨H⟩ (Yeates 1988), from twin-related pairs selected by their summed intensity. Pseudo-symmetry raises it, noise lowers it — a start value for BASF, not a result.
- **Inversion twin:** σ is rescaled from repeated measurements first; expected signal ⟨ΔI²⟩/⟨I²⟩ = 4 Σf″² / Σf₀²(s).
- **Pseudo-translations:** mean and median intensity of each candidate weak class relative to the rest, ignoring reflections with a zero index. A strong pseudo-translation is flagged on the E² tab (it makes data look more centrosymmetric) and on tab 3 (it raises the twin fraction).
- **CC½** by the σ–τ method (Assmann, Brehm & Diederichs 2016). **R<sub>σ</sub>** = Σσ(⟨I⟩) / Σ⟨I⟩ over the merged reflections.
- **Notation:** space-group and point-group symbols on screen follow International Tables (italic letters, upright digits, subscript screw components, overbars). Files keep the plain ASCII symbol.

## Limitations

- HKLF 5 reflection files are not supported as input yet; the app cannot write them either. A non-merohedral twin (index 1) is reported, but needs an HKLF 5 file from elsewhere (e.g. PLATON, TWINABS).
- No hexagonal ↔ rhombohedral axes transformation.
- The space-group diagram covers triclinic, monoclinic and orthorhombic groups only.
- Absences cannot distinguish centrosymmetric from acentric groups with the same conditions, nor enantiomorphic pairs (*P*4₁ / *P*4₃); E² and the refinement have to decide.
- μ is computed from NIST FFAST and may differ from SHELXL by about 1 %.

## Data sources

- Hall symbols: S. R. Hall, *Acta Cryst.* A37, 517 (1981); International Tables Vol. B.
- f″: International Tables Vol. C, Table 4.2.6.8 (ten Kα lines, Ti to U), as used by SHELXL. Two misprints in the printed table are corrected (Ca at Co Kα, U at Au Kα). For other elements or wavelengths: NIST FFAST.
- μ/ρ: NIST FFAST (C. T. Chantler et al., version 2.1, 2005).
- Atomic weights: IUPAC standard values.
- Lattice symmetry: Y. Le Page, *J. Appl. Cryst.* 15, 255 (1982). H test: T. O. Yeates, *Acta Cryst.* A44, 142 (1988).

## Project structure

```
index.html
static/css/style.css          styling, light and dark mode
static/js/
  hkl-parser.js, cif-reflections.js, cell-parser.js, ins-parser.js, p4p-parser.js
                              file input
  cell-metric.js, d-spacing.js, recip-geometry.js
                              cell metric, d-spacings, reciprocal geometry
  hall-symbols-data.js, hall-symbol-parser.js, reflection-conditions.js
                              space-group symmetry and reflection conditions
  laue-merge.js, shell-stats.js
                              Laue class, merging statistics
  absences.js, pseudo-translation.js
                              systematic absences, pseudo-translations
  candidate-matcher.js, setting-transform.js
                              candidate ranking, setting transformations
  twin-check.js, anomalous-check.js
                              twin laws, H test, inversion twin
  e-statistics.js             Wilson plot, E² statistics, N(z)
  crystal-data.js, absorption-data.js, dispersion-table.js, scattering-factors.js
                              crystal data, μ, f″, form factors
  ins-writer.js, hkl-writer.js
                              SHELX output
  *-viz.js, layer-viewer.js, laue-3d.js, sg-diagram.js, reflection-marks.js,
  plot-theme.js, plot-zoom.js, view-modal.js
                              plots and diagrams
  sg-format.js                symbol notation
  report-writer.js            report tab: HTML report and text listing
  tabs.js, app.js             user interface
```
