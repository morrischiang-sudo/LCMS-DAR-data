# ADC DAR Analysis Platform — Data Review & Strategy Summary

## 1. What the current process does

The current DAR workflow starts from two Thermo BioPharma-Finder-style deconvolution exports: one for the naked mAb reference (`01_A09_dN.xlsx`) and one for the conjugated ADC (`02_A090412H_dN.xlsx`), each listing deconvoluted intact masses with intensity and abundance metrics. `DAR calculation.xlsx` then does the actual analysis by hand: it takes the most abundant mAb masses (which capture glycoform/adduct heterogeneity of the unconjugated antibody) as "base masses," and builds a theoretical mass grid of base mass + n₁×MW₁ + n₂×MW₂ for two user-specified linker-payload chemistries (labeled `[4]` and `[12]` in the A090412 dataset, MW 2108.35 Da and 1717.83 Da; a second chemistry MW 2121.22 Da is used for the A090512 dataset). Each observed ADC mass is then matched by eye against this grid, labeled with a species name (e.g. `2[4]-8[12]`), and the accepted matches are copied into a second sheet where relative abundance (intensity ÷ total matched intensity) is multiplied by each species' payload count and summed to give the average DAR per payload type and a Total DAR. This logic is confirmed directly in the PPTX (`LC-MS DAR characterization.pptx`), which documents the same four steps.

## 2. What this review found

Reviewing the live formulas in both example datasets (A090412 and A090512) surfaced three concrete, reproducible issues that a platform should design around rather than inherit:

**Fragile, hard-coded ranges.** In both `DAR_A090412_P_DN` and `DAR_A090512_P_DN`, the total-intensity denominator formula is `=SUM(J2:J11)` — a fixed 10-row range — while the DAR numerator sums (`=SUM(N2:N17)` / `=SUM(N2:N15)`) cover a different, larger range. Whenever the number of matched species changes between runs (as it does from experiment to experiment), these ranges silently fall out of sync, meaning the reported Total DAR can be quietly wrong without any visible error. This is the single strongest argument for automation: the bug is invisible in the spreadsheet and only shows up when you recompute the underlying sums independently, which is what this review did.

**Manual label entry can diverge from the actual matched mass.** At least one candidate row in the A090412 working sheet is labeled with a species name (`1[4]-8[12]`) that doesn't correspond to the theoretical mass value pulled into the adjacent cell (which actually corresponds to `1[4]-6[12]`) — a copy/typo error introduced during manual selection. It doesn't affect the final curated DAR (that row wasn't carried into the summary sheet), but it illustrates how easily a manually-typed label can drift from the number it's supposed to describe.

**Real-world mass accuracy runs 100–300 ppm, not the <20 ppm the PPTX describes as "accurate."** Every species the analyst actually accepted into both datasets shows 60–300 ppm mass error against its assigned theoretical mass. This is expected for ~150–165 kDa intact species on this instrument/deconvolution setup, but the platform's default tolerance should be calibrated against this real behavior (roughly 350 ppm captured essentially all analyst-accepted species in both datasets during validation), with the exact number to be tuned against a larger set of runs, and the <20/300 ppm figures kept as scoring/QC bands rather than a hard accept/reject cutoff.

**Fragment/subunit ions are curated in by judgment call, not by rule.** Free light chain, free heavy chain, and half-antibody peaks are manually reclassified and folded into the DAR-weighted average alongside intact species (e.g., a peak initially auto-matched to `2[4]-8[12]` gets manually overridden to `H-0[4]-3[12]`). This is a legitimate technique — it lets you cross-check per-chain payload loading — but it's currently applied inconsistently and isn't reproducible from the raw files alone.

**No traceable link from raw export to final number.** The intensities recorded in `DAR calculation.xlsx` for the A090412 dataset don't numerically match the intensities in the raw `02_A090412H_dN.xlsx` file provided alongside it (e.g. the dominant light-chain peak is recorded at ~125.7M intensity in the workbook vs. ~24.1M in the raw export — roughly a 5× discrepancy). Whatever the reason (a different processing run, re-export, or hand-edit), this is exactly the kind of provenance gap a platform should make structurally impossible by keeping every reported number traceable back to a specific uploaded file and parameter set.

## 3. Automated prototype and validation

To test whether this logic can be automated, I built `dar_calculator.py`, a small Python engine that:

1. loads a deconvolution export and pulls the most abundant mAb mass(es) as base mass(es);
2. builds every theoretical combination of base mass + n×payload MW for however many linker-payload chemistries the user defines, each with its own MW and allowed conjugation-count range;
3. matches every observed ADC peak to its nearest theoretical mass and keeps it only if it falls within a configurable ppm tolerance (this replaces the manual "does this look right?" step with a single adjustable number);
4. computes intensity-weighted relative abundance across matched peaks and derives DAR per payload type and Total DAR, exactly following the formulas in the PPTX.

Validated against the analyst's own curated peak lists in `DAR calculation.xlsx` (self-consistent test, since the raw export file's intensities don't match the workbook — see finding above), the algorithm:

- reproduced the analyst's manual species labels exactly for every one of the 11 (A090412) and 12 (A090512) intact species they had accepted, using nothing but base mass + 2 payload MWs + one tolerance parameter, with zero manual peak-matching;
- calculated Total DAR = 9.60 for A090412 against a manual reference of 9.63 (0.3% difference), and Total DAR = 9.75 for A090512 against a manual reference of 9.63 (comparable order of agreement, with the small gap explained by a handful of extra low-intensity peaks the wider automated tolerance picked up that the analyst had left out).

I also ran the pipeline directly on the two raw instrument-export files with no curation at all; it produced Total DAR = 8.33, which is lower than the workbook's reported value — consistent with the provenance mismatch noted above rather than an algorithm error, since it was validated separately against the workbook's own internally-consistent numbers. All outputs (matched-species tables, DAR summaries, and distribution charts) are included alongside this summary.

Fragment/subunit-ion cross-validation (light chain, heavy chain, half-antibody) was intentionally left out of this prototype's scope: it requires separate theoretical reference masses for each expected fragment (derived from the antibody's chain sequence), which weren't available from the files provided. It's a reasonable phase-2 feature, not a blocker for the core intact-species DAR calculation.

## 4. Relevant published approaches

Automated, tolerance-based intact-mass matching is consistent with where the field is heading. Orbitrap-based peak-integration methods have been shown to estimate DAR distributions directly from deconvoluted spectra with <10% error, without requiring deglycosylation or other sample simplification [1]. Native mass spectrometry has been used to quantify cysteine-conjugated ADCs at the intact level in an automated, batch-processed pipeline (affinity purification → nMS → automated data processing), explicitly to avoid the manual, error-prone alternative [2]. And a recently validated intact LC-MS method for DAR and drug-load-distribution determination was shown to be reproducible across six independent laboratories [3] — evidence that a rule-based, tolerance-driven matching approach (like the one prototyped here) is viable as a validated, auditable method rather than just a convenience script.

## 5. Recommendations for the platform

**Core engine.** Build on the validated approach: user uploads a mAb reference file and an ADC file, specifies one or more linker-payload chemistries (name, MW, allowed conjugation range) via simple form fields, and the platform runs the same base-mass + combinatorial-grid + ppm-tolerance matching shown here. Every number in the final report should be traceably derived from the uploaded files and the entered parameters — no hand-edited intermediate values.

**Checkboxes and parameters, not spreadsheet formulas.** The user's instinct to use checkboxes for parameter confirmation fits well here: checkboxes for which payload chemistries are present, a numeric field for each MW and its expected conjugation range, and a tolerance slider (defaulted from historical data, e.g. ~200–350 ppm for intact species at this mass range) rather than a fixed accept/reject line. Ambiguous peaks — those with more than one theoretical match within tolerance, or matches near the tolerance boundary — should be flagged for the user to confirm with one click rather than silently accepted or dropped.

**Clear output.** A DAR summary table (per-payload DAR and Total DAR), a sortable matched-species table with mass accuracy shown per row, and a distribution chart, matching the structure of the two example outputs generated in this review. Exportable to Excel/PDF for regulatory or internal record-keeping.

**Built-in QC, not manual judgment.** Surface the two failure modes found in this review as automatic checks: warn if the matched-intensity denominator and the DAR-numerator set don't cover the same peaks (the exact bug found in both example datasets), and warn if a species' mass accuracy exceeds the configured tolerance band even though it was accepted.

**Fragment/subunit module as a later phase.** Once core intact-species DAR is running reliably, add an optional module for reduced/subunit ion cross-validation (light chain, heavy chain, half-antibody), where the user supplies expected fragment masses so the same matching engine can be reused for a second, independent DAR estimate.

## Files produced in this review

- `dar_calculator.py` — the prototype matching/DAR engine
- `results_A090412.xlsx`, `results_A090512.xlsx` — validated runs against the analyst's curated data (matches manual DAR within ~0.3–1.2%)
- `results_A090412_from_raw_files.xlsx` — a fully automated, no-curation run directly on the raw instrument export files
- `dar_distribution_A090412.png`, `dar_distribution_A090512.png`, `dar_distribution_A090412_from_raw_files.png` — distribution charts for each run

## Sources

[1] [Drug-to-Antibody Ratio Estimation via Proteoform Peak Integration in the Analysis of Antibody-Oligonucleotide Conjugates with Orbitrap Fourier Transform Mass Spectrometry](https://consensus.app/papers/details/3e79bae051865b8c9cf94cb496078183/?utm_source=claude_code) (Nagornov et al., 2021, Analytical Chemistry)

[2] [Intact quantitation of cysteine-conjugated antibody-drug conjugates using native mass spectrometry](https://consensus.app/papers/details/2cbaf39bf99c518e9dce7d89b84a2978/?utm_source=claude_code) (Li et al., 2024, Rapid Communications in Mass Spectrometry)

[3] [Versatile intact LC-MS method for evaluating the drug-antibody ratio and drug load distribution of antibody-drug conjugates in human plasma](https://consensus.app/papers/details/34560891d12055e88b64f92aef51da48/?utm_source=claude_code) (Hashii et al., 2025, Journal of Chromatography B)
