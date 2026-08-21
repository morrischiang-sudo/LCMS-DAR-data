# AW7310 (RP LC-MS / IdeZ digestion): DAR Compass vs. manual analysis — discrepancy diagnosis

## Bottom line

The P113 discrepancy is real, exactly the size you flagged, and traced to a
single configuration mistake: **P113's mass variant was entered as its
*pair* mass (2971 Da) while the count step was left at 1 (per-molecule).**
Since DAR Compass always treats the configured mass as "Da added per unit of
count," a unit of count came out representing *two* real P113 molecules
instead of one — so every reported P113 count, and DAR[P113], read almost
exactly half its true value. I reproduced this by rebuilding the exact
theoretical grid from your raw files with a guess at your settings, then
re-ran with just that field corrected (and, separately, with the MMAE
breakage variants removed) — the corrected run reproduces your manual
numbers essentially exactly, on real data, not a hand-wave.

| | MMAE mode (%) | P113 mode (%) | DAR [MMAE] | DAR [P113] | Total DAR | Ambiguous |
|---|---|---|---|---|---|---|
| **App, as configured** (P113 = 2971 Da, step 1) | 2 → 56.2% | **3 → 83.2%** | 2.08 | 2.94 | 5.02 | 0 |
| **P113 fixed only** (P113 = 1485.5 Da, step 2; MMAE unchanged) | 2 → 56.2% | **6 → 83.4%** | 2.09 | 5.86 | 7.95 | 0 |
| **P113 fixed + MMAE breakage variants removed** | 2 → 58.6% | 6 → 89.1% | 1.96 | 5.85 | 7.81 | 0 |
| **Your manual reference (AW7310_internal data.xlsx)** | 2 → 58.8% | 6 → 89.1% | 1.96 | 5.88 | 7.83 | n/a |

(Numbers in row 1 are my best reconstruction of your app run from the
description in your message — MMAE at 2460/1532/604 Da, P113 at 2971 Da,
300 ppm — so the exact percentages may differ slightly from your screenshot,
but the pattern that matters — P113's count and DAR reading at almost
exactly half — reproduces cleanly and doesn't depend on getting those
details exactly right.)

## Root cause: P113's mass value and count step don't agree on what "1 unit" means

DAR Compass's rule is simple and the same everywhere in the app: the number
you enter for a chemistry's mass variant is **the Da added by one
conjugated molecule of that payload**, and the count shown/summed for DAR
(`n[label]`) is **the actual number of molecules**. "Allowed count step"
only restricts *which* totals are physically possible (step 2 = only even
totals, for chemistries that conjugate two sites at once) — it does not
change how the mass value is scaled.

Your P113 conjugate only ever shows mass shifts in increments of ~2971 Da
(never ~1485 Da alone) — which is exactly what you'd see if P113 attaches to
both members of a reduced site (e.g. an interchain disulfide) at once, two
molecules per event. That's a real, useful observation. But "2971 Da" is
the mass of *two* P113 molecules, not one. Entering 2971 Da directly as the
MW with step left at 1 means the app's grid treats each unit of count as
2971 Da — so a peak that's actually `6[P113]` (6 real molecules, mass shift
3 × 2971 = 8913 Da) gets computed and labeled as `3[P113]` instead, because
the app is dividing the same 8913 Da by 2971 instead of by 1485.5.

I confirmed this directly on your dominant ADC F(ab')2 peak
(112390.30 Da, 45.5% of raw intensity):

- **As configured** (P113 = 2971 Da, step 1): matches `2[MMAE]-3[P113]`,
  theoretical mass 112385.97 Da (38.5 ppm).
- **P113 corrected** (P113 = 1485.5 Da, step 2): matches
  `2[MMAE]-6[P113]`, theoretical mass **112385.97 Da — identical** (38.5 ppm,
  identical relative abundance ~46.4%).

Same peak, same fit quality, same underlying chemistry call — the only
difference is whether the count label (and therefore the DAR contribution)
says 3 or 6. Every other P113-containing peak follows the same pattern,
which is why DAR[P113] roughly doubles (2.94 → 5.86) once corrected, while
DAR[MMAE] doesn't move at all (2.08 → 2.09) — this is specific to P113's
configuration, not a general matching problem.

**Fix:** in the RP page sidebar, open P113's chemistry expander:
1. Change the mass variant from **2971 Da → 1485.5 Da** (half of what you
   entered — the true per-molecule mass).
2. Set **"Allowed count step" → 2** for both the F(ab')2 and Fc conjugation
   ranges (it's configured per fragment) — this keeps the grid restricted to
   even totals only (0, 2, 4, 6, 8), matching the real chemistry and
   avoiding any risk of a spurious match at an odd count.

Do **not** leave the mass at 2971 Da even with step set to 2 — that
double-counts (n × 2971 instead of n × 1485.5) and would overstate
DAR[P113] by a further 2×.

## Secondary, smaller effect: MMAE's two breakage variants blur its distribution slightly

Your manual reference file only ever labels MMAE species as plain `n[3]` —
it doesn't appear to model the 1532 Da / 604 Da breakage masses you
mentioned at all. With those two variants included, DAR Compass's MMAE mode
comes out at 56.2% (count 2) instead of your reference's 58.8%; removing
them (single 2460 Da mass only) reproduces 58.6% — matching your reference
to within 0.2 points, and also nudges P113 up to 89.1% (from 83.4%),
matching your reference's 89.1% almost exactly.

This is a much smaller effect than the P113 issue (a couple of percentage
points, not a 2× factor) and isn't necessarily wrong to include — if you
have independent evidence that MMAE genuinely fragments to 1532/604 Da
under your sample prep, keep them and expect a small, known blur relative
to a simpler intact-only model. If that evidence is soft, dropping down to
a single 2460 Da mass reproduces your manual numbers most closely, the same
conclusion as the earlier ACE723 diagnosis for a different multi-variant
chemistry.

## What to change in the app

1. P113 chemistry → mass variant: **2971 Da → 1485.5 Da**.
2. P113 chemistry → **Allowed count step**: 1 → **2** (both F(ab')2 and Fc
   conjugation-range sections, since they're independent per fragment).
3. Optional, smaller effect: consider dropping MMAE's 1532/604 Da breakage
   variants back to a single 2460 Da mass unless you have specific evidence
   for that fragmentation under your RP sample prep.
4. Re-run. Expect P113's drug-load mode to land at count 6 (not 3) at
   ~83-89% depending on whether you keep the MMAE breakage variants, and
   Total DAR (MMAE + P113, reported separately, no combined figure) close
   to your manual 1.96 + 5.88 = 7.83.

No code changes are needed — this is a configuration issue specific to how
P113's mass-per-count was entered, not a bug in the matching/DAR engine.
