# ACE723: DAR Compass vs. Manual Analysis — Discrepancy Diagnosis

## Bottom line

The gap is real but almost entirely explained by one setting: **chemistry "M"'s
"Max conjugation count" was set to 1** in the app run behind `dar_report (2).xlsx`,
when the actual sample needs it to go up to at least 4. I reproduced your exact
app output byte-for-byte with that setting, then re-ran with just that one field
corrected — Total DAR jumps from 8.77 to 9.98, within 1% of your manual 10.08,
with **zero ambiguous matches**. Three smaller, secondary factors explain the
remaining gap. None of this points to a bug in the matching/DAR math itself —
the engine is doing exactly what it's told; it was told the wrong thing for M.

| | DAR [M] | DAR [D] | Total DAR | Matched peaks | Ambiguous |
|---|---|---|---|---|---|
| **App, as run** (Max conjugation count for M = 1) | 0.99 | 7.78 | 8.77 | 5 / 45 | 0 |
| **App, with M's max fixed to 5** (single mAb base mass, no breakage variant) | 2.01 | 7.98 | 9.98 | 10 / 45 | 0 |
| **Your manual analysis** | 2.16 | 7.92 | 10.08 | 20 / 45 (17 systematic + 3 hand-typed) | n/a |

## Root cause 1 (primary): Max conjugation count for M = 1

Your raw ACE723 cAb file's single most abundant peak (157100.96 Da, 35.0% of
total intensity) is `2[M]-8[D]` — it needs **two** copies of M. Every peak
needing 2, 3, or 4 copies of M (which is most of the high-abundance signal)
has no nearby theoretical mass to match against once M's range is capped at
`{0, 1}`, so it's silently dropped — not flagged as an error, just absent.
That's why only 5 of 45 peaks matched, and why DAR [M] came out at roughly
half its real value: the calculation is only ever allowed to see M-loads of
0 or 1.

I confirmed this by rebuilding the exact grid from your two raw files with
M capped at 1 and the breakage variant you'd configured (MW ≈ 993.17, see
Root cause 2) — it reproduces the app's actual 5 matched species, ppm
errors, and DAR values exactly. Re-running with only that field changed
(Max conjugation count: 1 → 5) jumps matched peaks to 10 and Total DAR to
9.98.

**Fix:** in the sidebar, open chemistry M's expander and set **Max
conjugation count** to at least 4 (5 to match the range you used manually).

## Root cause 2: the M "breakage" mass variant doesn't match how you're handling it manually

The app run had a second mass variant on M (MW ≈ 993.17 Da, a ~762 Da loss
from intact M's 1755.17 Da). Your manual analysis handles apparent M-loss
differently: three peaks (154595.79, 156350.39, 158100.47 Da — together
~13.4% of matched fractional abundance) are labeled `...(-M)` with
hand-typed theoretical masses close to the *observed* mass itself (implied
loss ≈ 725–735 Da, not 762 Da, and not applied systematically). These are
one-off expert calls, not a fixed second mass applied everywhere M appears —
they're a different kind of judgment than a configured mass variant, so the
app can't and won't reproduce them by matching MW values more precisely.

With Max conjugation count corrected to 5, re-adding that same breakage
variant matches more peaks (28 vs. 10) but pushes DAR [M] to 2.79 and DAR
[D] down to 6.93 — *further* from your manual 2.16 / 7.92, and introduces 2
ambiguous flags. In this dataset, the variant doesn't reproduce your manual
`(-M)` calls; it mostly finds lower-confidence alternative matches for peaks
that already had a fine explanation without it.

**Recommendation:** turn the M breakage variant off (set "Number of mass
variants" back to 1) unless you have independent evidence pinning down a
specific, systematic breakage mass. Treat the handful of outlier peaks
case-by-case, the way you already do manually — that's a legitimate
approach the app doesn't try to replace (see the "known limitations" note
in the README about fragment/subunit calls).

## Root cause 3 (minor): one vs. two mAb base masses

Your manual grid anchors off **two** mAb reference masses — 145290.97 Da
(the dominant glycoform) and 145465.71 Da (a second glycoform at 14.8%
fractional abundance) — while the app defaults to "Most abundant only"
(one base mass). About 40% of your manually-accepted systematic matches
(7 of 17, excluding the 3 hand-typed ones) resolve only through the second
base mass.

In this specific dataset, switching the app to "Top N glycoform/adduct
variants" with N=2 didn't clearly help once M's max count was already
fixed — it matched more peaks (19 vs. 10) but drifted DAR [M] to 2.30 and
Total DAR to 9.60, slightly further from your manual number, not closer.
Worth knowing about, but not the lever to pull here — M's max conjugation
count is.

## Root cause 4: three fragment calls with no systematic equivalent

The three `(-M)` peaks mentioned above use theoretical masses typed by hand
per-peak, not derived from any `base + n×MW` formula. No setting change
makes the app reproduce these — they're exactly the kind of ad hoc,
expert-judgment fragment call that's outside what automated combinatorial
matching does. This is very likely most of the residual ~1% gap between the
corrected app run (9.98) and your full manual total (10.08).

## What to change in the app

1. Chemistry M → **Max conjugation count**: 1 → 5 (or at least 4).
2. Chemistry M → **Number of mass variants**: 2 → 1 (drop the breakage
   variant, at least as a starting point).
3. Leave **mAb base mass(es) to use** at "Most abundant only" — the second
   base mass isn't earning its complexity here.
4. Re-run. Expect ~9-10 matched species, zero or near-zero ambiguous flags,
   and a Total DAR close to 9.98 — within about 1% of your manual result,
   with the remainder attributable to the three hand-called fragment peaks
   this method doesn't try to capture.
