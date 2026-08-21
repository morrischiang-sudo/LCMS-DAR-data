"""
dar_calculator.py
------------------
Prototype automation engine for ADC Drug-to-Antibody Ratio (DAR) analysis
from LC-MS deconvoluted mass lists (e.g. Thermo BioPharma Finder exports).

This reproduces, as an automated/reproducible algorithm, the manual Excel
workflow documented in "DAR calculation.xlsx" and "LC-MS DAR
characterization.pptx":

  1. Take the naked mAb (reference) deconvolution export and use its
     observed intact masses as "base masses" (these capture glycoform /
     adduct heterogeneity of the unconjugated antibody).
  2. For each linker-payload chemistry the user specifies (name, one or
     more possible masses per conjugation event, and an allowed
     conjugation-number range), build the full combinatorial grid of
     theoretical intact masses. Each chemistry may have multiple "mass
     variants" (e.g. an intact linker-payload mass and one or more
     lower-mass forms from breakage during harsh sample processing);
     variants are mixed independently per attachment site within one
     molecule, since breakage is expected to act per-site rather than
     uniformly across a whole molecule:
         theoretical_mass = base_mass + sum over chemistries of
                             sum over that chemistry's variants of
                             (n_variant * MW_variant)
  3. Match each observed ADC deconvolution peak to its nearest theoretical
     mass. Keep the match only if the mass accuracy is within a ppm
     tolerance (the PPTX states <20 ppm = accurate, ~300 ppm = inaccurate;
     exposed here as a user-adjustable parameter).
  4. Compute intensity-weighted relative abundance across all matched
     species, then per-payload DAR = sum(relative_abundance * n_variant *
     dar_weight_variant), and Total DAR = sum across payload chemistries.
     `dar_weight` lets a "broken" mass variant count toward DAR the same
     as an intact one (default), partially, or not at all (0.0) - this is
     a chemistry-specific judgment call, not something this module decides
     on its own.

This module is meant to be the computational core behind an upload-and-click
web UI: user uploads mAb + ADC deconvolution files, ticks which linker-payload
chemistries apply and enters their mass variant(s)/valence, and gets back a
DAR report.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from math import comb
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Data loading
# --------------------------------------------------------------------------

def load_deconvolution_file(path: str | Path, sheet_name=0) -> pd.DataFrame:
    """Load a BioPharma-Finder-style deconvolution export (.xlsx).

    Expected columns include at least:
        'Average Mass', 'Sum Intensity', 'Relative Abundance',
        'Fractional Abundance'
    """
    df = pd.read_excel(path, sheet_name=sheet_name)
    required = {"Average Mass", "Sum Intensity"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path}: missing expected column(s): {missing}")
    return df


def base_masses_from_mab(
    mab_df: pd.DataFrame,
    top_n: int | None = None,
    min_fractional_abundance: float | None = None,
) -> list[float]:
    """Return candidate unconjugated-mAb masses (glycoforms/adducts),
    ranked by abundance, to use as the base of the theoretical mass ladder.
    """
    df = mab_df.copy()
    abundance_col = "Fractional Abundance" if "Fractional Abundance" in df.columns else "Relative Abundance"
    df = df.sort_values(abundance_col, ascending=False)
    if min_fractional_abundance is not None:
        df = df[df[abundance_col] >= min_fractional_abundance]
    if top_n is not None:
        df = df.head(top_n)
    return df["Average Mass"].tolist()


def filter_by_fractional_abundance(
    df: pd.DataFrame,
    min_fractional_abundance: float = 0.0,
    abundance_col: str | None = None,
) -> pd.DataFrame:
    """Keep only rows at or above a fractional-abundance threshold.

    This is the automated equivalent of an analyst visually skipping very
    low-abundance/noise peaks when picking "top candidates" out of a raw
    deconvolution export before matching them against the theoretical mass
    ladder. Typically applied to the ADC file (candidate species selection)
    but works on any deconvolution dataframe.

    A threshold of 0 (the default) is a no-op and returns the input
    unchanged, so existing callers behave exactly as before.
    """
    if min_fractional_abundance is None or min_fractional_abundance <= 0:
        return df

    col = abundance_col
    if col is None:
        col = "Fractional Abundance" if "Fractional Abundance" in df.columns else "Relative Abundance"
    if col not in df.columns:
        raise ValueError(
            f"Cannot filter by fractional abundance: column '{col}' not found "
            f"(available columns: {list(df.columns)})."
        )
    return df[df[col] >= min_fractional_abundance].reset_index(drop=True)


# --------------------------------------------------------------------------
# Theoretical mass grid
# --------------------------------------------------------------------------

@dataclass
class MassVariant:
    """One possible per-site mass for a linker-payload chemistry.

    Most chemistries need only one variant (the intact linker-payload mass).
    Add more when breakage during harsh sample processing produces
    additional possible masses at a conjugation site (e.g. a lower-mass
    degraded form).
    """
    mw: float                      # Da, added mass per conjugation event for this variant
    dar_weight: float = 1.0        # how much one site of this variant counts toward DAR
                                    # (1.0 = counts as a full payload, same as today;
                                    #  0.0 = doesn't count, e.g. total payload loss)
    variant_label: str = ""        # filled in automatically by PayloadDef if left blank


@dataclass
class PayloadDef:
    """One linker-payload chemistry the user specifies for the platform.

    Backward compatible: existing calls like
        PayloadDef(label="4", mw=2108.35, n_values=range(0, 6))
    keep working unchanged (a single implicit MassVariant with dar_weight=1.0).
    To add breakage-derived mass variants, pass `variants` instead:
        PayloadDef(label="MMAE", n_values=range(0, 5), variants=[
            MassVariant(mw=2460.88),                  # intact
            MassVariant(mw=1718.85, dar_weight=0.0),  # payload fully lost
        ])
    """
    label: str                                    # e.g. "4", "12", "DXd", "MMAE"
    mw: float | None = None                       # Da; ignored if `variants` is given
    n_values: Sequence[int] = field(default_factory=lambda: range(0, 9))  # allowed TOTAL counts
    variants: Sequence[MassVariant] | None = None  # one or more possible per-site masses

    def __post_init__(self):
        if self.variants:
            variants = list(self.variants)
        elif self.mw is not None:
            variants = [MassVariant(mw=float(self.mw), dar_weight=1.0)]
        else:
            raise ValueError(f"PayloadDef '{self.label}': provide either `mw` or `variants`.")

        if len(variants) == 1:
            if not variants[0].variant_label:
                variants[0].variant_label = self.label
        else:
            for i, v in enumerate(variants):
                if not v.variant_label:
                    v.variant_label = self.label if i == 0 else f"{self.label}_b{i}"

        self.variants = variants


def _partitions(n: int, k: int):
    """Yield every k-tuple of non-negative ints summing to n.

    This is how one chemistry's total occupied-site count (n) is split
    across its mass variants, since breakage is modeled as acting
    independently per attachment site (a single molecule can have some
    sites intact and others broken at the same time).
    """
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for rest in _partitions(n - i, k - 1):
            yield (i,) + rest


def _chemistry_label(combo: tuple[int, ...], variants: Sequence[MassVariant]) -> str:
    if len(variants) == 1:
        return f"{combo[0]}[{variants[0].variant_label}]"
    nonzero = [(c, v) for c, v in zip(combo, variants) if c > 0]
    if not nonzero:
        return f"0[{variants[0].variant_label}]"
    return "+".join(f"{c}[{v.variant_label}]" for c, v in nonzero)


def estimate_theoretical_grid_size(base_masses: Sequence[float], payload_defs: Sequence[PayloadDef]) -> int:
    """Roughly how many theoretical species build_theoretical_table will produce.

    Useful as a cheap sanity check before running a large multi-variant
    grid - the count grows quickly with more mass variants per chemistry.
    """
    total = max(len(base_masses), 1)
    for p in payload_defs:
        k = len(p.variants)
        chem_count = 0
        for n in p.n_values:
            chem_count += comb(n + k - 1, k - 1) if k > 1 else 1
        total *= max(chem_count, 1)
    return total


def build_theoretical_table(base_masses: Sequence[float], payload_defs: Sequence[PayloadDef]) -> pd.DataFrame:
    """Build every combination of (base mass) x (per-chemistry variant mixture).

    For a chemistry with a single mass variant, this reduces exactly to the
    original behavior: one count `n` per chemistry, label "n[label]".

    For a chemistry with multiple mass variants, every total count `n` in
    `n_values` is split every possible way across that chemistry's variants
    (e.g. n=3 with 2 variants -> (3,0), (2,1), (1,2), (0,3)), modeling
    breakage as an independent per-site event within one molecule.
    """
    per_chem_options = []
    for p in payload_defs:
        variants = p.variants
        options = []
        for n in p.n_values:
            for combo in _partitions(n, len(variants)):
                mass_contrib = sum(c * v.mw for c, v in zip(combo, variants))
                label_piece = _chemistry_label(combo, variants)
                n_by_variant = {v.variant_label: c for c, v in zip(combo, variants)}
                options.append({
                    "total_n": n,
                    "mass_contrib": mass_contrib,
                    "label_piece": label_piece,
                    "n_by_variant": n_by_variant,
                })
        per_chem_options.append(options)

    rows = []
    for base in base_masses:
        for combo_opts in itertools.product(*per_chem_options):
            mass = base + sum(o["mass_contrib"] for o in combo_opts)
            label = "-".join(o["label_piece"] for o in combo_opts)
            row = {"base_mass": base, "theoretical_mass": mass, "label": label}
            for p, o in zip(payload_defs, combo_opts):
                for var_label, c in o["n_by_variant"].items():
                    row[f"n_{var_label}"] = c
                if len(p.variants) > 1:
                    row[f"n_{p.label}_total"] = o["total_n"]
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Matching
# --------------------------------------------------------------------------

def match_species(
    observed_df: pd.DataFrame,
    theoretical_df: pd.DataFrame,
    ppm_tolerance: float = 150.0,
    mass_col: str = "Average Mass",
) -> pd.DataFrame:
    """Match every observed peak to its closest theoretical species.

    A peak is KEPT only if the closest theoretical mass is within
    ppm_tolerance; otherwise it is treated as a fragment/adduct/noise peak
    and dropped from the intact-species DAR calculation (this is the
    automated replacement for the manual "is this peak real?" judgment call).

    Each kept match is also checked for ambiguity: if a *different* species
    label also falls within ppm_tolerance, the match is flagged
    (``ambiguous=True``) with the runner-up species and its ppm error, so a
    UI can surface it for human review instead of silently picking the
    nearest one. This directly automates the "does this look right?" step
    that was previously a manual judgment call.
    """
    theo_masses = theoretical_df["theoretical_mass"].to_numpy()
    theo_labels = theoretical_df["label"].to_numpy()
    n_cols = [c for c in theoretical_df.columns if c.startswith("n_")]

    matched_rows = []
    for _, row in observed_df.iterrows():
        om = row[mass_col]
        ppm_errors = np.abs((theo_masses - om) / theo_masses) * 1e6
        order = np.argsort(ppm_errors)
        best_idx = int(order[0])
        if ppm_errors[best_idx] > ppm_tolerance:
            continue

        trow = theoretical_df.iloc[best_idx]
        rec = row.to_dict()
        rec["theoretical_mass"] = trow["theoretical_mass"]
        rec["species"] = trow["label"]
        rec["ppm_error"] = ppm_errors[best_idx]
        for c in n_cols:
            rec[c] = trow[c]

        # Ambiguity check: does a different-labeled species also fall within tolerance?
        rec["ambiguous"] = False
        rec["runner_up_species"] = None
        rec["runner_up_ppm_error"] = None
        for idx2 in order[1:]:
            err2 = ppm_errors[idx2]
            if err2 > ppm_tolerance:
                break  # order is sorted ascending, nothing further can qualify
            if theo_labels[idx2] != trow["label"]:
                rec["ambiguous"] = True
                rec["runner_up_species"] = theo_labels[idx2]
                rec["runner_up_ppm_error"] = err2
                break

        matched_rows.append(rec)

    matched = pd.DataFrame(matched_rows)
    return matched.sort_values("Sum Intensity", ascending=False).reset_index(drop=True) if len(matched) else matched


def build_verification_table(
    observed_df: pd.DataFrame,
    theoretical_df: pd.DataFrame,
    ppm_tolerance: float,
    mass_col: str = "Average Mass",
    abundance_threshold: float = 0.0,
    abundance_col: str | None = None,
) -> pd.DataFrame:
    """Full audit trail: every peak in `observed_df` (matched or not), with
    its closest theoretical candidate, delta mass, ppm error, and pass/fail
    flags - so a researcher can verify (or debug) exactly why any given peak
    was or wasn't counted, including peaks excluded by the fractional
    abundance threshold before matching even started.

    Unlike `match_species`, nothing here is filtered out: a peak whose
    closest theoretical mass is thousands of ppm away still gets a row,
    showing that value plainly. This is deliberate - it's what makes it
    possible to notice, for example, that a configured chemistry's max
    conjugation count was set too low to reach a peak's true species,
    rather than silently dropping that peak with no trace.
    """
    if abundance_col is None:
        abundance_col = "Fractional Abundance" if "Fractional Abundance" in observed_df.columns else "Relative Abundance"

    theo_masses = theoretical_df["theoretical_mass"].to_numpy()
    theo_labels = theoretical_df["label"].to_numpy()

    rows = []
    for _, row in observed_df.iterrows():
        om = row[mass_col]
        ppm_errors = np.abs((theo_masses - om) / theo_masses) * 1e6
        idx = int(np.argmin(ppm_errors))

        rec = row.to_dict()
        rec["closest_theoretical_mass"] = theo_masses[idx]
        rec["closest_species"] = theo_labels[idx]
        rec["delta_mass"] = om - theo_masses[idx]
        rec["ppm_error"] = ppm_errors[idx]

        abundance_value = row.get(abundance_col)
        rec["passed_abundance_threshold"] = bool(abundance_value is None or abundance_value >= abundance_threshold)
        rec["within_ppm_tolerance"] = bool(ppm_errors[idx] <= ppm_tolerance)
        rec["matched"] = bool(rec["passed_abundance_threshold"] and rec["within_ppm_tolerance"])
        rows.append(rec)

    return pd.DataFrame(rows).sort_values(mass_col).reset_index(drop=True)


def build_selection_summary(verification_df: pd.DataFrame, abundance_col: str | None = None) -> pd.DataFrame:
    """Quick funnel overview: how many ADC peaks - and how much of the total
    signal - survive each stage (fractional abundance threshold, then ppm
    matching).

    Peak *count* and *intensity* percentages are shown side by side because
    they can tell very different stories: matching only a small fraction of
    peaks by count can still mean capturing nearly all of the real signal,
    since low-abundance noise peaks are expected to go unmatched. Built
    directly from `build_verification_table`'s output, so it always agrees
    with what that table shows.
    """
    if verification_df.empty:
        return pd.DataFrame()

    if abundance_col is None:
        abundance_col = "Fractional Abundance" if "Fractional Abundance" in verification_df.columns else "Relative Abundance"

    n_total = len(verification_df)
    passed = verification_df["passed_abundance_threshold"]
    matched = verification_df["matched"]
    unmatched_candidates = passed & ~matched

    total_abundance = verification_df[abundance_col].sum() if abundance_col in verification_df.columns else None

    def pct_of_abundance(mask):
        if total_abundance in (None, 0):
            return None
        return float(verification_df.loc[mask, abundance_col].sum() / total_abundance * 100)

    stages = [
        ("Total peaks in ADC file", pd.Series([True] * n_total, index=verification_df.index)),
        ("Excluded by fractional abundance threshold", ~passed),
        ("Candidate peaks (passed threshold)", passed),
        ("Matched to a theoretical species", matched),
        ("Unmatched (passed threshold, no species within tolerance)", unmatched_candidates),
    ]

    rows = []
    for name, mask in stages:
        count = int(mask.sum())
        rows.append({
            "Stage": name,
            "Peak count": count,
            "% of peaks": round(count / n_total * 100, 1) if n_total else None,
            "% of total signal (Fractional Abundance)": (
                round(pct_of_abundance(mask), 1) if pct_of_abundance(mask) is not None else None
            ),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# DAR calculation
# --------------------------------------------------------------------------

def calculate_dar(matched_df: pd.DataFrame, payload_defs: Sequence[PayloadDef], intensity_col: str = "Sum Intensity") -> tuple[dict, pd.DataFrame]:
    """Intensity-weighted DAR per payload chemistry + total, following the
    formulas in the PPTX, generalized for multiple mass variants per
    chemistry:
        relative_abundance(species) = intensity(species) / sum(intensity, matched species)
        DAR(chemistry) = sum over matched species of
                          relative_abundance(species) *
                          sum over that chemistry's variants of (n_variant(species) * dar_weight_variant)
        Total DAR      = sum(DAR(chemistry) for all chemistries)

    A chemistry with a single mass variant (dar_weight=1.0, the default)
    reduces exactly to the original formula: DAR = sum(relative_abundance * n_payload).
    """
    if matched_df.empty:
        return {**{p.label: 0.0 for p in payload_defs}, "total": 0.0}, matched_df

    out = matched_df.copy()
    total_intensity = out[intensity_col].sum()
    out["relative_abundance"] = out[intensity_col] / total_intensity

    dar = {}
    for p in payload_defs:
        contrib = pd.Series(0.0, index=out.index)
        for v in p.variants:
            col = f"n_{v.variant_label}"
            if col in out.columns:
                contrib = contrib + out[col] * v.dar_weight
        out[f"dar_contrib_{p.label}"] = out["relative_abundance"] * contrib
        dar[p.label] = out[f"dar_contrib_{p.label}"].sum()
    dar["total"] = float(sum(dar.values()))
    return dar, out


def consolidate_by_total_count(
    matched_df: pd.DataFrame,
    payload_defs: Sequence[PayloadDef],
    intensity_col: str = "Sum Intensity",
) -> pd.DataFrame:
    """Collapse matched species down to one row per unique combination of
    TOTAL occupied-site count per chemistry, discarding which specific mass
    variant(s) made up that total.

    This exists purely to make a simplified chart/summary possible once a
    chemistry has multiple mass variants (e.g. intact + a breakage product):
    without it, every different intact/broken split at the same total count
    shows up as its own bar, which gets cluttered fast. This does NOT change
    any DAR number - it only re-groups the already-matched, already-weighted
    rows for a cleaner display. Run `calculate_dar` first; pass its output
    (which already has `relative_abundance` and `ambiguous` columns) in here.

    For chemistries with a single mass variant, the "total count" is just
    that chemistry's regular n_<label> column, so this is a no-op relabeling
    in the common case (every group has exactly one row).
    """
    if matched_df.empty:
        return matched_df

    total_cols = []
    for p in payload_defs:
        col = f"n_{p.label}_total" if len(p.variants) > 1 else f"n_{p.variants[0].variant_label}"
        total_cols.append((p.label, col))

    df = matched_df.copy()
    group_cols = [col for _, col in total_cols]
    agg_kwargs = {intensity_col: (intensity_col, "sum")}
    if "relative_abundance" in df.columns:
        agg_kwargs["relative_abundance"] = ("relative_abundance", "sum")
    if "ambiguous" in df.columns:
        agg_kwargs["ambiguous"] = ("ambiguous", "any")
    grouped = df.groupby(group_cols, as_index=False).agg(**agg_kwargs)

    grouped["species"] = grouped.apply(
        lambda row: "-".join(f"{int(row[col])}[{label}]" for label, col in total_cols),
        axis=1,
    )
    return grouped.sort_values(intensity_col, ascending=False).reset_index(drop=True)


def marginal_distribution_by_chemistry(
    matched_df: pd.DataFrame,
    payload_defs: Sequence[PayloadDef],
    intensity_col: str = "Sum Intensity",
) -> pd.DataFrame:
    """Per-chemistry drug-load distribution: the % of total matched
    intensity at each total occupied-site count, for one chemistry at a
    time, independent of every other chemistry's count.

    This is the standard "drug-load distribution" report format (one row
    per linker-payload, one column per count 0..max, plus an average) -
    distinct from `consolidate_by_total_count`, which groups by the joint
    combination of every chemistry's count together rather than one
    chemistry at a time.

    Every count configured in a chemistry's `n_values` gets an explicit row
    (0.0 if nothing matched there), so a count that's in range but simply
    unobserved still shows as 0, while a count outside that chemistry's
    configured range doesn't appear at all - matching how a report table
    would leave truly out-of-range cells blank rather than showing 0.

    Returns a long-format dataframe: columns
    ['chemistry', 'count', 'relative_abundance_pct'].
    """
    rows = []
    for p in payload_defs:
        col = f"n_{p.label}_total" if len(p.variants) > 1 else f"n_{p.variants[0].variant_label}"
        if matched_df.empty or col not in matched_df.columns:
            grouped = {}
        else:
            grouped = matched_df.groupby(col)["relative_abundance"].sum().to_dict()
        for n in p.n_values:
            rows.append({
                "chemistry": p.label,
                "count": int(n),
                "relative_abundance_pct": grouped.get(n, 0.0) * 100,
            })
    return pd.DataFrame(rows)


def build_drug_load_summary_table(
    matched_df: pd.DataFrame,
    payload_defs: Sequence[PayloadDef],
    dar: dict,
) -> pd.DataFrame:
    """Wide-format drug-load distribution table: one row per chemistry,
    one column per count, plus an 'Average' column.

    'Average' is taken directly from `calculate_dar`'s DAR value for that
    chemistry (the same number shown in the app's top-line metrics), not
    recomputed from this table's own counts - the two agree exactly when
    every mass variant's `dar_weight` is 1.0 (the default), and can differ
    when a variant has a lower/zero DAR weight, since this table's counts
    are total occupied sites regardless of variant, while DAR discounts
    variants that don't fully count as a payload.
    """
    long_df = marginal_distribution_by_chemistry(matched_df, payload_defs)
    if long_df.empty:
        return pd.DataFrame()

    pivot = long_df.pivot_table(index="chemistry", columns="count", values="relative_abundance_pct")
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)
    pivot = pivot.reindex([p.label for p in payload_defs])
    pivot.columns = [str(c) for c in pivot.columns]  # numeric-sorted first, then stringified for display/export
    pivot["Average"] = [dar.get(label, float("nan")) for label in pivot.index]
    return pivot


# --------------------------------------------------------------------------
# High-level convenience wrapper
# --------------------------------------------------------------------------

def run_dar_analysis(
    mab_path: str | Path,
    adc_path: str | Path,
    payload_defs: Sequence[PayloadDef],
    ppm_tolerance: float = 150.0,
    base_mass_top_n: int | None = None,
    adc_min_fractional_abundance: float = 0.0,
) -> dict:
    mab_df = load_deconvolution_file(mab_path)
    adc_df_all = load_deconvolution_file(adc_path)
    adc_df = filter_by_fractional_abundance(adc_df_all, adc_min_fractional_abundance)

    base_masses = base_masses_from_mab(mab_df, top_n=base_mass_top_n)
    theoretical = build_theoretical_table(base_masses, payload_defs)
    matched = match_species(adc_df, theoretical, ppm_tolerance=ppm_tolerance)
    dar, matched_with_contrib = calculate_dar(matched, payload_defs)

    return {
        "base_masses": base_masses,
        "theoretical_table": theoretical,
        "matched_species": matched_with_contrib,
        "n_observed_peaks": len(adc_df_all),
        "n_candidate_peaks": len(adc_df),
        "n_matched_peaks": len(matched_with_contrib),
        "dar": dar,
    }


if __name__ == "__main__":
    import sys
    print(__doc__)
