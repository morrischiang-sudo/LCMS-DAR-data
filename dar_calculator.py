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
  2. For each linker-payload chemistry the user specifies (name, monoisotopic/
     average MW, allowed conjugation-number range), build the full
     combinatorial grid of theoretical intact masses:
         theoretical_mass = base_mass + sum(n_i * MW_i)
  3. Match each observed ADC deconvolution peak to its nearest theoretical
     mass. Keep the match only if the mass accuracy is within a ppm
     tolerance (the PPTX states <20 ppm = accurate, ~300 ppm = inaccurate;
     exposed here as a user-adjustable parameter).
  4. Compute intensity-weighted relative abundance across all matched
     species, then per-payload DAR = sum(relative_abundance * n_payload),
     and Total DAR = sum across payload types.

This module is meant to be the computational core behind an upload-and-click
web UI: user uploads mAb + ADC deconvolution files, ticks which linker-payload
chemistries apply and enters their MW/valence, and gets back a DAR report.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
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
class PayloadDef:
    """One linker-payload chemistry the user specifies for the platform."""
    label: str            # e.g. "4", "12", "DXd", "MMAE"
    mw: float              # Da, added mass per conjugation event
    n_values: Sequence[int] = field(default_factory=lambda: range(0, 9))  # allowed counts


def build_theoretical_table(base_masses: Sequence[float], payload_defs: Sequence[PayloadDef]) -> pd.DataFrame:
    """Build every combination of (base mass) x (n per payload type)."""
    rows = []
    n_value_lists = [list(p.n_values) for p in payload_defs]
    for base in base_masses:
        for combo in itertools.product(*n_value_lists):
            mass = base + sum(n * p.mw for n, p in zip(combo, payload_defs))
            label = "-".join(f"{n}[{p.label}]" for n, p in zip(combo, payload_defs))
            row = {"base_mass": base, "theoretical_mass": mass, "label": label}
            for n, p in zip(combo, payload_defs):
                row[f"n_{p.label}"] = n
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


# --------------------------------------------------------------------------
# DAR calculation
# --------------------------------------------------------------------------

def calculate_dar(matched_df: pd.DataFrame, payload_labels: Sequence[str], intensity_col: str = "Sum Intensity") -> tuple[dict, pd.DataFrame]:
    """Intensity-weighted DAR per payload type + total, following the
    formulas in the PPTX:
        relative_abundance(species) = intensity(species) / sum(intensity, matched species)
        DAR(payload)  = sum(relative_abundance * n_payload)
        Total DAR     = sum(DAR(payload) for all payload types)
    """
    if matched_df.empty:
        return {**{lbl: 0.0 for lbl in payload_labels}, "total": 0.0}, matched_df

    out = matched_df.copy()
    total_intensity = out[intensity_col].sum()
    out["relative_abundance"] = out[intensity_col] / total_intensity

    dar = {}
    for lbl in payload_labels:
        col = f"n_{lbl}"
        out[f"dar_contrib_{lbl}"] = out["relative_abundance"] * out[col]
        dar[lbl] = out[f"dar_contrib_{lbl}"].sum()
    dar["total"] = float(sum(dar.values()))
    return dar, out


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
    dar, matched_with_contrib = calculate_dar(matched, [p.label for p in payload_defs])

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
