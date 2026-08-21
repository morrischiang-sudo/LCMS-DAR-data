"""
Synthetic "try it now" example dataset for the SEC LC-MS Analysis page.

Not real instrument data - a small, fabricated deconvolution export built
so a brand-new user can click through the whole app once, with no files of
their own, before deciding whether to trust it with real data. The
generated files and the `PayloadDef`/settings below are built from the same
constants, so re-running the SEC page's analysis on these files with these
settings reproduces a clean, unambiguous result every time.
"""

from __future__ import annotations

import io
from math import comb

import numpy as np
import pandas as pd

from dar_calculator import PayloadDef

# --------------------------------------------------------------------------
# Constants shared between the generated files and the example run settings
# --------------------------------------------------------------------------
NAKED_MAB_MASS = 148_000.0    # Da, a plausible intact IgG1 mass
PAYLOAD_LABEL = "Example payload"
PAYLOAD_MW = 1_000.0          # Da per conjugation event
MAX_COUNT = 8
BINOMIAL_P = 0.45             # skews the fabricated distribution toward DAR ~3-4
PPM_NOISE = 40.0              # small, realistic per-peak mass noise
RNG_SEED = 42

EXAMPLE_PAYLOAD_DEFS = [PayloadDef(label=PAYLOAD_LABEL, mw=PAYLOAD_MW, n_values=list(range(0, MAX_COUNT + 1)))]
EXAMPLE_PPM_TOLERANCE = 150
EXAMPLE_ABUNDANCE_THRESHOLD = 0.0
EXAMPLE_BASE_MASS_MODE = "Most abundant only (recommended)"
EXAMPLE_BASE_MASS_TOP_N = 1


def _to_named_buffer(df: pd.DataFrame, name: str) -> io.BytesIO:
    buf = io.BytesIO()
    df.to_excel(buf, index=False)
    buf.seek(0)
    buf.name = name  # read_uploaded_excel() only reads .name in its error path
    return buf


def build_example_files() -> tuple[io.BytesIO, io.BytesIO]:
    """Build (naked_mab_buffer, adc_buffer) - in-memory .xlsx files with the
    same columns a real BioPharma-Finder-style deconvolution export has.
    """
    rng = np.random.default_rng(RNG_SEED)

    mab_df = pd.DataFrame({
        "Average Mass": [NAKED_MAB_MASS],
        "Sum Intensity": [1.0e9],
        "Relative Abundance": [100.0],
        "Fractional Abundance": [100.0],
    })

    counts = list(range(0, MAX_COUNT + 1))
    weights = np.array([comb(MAX_COUNT, n) * BINOMIAL_P**n * (1 - BINOMIAL_P) ** (MAX_COUNT - n) for n in counts])
    weights = weights / weights.sum()

    total_intensity = 5.0e8
    intensities = weights * total_intensity
    ppm_offsets = rng.uniform(-PPM_NOISE, PPM_NOISE, size=len(counts))
    masses = [
        NAKED_MAB_MASS + n * PAYLOAD_MW + (NAKED_MAB_MASS + n * PAYLOAD_MW) * ppm / 1e6
        for n, ppm in zip(counts, ppm_offsets)
    ]

    adc_df = pd.DataFrame({
        "Average Mass": masses,
        "Sum Intensity": intensities,
    })
    adc_df["Relative Abundance"] = adc_df["Sum Intensity"] / adc_df["Sum Intensity"].sum() * 100
    adc_df["Fractional Abundance"] = adc_df["Relative Abundance"]
    adc_df = adc_df.sort_values("Sum Intensity", ascending=False).reset_index(drop=True)

    return (
        _to_named_buffer(mab_df, "example_naked_mAb.xlsx"),
        _to_named_buffer(adc_df, "example_ADC.xlsx"),
    )
