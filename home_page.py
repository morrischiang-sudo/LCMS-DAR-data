"""
Home / chooser page.

First thing a new user sees: two large, descriptive cards explaining when to
use each workflow, with a direct link into each. `app.py` builds the actual
st.Page objects (st.page_link needs the same Page object that was registered
with st.navigation for callable-based pages) and hands them to this module
via st.session_state before calling st.navigation.
"""

from __future__ import annotations

from pathlib import Path

import streamlit as st

from dar_ui_helpers import APP_TAGLINE, APP_TITLE, run_fragment_analysis
from example_data import (
    EXAMPLE_ABUNDANCE_THRESHOLD,
    EXAMPLE_BASE_MASS_MODE,
    EXAMPLE_BASE_MASS_TOP_N,
    EXAMPLE_PAYLOAD_DEFS,
    EXAMPLE_PPM_TOLERANCE,
    build_example_files,
)

_ASSET_DIR = Path(__file__).parent
SEC_CARTOON = _ASSET_DIR / "SEC LC-MS section cartoon.png"
IDEZ_CARTOON = _ASSET_DIR / "IdelZ section cartoon.png"


def render() -> None:
    st.title(f"\U0001F9EA {APP_TITLE}")
    st.caption(APP_TAGLINE)
    st.markdown("#### Two independent workflows — pick the one that matches your sample prep")

    pages = st.session_state.get("_pages", {})
    sec_page_obj = pages.get("sec")
    rp_page_obj = pages.get("rp")

    col1, col2 = st.columns(2, gap="large")

    with col1:
        with st.container(border=True):
            st.markdown("### \U0001F9EA SEC LC-MS Analysis")
            st.markdown(
                "Intact-mass DAR analysis for a **whole, undigested ADC**. Upload one naked-mAb "
                "reference file and one ADC file; DAR is calculated directly from the intact "
                "mass shift."
            )
            st.markdown(
                "**Use this when:** your sample wasn't enzymatically digested, or you don't "
                "need to know *where* on the antibody the payload sits — just the overall DAR."
            )
            if sec_page_obj is not None:
                st.page_link(sec_page_obj, label="Open SEC LC-MS Analysis →", icon="\U0001F9EA", use_container_width=True)
        if SEC_CARTOON.exists():
            st.image(
                str(SEC_CARTOON),
                caption="SEC separates intact ADC species by drug load — higher-DAR species elute earlier.",
                use_container_width=True,
            )

    with col2:
        with st.container(border=True):
            st.markdown("### \U0001F52C RP LC-MS (IdeZ digestion)")
            st.markdown(
                "DAR analysis for samples treated with **IdeS/IdeZ digestion**, which splits the "
                "antibody into F(ab')2 and Fc fragments before RP-LC/MS. Needs four files (naked "
                "+ ADC, for each fragment)."
            )
            st.markdown(
                "**Use this when:** you need to localize *where* the payload attached (Fab vs "
                "Fc), or confirm there's no unexpected off-target Fc conjugation."
            )
            if rp_page_obj is not None:
                st.page_link(rp_page_obj, label="Open RP LC-MS (IdeZ digestion) →", icon="\U0001F52C", use_container_width=True)
        if IDEZ_CARTOON.exists():
            st.image(
                str(IDEZ_CARTOON),
                caption="IdeZ digestion splits the antibody/ADC into F(ab')2 (~100 kDa) and Fc (~50 kDa) fragments, analyzed independently.",
                use_container_width=True,
            )

    st.divider()
    st.caption(
        "Not sure which one applies? If your deconvolution export is a single intact-mass file "
        "per sample, use SEC LC-MS. If your method note mentions IdeS, IdeZ, or F(ab')2/Fc "
        "fragments, use RP LC-MS."
    )

    st.markdown("#### New here?")
    st.caption(
        "Don't have a deconvolution file handy? Run the SEC LC-MS workflow on a small fabricated "
        "example dataset to see the whole app - metrics, charts, tables, downloads - before "
        "uploading anything of your own."
    )
    if st.button("▶ Try it with example data (SEC LC-MS)", key="home_try_example"):
        mab_buf, adc_buf = build_example_files()
        result = run_fragment_analysis(
            mab_buf, adc_buf, EXAMPLE_PAYLOAD_DEFS, EXAMPLE_ABUNDANCE_THRESHOLD,
            EXAMPLE_PPM_TOLERANCE, EXAMPLE_BASE_MASS_MODE, EXAMPLE_BASE_MASS_TOP_N,
        )
        if result is not None:
            st.session_state["sec_results"] = result
            st.session_state["_sec_results_are_example"] = True
            if sec_page_obj is not None:
                st.switch_page(sec_page_obj)
