"""
Home / chooser page.

First thing a new user sees: two large, descriptive cards explaining when to
use each workflow, with a direct link into each. `app.py` builds the actual
st.Page objects (st.page_link needs the same Page object that was registered
with st.navigation for callable-based pages) and hands them to this module
via st.session_state before calling st.navigation.
"""

from __future__ import annotations

import streamlit as st

from dar_ui_helpers import APP_TAGLINE, APP_TITLE


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

    st.divider()
    st.caption(
        "Not sure which one applies? If your deconvolution export is a single intact-mass file "
        "per sample, use SEC LC-MS. If your method note mentions IdeS, IdeZ, or F(ab')2/Fc "
        "fragments, use RP LC-MS."
    )
