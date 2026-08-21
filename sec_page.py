"""
SEC LC-MS Analysis page.

Intact-mass DAR analysis from a single naked-mAb reference file and a
single ADC file (the original DAR Compass workflow). See
DAR_platform_strategy_summary.md for the methodology this implements.
"""

from __future__ import annotations

import streamlit as st

from dar_calculator import PayloadDef
from dar_ui_helpers import (
    APP_TAGLINE,
    APP_TITLE,
    render_chemistry_variants_ui,
    render_conjugation_range_ui,
    render_results_section,
    run_fragment_analysis,
)


def render() -> None:
    # ----------------------------------------------------------------------
    # Sidebar - inputs
    # ----------------------------------------------------------------------
    st.sidebar.title(f"\U0001F9EA {APP_TITLE}")
    st.sidebar.caption("SEC LC-MS Analysis")

    st.sidebar.header("1. Upload deconvolution files")
    mab_file = st.sidebar.file_uploader(
        "Naked mAb reference export (.xlsx)",
        type=["xlsx"], key="sec_mab_file",
        help="Deconvoluted intact-mass export for the unconjugated antibody. Used to determine base masses (glycoform/adduct variants).",
    )
    adc_file = st.sidebar.file_uploader(
        "ADC export (.xlsx)",
        type=["xlsx"], key="sec_adc_file",
        help="Deconvoluted intact-mass export for the conjugated ADC sample.",
    )

    st.sidebar.header("2. Linker-payload chemistries")
    n_payloads = st.sidebar.number_input(
        "Number of distinct linker-payload chemistries",
        min_value=1, max_value=4, value=2, step=1, key="sec_npayloads",
    )

    payload_defs: list[PayloadDef] = []
    for i in range(int(n_payloads)):
        with st.sidebar.expander(f"Chemistry {i + 1}", expanded=(i < 2)):
            include = st.checkbox("Include in analysis", value=True, key=f"sec_include_{i}")
            label, variants = render_chemistry_variants_ui("sec", i)
            max_n, step = render_conjugation_range_ui("sec", i)
            if include:
                payload_defs.append(PayloadDef(
                    label=label.strip() or str(i + 1),
                    variants=variants,
                    n_values=list(range(0, max_n + 1, step)),
                ))

    st.sidebar.header("3. ADC candidate selection")
    adc_min_fractional_abundance = st.sidebar.number_input(
        "ADC fractional abundance threshold (%)",
        min_value=0.0, max_value=100.0, value=0.0, step=0.01, format="%.2f",
        key="sec_abundance_threshold",
        help=(
            "Only ADC peaks with Fractional Abundance at or above this value are treated as "
            "candidate species and passed on to matching/DAR calculation. Raise this to drop "
            "very low-abundance peaks (noise, minor adducts) before matching against the "
            "theoretical mass ladder. Leave at 0 to consider every peak in the file."
        ),
    )

    st.sidebar.header("4. Matching parameters")
    ppm_tolerance = st.sidebar.slider(
        "Mass accuracy tolerance (ppm)",
        min_value=10, max_value=500, value=300, step=10, key="sec_ppm",
        help="Peaks are matched to the nearest theoretical species and kept only if within this tolerance. Validated against historical runs at ~250-350 ppm for ~150-165 kDa intact species; tune per instrument/method.",
    )
    base_mass_mode = st.sidebar.radio(
        "mAb base mass(es) to use",
        options=["Most abundant only (recommended)", "Top N glycoform/adduct variants"],
        key="sec_base_mass_mode",
    )
    base_mass_top_n = 1
    if base_mass_mode.startswith("Top N"):
        base_mass_top_n = st.sidebar.number_input("N", min_value=1, max_value=10, value=3, step=1, key="sec_base_mass_topn")

    run_clicked = st.sidebar.button("Run DAR analysis", type="primary", use_container_width=True, key="sec_run")

    # ----------------------------------------------------------------------
    # Main panel
    # ----------------------------------------------------------------------
    # Results are computed only when "Run DAR analysis" is clicked, then
    # stashed in st.session_state and rendered from there on every
    # subsequent script run - otherwise results would disappear the moment
    # you touched any other widget (like the chart-detail toggle below),
    # since st.button() only reads True on the exact run it was clicked on.

    st.title(f"\U0001F9EA {APP_TITLE}")
    st.caption(f"SEC LC-MS Analysis — {APP_TAGLINE}")

    if run_clicked:
        result = run_fragment_analysis(
            mab_file, adc_file, payload_defs, adc_min_fractional_abundance,
            ppm_tolerance, base_mass_mode, base_mass_top_n,
        )
        if result is not None:
            st.session_state["sec_results"] = result

    if "sec_results" not in st.session_state:
        st.info(
            "Upload a mAb reference file and an ADC file in the sidebar, confirm your linker-payload "
            "chemistries and tolerance, then click **Run DAR analysis**.\n\n"
            "Expected file format: a deconvoluted intact-mass export (e.g. Thermo BioPharma Finder) "
            "with at least `Average Mass` and `Sum Intensity` columns."
        )
        st.stop()

    render_results_section(st.session_state["sec_results"], key_prefix="sec")
