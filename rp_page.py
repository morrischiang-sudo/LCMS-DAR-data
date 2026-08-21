"""
RP LC-MS (IdeZ digestion) Analysis page.

IdeS (from Streptococcus pyogenes) and IdeZ (its counterpart, optimized for
mouse IgG) are cysteine proteases that cleave IgG at a single, highly
specific site below the hinge, producing an F(ab')2 fragment (~100 kDa) and
Fc/2 fragment(s) per antibody (Sjoegren et al. 2016, The Analyst). This is
now a standard technique for ADC characterization specifically because it
localizes where a payload sits: the smaller ~100 kDa F(ab')2 fragment gives
better MS sensitivity/resolution than the ~150 kDa intact ADC, which matters
most for site-specific conjugates where the payload sits in the Fab region
(Su et al. 2016, Analytical Chemistry). A companion workflow independently
analyzes the Fc fragment to confirm there's no unexpected off-target
conjugation there and to check glycosylation (Rouviere et al. 2013, mAbs).

This page runs the exact same matching/DAR engine as the SEC LC-MS page,
independently for each fragment: naked Ab F(ab')2 is the baseline for ADC
F(ab')2, and naked Ab Fc is the baseline for ADC Fc. Per an explicit design
decision, F(ab')2 and Fc are reported as two fully independent results with
no combined "whole antibody" Total DAR - combining them correctly would
require knowing whether the uploaded Fc masses represent a full Fc dimer or
a single Fc/2 subunit (one heavy chain's worth), which varies by sample
prep and isn't something this page assumes.
"""

from __future__ import annotations

import streamlit as st

from dar_calculator import MassVariant, PayloadDef
from dar_ui_helpers import (
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
    st.sidebar.caption("RP LC-MS (IdeZ digestion) Analysis")

    st.sidebar.header("1. Upload deconvolution files")
    st.sidebar.markdown("**Naked antibody (reference)**")
    mab_fab2_file = st.sidebar.file_uploader(
        "Naked Ab — F(ab')2 (.xlsx)", type=["xlsx"], key="rp_mab_fab2",
        help="Deconvoluted intact-mass export for the unconjugated antibody's F(ab')2 fragment after IdeS/IdeZ digestion.",
    )
    mab_fc_file = st.sidebar.file_uploader(
        "Naked Ab — Fc (.xlsx)", type=["xlsx"], key="rp_mab_fc",
        help="Deconvoluted intact-mass export for the unconjugated antibody's Fc fragment after IdeS/IdeZ digestion.",
    )
    st.sidebar.markdown("**ADC (sample)**")
    adc_fab2_file = st.sidebar.file_uploader(
        "ADC — F(ab')2 (.xlsx)", type=["xlsx"], key="rp_adc_fab2",
        help="Deconvoluted intact-mass export for the ADC's F(ab')2 fragment after IdeS/IdeZ digestion.",
    )
    adc_fc_file = st.sidebar.file_uploader(
        "ADC — Fc (.xlsx)", type=["xlsx"], key="rp_adc_fc",
        help="Deconvoluted intact-mass export for the ADC's Fc fragment after IdeS/IdeZ digestion.",
    )

    st.sidebar.header("2. Linker-payload chemistries")
    st.sidebar.caption(
        "MW and mass variants are shared across both fragments (it's the same physical payload) - "
        "the max conjugation count is set separately per fragment, since where the payload ends up "
        "is exactly what this analysis investigates."
    )
    n_payloads = st.sidebar.number_input(
        "Number of distinct linker-payload chemistries",
        min_value=1, max_value=4, value=2, step=1, key="rp_npayloads",
    )

    payload_defs_fab2: list[PayloadDef] = []
    payload_defs_fc: list[PayloadDef] = []
    for i in range(int(n_payloads)):
        with st.sidebar.expander(f"Chemistry {i + 1}", expanded=(i < 2)):
            include = st.checkbox("Include in analysis", value=True, key=f"rp_include_{i}")
            label, variants = render_chemistry_variants_ui("rp", i)

            st.markdown("**Conjugation range — F(ab')2**")
            max_n_fab2, step_fab2 = render_conjugation_range_ui("rp_fab2", i, fragment_label="F(ab')2")

            st.markdown("**Conjugation range — Fc**")
            max_n_fc, step_fc = render_conjugation_range_ui("rp_fc", i, fragment_label="Fc")

            if include:
                lbl = label.strip() or str(i + 1)
                payload_defs_fab2.append(PayloadDef(
                    label=lbl,
                    variants=[MassVariant(mw=v.mw, dar_weight=v.dar_weight) for v in variants],
                    n_values=list(range(0, max_n_fab2 + 1, step_fab2)),
                ))
                payload_defs_fc.append(PayloadDef(
                    label=lbl,
                    variants=[MassVariant(mw=v.mw, dar_weight=v.dar_weight) for v in variants],
                    n_values=list(range(0, max_n_fc + 1, step_fc)),
                ))

    st.sidebar.header("3. ADC candidate selection")
    adc_min_fractional_abundance = st.sidebar.number_input(
        "ADC fractional abundance threshold (%)",
        min_value=0.0, max_value=100.0, value=0.0, step=0.01, format="%.2f",
        key="rp_abundance_threshold",
        help="Applied to both the F(ab')2 and Fc ADC files. Leave at 0 to consider every peak in each file.",
    )

    st.sidebar.header("4. Matching parameters")
    ppm_tolerance = st.sidebar.slider(
        "Mass accuracy tolerance (ppm)",
        min_value=10, max_value=500, value=300, step=10, key="rp_ppm",
        help="Applied to both fragments. Validated against historical SEC LC-MS runs at ~250-350 ppm; tune per instrument/method.",
    )
    base_mass_mode = st.sidebar.radio(
        "Base mass(es) to use",
        options=["Most abundant only (recommended)", "Top N glycoform/adduct variants"],
        key="rp_base_mass_mode",
    )
    base_mass_top_n = 1
    if base_mass_mode.startswith("Top N"):
        base_mass_top_n = st.sidebar.number_input("N", min_value=1, max_value=10, value=3, step=1, key="rp_base_mass_topn")

    run_clicked = st.sidebar.button("Run RP LC-MS analysis", type="primary", use_container_width=True, key="rp_run")

    # ----------------------------------------------------------------------
    # Main panel
    # ----------------------------------------------------------------------
    st.title(f"\U0001F9EA {APP_TITLE}")
    st.caption("RP LC-MS (IdeZ digestion) Analysis")
    st.info(
        "IdeS/IdeZ digestion cleaves the antibody below the hinge into an F(ab')2 fragment and Fc "
        "fragment, analyzed here independently using the exact same matching/DAR engine as the SEC "
        "LC-MS page. Naked Ab F(ab')2 is the baseline for ADC F(ab')2; naked Ab Fc is the baseline "
        "for ADC Fc. There is no combined \"whole antibody\" total shown here - F(ab')2 and Fc are "
        "reported as fully independent results.",
        icon="\U0001F9EA",
    )

    if run_clicked:
        result_fab2 = run_fragment_analysis(
            mab_fab2_file, adc_fab2_file, payload_defs_fab2, adc_min_fractional_abundance,
            ppm_tolerance, base_mass_mode, base_mass_top_n,
        )
        result_fc = run_fragment_analysis(
            mab_fc_file, adc_fc_file, payload_defs_fc, adc_min_fractional_abundance,
            ppm_tolerance, base_mass_mode, base_mass_top_n,
        )
        if result_fab2 is not None:
            st.session_state["rp_results_fab2"] = result_fab2
        if result_fc is not None:
            st.session_state["rp_results_fc"] = result_fc

    has_fab2 = "rp_results_fab2" in st.session_state
    has_fc = "rp_results_fc" in st.session_state

    if not has_fab2 and not has_fc:
        st.info(
            "Upload all four files in the sidebar (naked Ab and ADC, for both F(ab')2 and Fc), confirm "
            "your linker-payload chemistries and tolerance, then click **Run RP LC-MS analysis**.\n\n"
            "Expected file format: a deconvoluted intact-mass export (e.g. Thermo BioPharma Finder) "
            "with at least `Average Mass` and `Sum Intensity` columns."
        )
        st.stop()

    st.header("F(ab')2 results")
    if has_fab2:
        render_results_section(st.session_state["rp_results_fab2"], key_prefix="rp_fab2", title_suffix=" (F(ab')2)")
    else:
        st.info("No F(ab')2 results yet — upload both F(ab')2 files and click Run.")

    st.divider()

    st.header("Fc results")
    if has_fc:
        render_results_section(st.session_state["rp_results_fc"], key_prefix="rp_fc", title_suffix=" (Fc)")
    else:
        st.info("No Fc results yet — upload both Fc files and click Run.")
