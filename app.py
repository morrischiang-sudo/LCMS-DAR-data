"""
DAR Compass - ADC DAR Analysis Platform
-----------------------------------------
Thin router. Actual page content lives in sec_page.py (SEC LC-MS Analysis)
and rp_page.py (RP LC-MS / IdeZ digestion Analysis); shared UI helpers and
the password gate live in dar_ui_helpers.py; the pure matching/DAR engine
lives in dar_calculator.py.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import streamlit as st

import rp_page
import sec_page
from dar_ui_helpers import APP_TITLE, require_password

st.set_page_config(page_title=APP_TITLE, page_icon="\U0001F9ED", layout="wide")

require_password()

pg = st.navigation([
    st.Page(sec_page.render, title="SEC LC-MS Analysis", icon="\U0001F9EA", url_path="sec", default=True),
    st.Page(rp_page.render, title="RP LC-MS (IdeZ digestion)", icon="\U0001F52C", url_path="rp-idez"),
])
pg.run()
