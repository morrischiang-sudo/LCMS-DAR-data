"""
DAR Compass - ADC DAR Analysis Platform
-----------------------------------------
Thin router. Actual page content lives in home_page.py (chooser/landing
page), sec_page.py (SEC LC-MS Analysis), and rp_page.py (RP LC-MS / IdeZ
digestion Analysis); shared UI helpers and the password gate live in
dar_ui_helpers.py; the pure matching/DAR engine lives in dar_calculator.py.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import streamlit as st

import home_page
import rp_page
import sec_page
from dar_ui_helpers import APP_TITLE, require_password

st.set_page_config(page_title=APP_TITLE, page_icon="\U0001F9ED", layout="wide")

require_password()

sec_st_page = st.Page(sec_page.render, title="SEC LC-MS Analysis", icon="\U0001F9EA", url_path="sec")
rp_st_page = st.Page(rp_page.render, title="RP LC-MS (IdeZ digestion)", icon="\U0001F52C", url_path="rp-idez")
home_st_page = st.Page(home_page.render, title="Home", icon="\U0001F3E0", url_path="home", default=True)

# home_page.render() needs the same Page objects registered below so
# st.page_link can jump to a callable-based page (st.page_link requires the
# exact Page object for callable pages - a bare path string won't work here).
st.session_state["_pages"] = {"sec": sec_st_page, "rp": rp_st_page}

pg = st.navigation([home_st_page, sec_st_page, rp_st_page])
pg.run()
