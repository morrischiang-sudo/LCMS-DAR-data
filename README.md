# DAR Compass — ADC DAR Analysis Platform

Automated Drug-to-Antibody Ratio (DAR) distribution analysis from LC-MS deconvoluted
intact-mass exports. See `DAR_platform_strategy_summary.md` for the background review
and validation this app is built on.

## Quick start (local, single user)

```bash
pip install -r requirements.txt
streamlit run app.py
```

This opens the app in your browser at `http://localhost:8501`.

## Quick start (shared internally, e.g. on a lab server)

```bash
pip install -r requirements.txt
streamlit run app.py --server.address 0.0.0.0 --server.port 8501
```

Anyone on the same network can then reach it at `http://<server-hostname>:8501`.
For anything beyond a same-network share (e.g. VPN-only access, login-gated access),
loop in IT/infrastructure to put it behind your usual reverse proxy / SSO setup.

## How to use it

1. Upload the naked mAb reference deconvolution export (`.xlsx`) and the ADC
   deconvolution export (`.xlsx`) in the sidebar. Both need at least
   `Average Mass` and `Sum Intensity` columns (standard BioPharma-Finder-style
   export format).
2. Set up each linker-payload chemistry you expect: a short label, its MW in Da,
   the maximum conjugation count, and whether counts step by 1 or 2 (use 2 if
   conjugation only happens in pairs, e.g. per interchain disulfide).
3. Optionally set the **ADC fractional abundance threshold (%)**. Only ADC peaks
   at or above this value are treated as candidate species before matching —
   useful for ignoring very low-abundance/noise peaks up front. Leave at 0 to
   consider every peak in the file.
4. Set the ppm mass-accuracy tolerance. 250-350 ppm reproduced the manually
   curated species assignments in both validation datasets; tune this per
   instrument/method rather than trusting a single fixed number.
5. Click **Run DAR analysis**. Review any peaks flagged **ambiguous** (more
   than one plausible species within tolerance) before trusting the reported
   DAR — this is the automated version of the "does this look right?" check
   that used to be a manual judgment call.
6. Download the DAR report (Excel) and chart (PNG) for your records.

## Files in this delivery

- `app.py` — the Streamlit UI
- `dar_calculator.py` — the underlying matching/DAR engine (also usable standalone/scripted)
- `requirements.txt` — pinned dependencies
- `DAR_platform_strategy_summary.md` — the workflow review, validation results, and
  platform design recommendations this app implements
- `results_A090412.xlsx`, `results_A090512.xlsx`, `results_A090412_from_raw_files.xlsx`,
  and their matching `dar_distribution_*.png` charts — the validation runs referenced
  in the strategy summary

## Known limitations / next steps

- Fragment/subunit-ion cross-validation (light chain, heavy chain, half-antibody) is
  not yet implemented — it needs user-supplied expected fragment masses. Planned as
  a phase-2 module (see strategy summary, section 5).
- No persistence/database yet — each run is stateless. If the team wants a history
  of past runs, that's a natural next feature (e.g. SQLite or a shared folder of
  saved reports).
- No authentication. Fine for a same-network internal tool; add if this grows beyond
  the immediate team or the data becomes more sensitive.
