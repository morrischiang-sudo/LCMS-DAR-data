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

## Deploying on Render (public URL, e.g. to share with a colleague outside the network)

Render can host this as a normal Python web service. Feasibility-wise there's nothing
unusual here — pandas/numpy/matplotlib/openpyxl are all plain wheels, and the app never
writes anything to disk (uploads and results are handled entirely in memory for that
session), so Render's ephemeral filesystem isn't a problem.

**Before deploying:** Render gives every web service a public `onrender.com` URL with
no built-in login — anyone with the link can open it. Since this handles proprietary
ADC data, the app includes an optional password gate (see below) that's strongly
recommended once it's on a public URL, even if you only ever share the link with one
colleague.

1. Push this folder to the GitHub (or GitLab/Bitbucket) repo you're using for the project.
2. In the [Render dashboard](https://dashboard.render.com), choose **New > Blueprint**
   and point it at that repo — Render will read `render.yaml` in this folder and
   pre-fill the service config (Python runtime, build/start commands, free plan).
   - Prefer clicking through manually instead? Choose **New > Web Service**, connect
     the repo, and set: **Build Command** `pip install -r requirements.txt`,
     **Start Command** `streamlit run app.py --server.port $PORT --server.address 0.0.0.0`.
3. When prompted for the `APP_PASSWORD` environment variable, enter a password to
   share with your colleague (Render stores it as a secret, not in the repo).
4. Deploy. You'll get a `https://dar-compass-xxxx.onrender.com`-style URL — share that
   plus the password with your colleague, and nothing else public.

**Plan choice:** the Free instance (in `render.yaml`) spins down after 15 minutes of
no traffic and takes about a minute to spin back up on the next request — fine for
occasional use, but noticeable if your colleague opens it cold. Change `plan: free` to
`plan: starter` in `render.yaml` ($7/month) if that cold start becomes annoying; it
removes the spin-down entirely with the same 512 MB RAM.

If your organization has a policy on hosting internal R&D data on third-party cloud
services (even transiently, in memory), it's worth confirming this is fine before
deploying — that's a compliance question, not a technical one.

## How to use it

1. Upload the naked mAb reference deconvolution export (`.xlsx`) and the ADC
   deconvolution export (`.xlsx`) in the sidebar. Both need at least
   `Average Mass` and `Sum Intensity` columns (standard BioPharma-Finder-style
   export format).
2. Set up each linker-payload chemistry you expect: a short label, the maximum
   conjugation count, and whether counts step by 1 or 2 (use 2 if conjugation only
   happens in pairs, e.g. per interchain disulfide). For each chemistry, enter one
   **mass variant** (its MW in Da) if it's a normal intact linker-payload. If you've
   seen unexpected linker-payload breakage under harsh sample processing, increase
   "Number of mass variants" and add the additional possible mass(es) — DAR Compass
   models breakage as happening independently per attachment site, so a single
   molecule can have some sites intact and others broken at once. Each variant has
   its own **DAR weight** (1.0 = counts as a full payload like today; lower or 0.0
   if that broken form represents partial or total payload loss — this is a
   chemistry-specific call only you can make).
3. Optionally set the **ADC fractional abundance threshold (%)**. Only ADC peaks
   at or above this value are treated as candidate species before matching —
   useful for ignoring very low-abundance/noise peaks up front. Leave at 0 to
   consider every peak in the file.
4. Set the ppm mass-accuracy tolerance. 250-350 ppm reproduced the manually
   curated species assignments in both validation datasets; tune this per
   instrument/method rather than trusting a single fixed number.

   **Using more than one mass variant on a chemistry pushes this a lot harder than
   it looks.** Adding a second mass variant doesn't just add one more theoretical
   mass — it adds every possible split of that chemistry's site count between the
   two variants, so the theoretical mass ladder gets much denser. In testing, adding
   one breakage variant with a wide 350 ppm tolerance flagged the large majority of
   matches as ambiguous — not a bug, but a sign the tolerance and/or max conjugation
   count need tightening once multiple variants are in play. Start narrower than you
   would for a single-variant chemistry and widen only as far as you can while
   keeping ambiguous flags manageable. The "Mass variants configured for this run"
   panel and the grid-size warning (shown if a configuration would build an unusually
   large number of theoretical species) are there to help you judge this.
5. Click **Run DAR analysis**. Review any peaks flagged **ambiguous** (more
   than one plausible species within tolerance) before trusting the reported
   DAR — this is the automated version of the "does this look right?" check
   that used to be a manual judgment call.
6. Check the **Drug-load distribution** table and chart: one row/panel per
   chemistry, one column per count (0, 1, 2, ...), showing the % of matched
   intensity at each count independent of the other chemistries, plus an
   Average column (the same DAR number as the metric above). This is the
   standard drug-load-distribution report format — one row per linker-payload,
   blank cells where a count is outside that chemistry's configured range,
   the modal count highlighted. If a chemistry has a mass variant with a
   DAR weight other than 1.0, a note explains why the row's simple mean
   won't exactly equal the Average column.
7. If any chemistry has more than one mass variant configured, a **Chart detail**
   toggle also appears above the species-level distribution chart further down:
   **Detailed** shows every intact/broken mixture as its own bar (precise, but
   can get crowded); **Consolidated** groups bars by total occupied-site count
   per chemistry only, hiding which specific variant(s) made up that total.
   This only changes that chart — the table and DAR numbers are identical
   either way. Results stay on screen while you switch this (or any other
   widget) back and forth; only clicking **Run DAR analysis** again recomputes
   them.
8. Download the DAR report (Excel, now including the drug-load distribution
   as its own sheet) and both charts (PNG) for your records.

## Files in this delivery

- `app.py` — the Streamlit UI (includes an optional `APP_PASSWORD`-gated login, active only
  when that environment variable/secret is set — see the Render section above)
- `dar_calculator.py` — the underlying matching/DAR engine (also usable standalone/scripted)
- `requirements.txt` — pinned dependencies
- `render.yaml` — Render Blueprint for one-step deployment (see Render section above)
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
- Authentication is a single shared password (via `APP_PASSWORD`), not per-user login.
  Fine for sharing with one or two named colleagues; if this grows to a bigger audience
  or needs individual accounts/audit trail, revisit with a real auth solution.
