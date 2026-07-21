# Foam Wash Column Model — Handoff

## Where this lives
- **Repo:** `mtlair/bigblue_mixing`
- **Branch:** `claude/foam-wash-column-review-pgfauf`
- **Model file:** `foam_wash_column_psd.R` (repo root)
- **Rendered output:** `output/foam_wash_column_psd.pdf` (3-panel plot) — ⚠️ re-render after the
  gas-state fix; the committed PDF still shows the old "gas-loss" panel 3. `Rscript foam_wash_column_psd.R`
- **Related (separate) model:** `morris_sensitivity_analysis.R` — upstream atomization/drying sensitivity study; not part of this column.

## Run it
```r
# R 4.3.3 with deSolve (1.40). Debian/Ubuntu: apt-get install r-base-core r-cran-desolve
Rscript foam_wash_column_psd.R          # prints diagnostics + writes plots to a device
```

## What it models
A **steady-state, 1-D (height z) preformed-foam wash column**. Particle-loaded foam is made
UPSTREAM (no gas sparger in this unit) and enters the base as a wet dispersion. The column
retains-vs-loses what arrives and washes entrained impurity out. Marched bottom→top as an
ODE system in `z` with `deSolve::ode`.

### Two zones (split at fixed decant depth `H_pool`)
- **Bottom — gravity-decant pool (`z < H_pool`):** liquid separates by density-differential
  **Stokes settling with hindered-settling (Richardson–Zaki) correction**. Pool velocity is
  computed, not fitted. `H_pool ≈ 0.45 m` (1–2 ft, stable interface).
- **Top — plug-flow foam (`z ≥ H_pool`):** throughput plug flow at the (very low) feed
  velocity `U_up`; washed by liquid added at the top.

### State variables
`Js_fine, Js_mid, Js_crs` (solid flux per size class), `C_imp` (impurity %), `eps_l` (liquid
holdup), `d_b` (mean bubble size), `J_g` (gas flux), `t_res` (cumulative residence).

### Mechanisms implemented
- **Residence:** throughput-based, ~1.5–2 h total (pool ~0.15 h + foam ~1.6 h). Interface is a
  fixed depth (observed stable in plug flow); feed rate sets plug residence.
- **Decant (impurity):** Stokes `U_settle = (ρ_liquid−ρ_foam)·g·d_sep²/(18μ) · (1−φ)^4.65`.
  Completes only if settling keeps up with upflow: `completion = min(1, U_settle/U_up)` — slow
  feed → clean separation. **This reproduces the plant decanter heuristic `t[h]=100μ/(ρA−ρB)`
  to within ~15%**, which de-magics that constant.
- **Bubble population:** `d_b(z)` grows by coalescence (faster in dry/unstable foam, ∝
  1/film_stability) minus breakage (restores toward inlet size). Feeds channeling + settling cut.
- **Gas state (conserved):** no sparger and no vent. Coalescence and bursting still happen in
  the column, but the gas is **not released here** — bursting past `d_b_burst` converts
  dispersed foam gas (`J_g_foam`) into large **retained slug gas** (`J_g_slug`) that can only
  escape in the next (downstream) unit up. Total gas `J_g_foam + J_g_slug` is conserved and
  carried out the top. **Bursting still coarsens bubbles and is the single mechanism for coarse
  dumping**, but it no longer removes gas from the column.
- **PSD loss:** buoyancy-failure detachment (`k_det`, coarse worst) + burst-driven collapse.
  Fines buoyant → retained; coarse stripped.
- **Regime map:** pool / snow-globe / plug. Above jamming, `plug_index = film_stability /
  (load_sens·loading)` — stable/lightly-loaded → plug; weak films/heavy load → snow-globe flakes.
- **Solids content:** `eps_s/(eps_s+eps_l)` — low top solids (~3–7%) is the **wet drainage
  equilibrium**, not incomplete drainage.
- **Pressure:** sets the inlet bubble size upstream (passed in as `d_b_in`).

## Calibrated baseline (current params)
| Quantity | Value |
|---|---|
| Residence | 1.74 h (0.16 pool + 1.58 foam) |
| Top solids content | 6.3% |
| Bubble size | 2.0 → 3.7 mm |
| Gas out top | 100% conserved (top split ~foam 70% / retained slug 30%) |
| Plug-flow height | 91% |
| Retention | fine 92% / mid 64% / coarse 9% |
| Impurity | 100% → 48% |

### Operational cliff (the key finding)
Film stability controls coalescence, which is make-or-break. As `film_stability` drops
1.5 → 0.5: `d_b` runs away 3 → 7 mm, the gas converts almost entirely to slug (foam-gas
fraction collapses ~100% → 2%), product and the plug regime disintegrate. Foam stability is the
variable the column lives or dies by. (Total gas is still conserved — it just leaves the top as
slug instead of dispersed foam.)

## Inputs that should come from the UPSTREAM unit (not fitted)
`mu_cont` (viscosity), `rho_foam`, `rho_liquid`, `rho_p`, particle sizes `d_fine/d_mid/d_crs`,
inlet bubble size `d_b_in`, decant-controlling droplet size `d_sep`, inlet holdups
`eps_l_in/eps_g_in/eps_s_in`, inlet loading `Js_*_in`.

## RESOLVED
- **✅ Gas is not vented in this column (fixed).** Gas is no longer a vent sink. `J_g` is split
  into `J_g_foam` (dispersed) + `J_g_slug` (retained), total conserved and carried out the top;
  bursting/coalescence convert foam gas → slug gas, which is only released in the next
  (downstream) unit up. Bursting still coarsens bubbles and dumps coarse particles. Implemented
  in `foam_wash_column_psd.R` section (3).

## OPEN ITEMS (priority order)
1. **Calibrate the bubble/gas kinetics** — `K_coal, K_break, K_burst, d_b_burst` are the
   least-anchored constants and now drive the dominant behavior. Highest-value data:
   measured bubble-size profile (inlet vs overflow) and actual gas split to the next stage.
2. **Pin `d_sep`** (~45 µm placeholder) — sets the decant timescale. Should be the real fine
   dispersed/emulsion droplet size from upstream.
3. **Couple gas state → holdup:** bursting (foam → slug) should make the collapsing foam wetter
   near the top (currently the gas split and `eps_l` evolve semi-independently).
4. **Bubble growth → drainage feedback:** bigger bubbles → wider Plateau borders → faster
   drainage (not yet coupled).
5. Optional fidelity: full population balance (bimodal distributions, fine-bubble tail) instead
   of a lumped mean diameter; counter-current wash liquid adding to holdup.

## Modeling notes / assumptions
- All closures are PLACEHOLDER forms calibrated to two plant anchors: **total residence
  1.5–2.5 h** and **top solids 3–7%**. Directional trends are trustworthy; absolute numbers
  need the data above.
- Verified end-to-end in R 4.3.3 / deSolve 1.40; single clean solve (zone switch is a smooth
  blend, no discontinuity).
- History of the derivation (why each mechanism was added) is in the git log on this branch.
