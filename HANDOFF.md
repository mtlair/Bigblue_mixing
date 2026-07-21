# Foam Wash Column Model — Handoff

## Where this lives
- **Repo:** `mtlair/bigblue_mixing`
- **Branch:** `claude/column-coalescence-bubble-burst-inrt41` (gas-state fix + surfactant/T/P
  chain; forked from `claude/foam-wash-column-review-pgfauf`)
- **Model file:** `foam_wash_column_psd.R` (repo root)
- **Rendered output:** `output/foam_wash_column_psd.pdf` (3-panel plot; panel 3 is the conserved
  gas-state split). Regenerate with `Rscript -e 'pdf("output/foam_wash_column_psd.pdf",width=12,height=4); source("foam_wash_column_psd.R"); dev.off()'`.
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
holdup), `d_b` (mean bubble size), `J_g_foam` + `J_g_slug` (gas state), `t_res` (residence).

### Mechanisms implemented
- **Residence:** throughput-based, ~1.5–2 h total (pool ~0.15 h + foam ~1.6 h). Interface is a
  fixed depth (observed stable in plug flow); feed rate sets plug residence.
- **Thermodynamic state (T, P):** `mu_cont(T)` (Andrade `μ = μ_ref·exp(B(1/T−1/T_ref))`),
  `σ(T)` linear, gas density `ρ_gas = P·MW/(RT)` (ideal gas). Pressure sets `ρ_gas` here and
  (upstream) `d_b_in`. All derived once in `derive_state_props()`.
- **Surfactant → film elasticity (steps 1–2):** Langmuir adsorption of `c_surf` → surface
  excess `Γ`; Szyszkowski `σ(Γ,T)`; Gibbs–Marangoni elasticity `E_gibbs = RT·Γ_inf·θ/(1−θ)`,
  rolled off above the CMC by micelle buffering. `film_stability := E_gibbs/E_stab_ref` — so
  surfactant type/dose now *drives* coalescence and the regime map instead of a hand-set knob.
- **Decant (impurity):** Stokes `U_settle = (ρ_liquid−ρ_foam)·g·d_sep²/(18·μ_eff) · (1−φ)^4.65`.
  Completes only if settling keeps up with upflow: `completion = min(1, U_settle/U_up)` — slow
  feed → clean separation. **This reproduces the plant decanter heuristic `t[h]=100μ/(ρA−ρB)`
  to within ~15%**, which de-magics that constant.
- **Krieger–Dougherty crowding:** as liquid drains, solids concentrate in the Plateau borders
  (`φ_cond = eps_s/(eps_s+eps_l)`), raising the local suspension viscosity
  `μ_eff = μ(T)·(1−φ_cond/φ_smax)^(−2.5·φ_smax)`, which further hinders settling on top of the
  Richardson–Zaki holdup term. High local border viscosity is a real retention mechanism.
- **Bubble population:** `d_b(z)` grows by coalescence (faster in dry/unstable foam, ∝
  1/film_stability, now surfactant-set) minus breakage (restores toward inlet size).
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
- **Film drainage (physical, step 3):** `eps_l` relaxes toward an equilibrium holdup set by
  film stability (`eps_l_dry_eff = eps_l_dry·film_stability^hold_exp`) with a timescale
  `tau_drain_eff = tau_drain·(μ(T)/μ_ref)·mobility` (viscous, slowed by rigid/high-surface-
  viscosity interfaces). The old hand-set `tau_drain`/`eps_l_dry` are now baseline anchors
  scaled by dimensionless factors that are 1.0 at the baseline state.
- **Solids content:** `eps_s/(eps_s+eps_l)` — low top solids (~3–7%) is the **wet drainage
  equilibrium**, not incomplete drainage.
- **Film thickness (reporting):** Plateau-border radius proxy `r_pb ~ d_b·√eps_l`, capillary
  suction `σ/r_pb`, and a DLVO equilibrium film thickness `h_eq = λ_D·ln(Π_charge·θ / P_cap)`
  (~49 nm baseline). Illustrative — not fed back into the ODE (see open item #4).

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
| Thermo state | T 298 K, P 1.5 bar → μ 2.0 mPa·s, ρ_gas 1.75 kg/m³, σ 21 mN/m |
| Surfactant | c 5 mol/m³ (CMC 6) → coverage 0.91, Gibbs E 99 mN/m → film_stability 1.00 |
| Film thickness | eq. `h_eq` ~49 nm (r_pb ~300 µm) |

### Operational cliff (the key finding)
Film stability controls coalescence, which is make-or-break — and it is now **set by the
surfactant** (coverage → Gibbs elasticity), not a free knob. Verified in the model:
- **Starve surfactant** to `c=1 mol/m³` (below CMC): elasticity collapses → `film_stability`
  0.20 → `d_b` runs away to ~90 mm, gas goes ~100% to slug, top solids fall to ~0.3%, plug
  regime disintegrates (91% → 40%).
- **Overdose** to `c=12 mol/m³` (well above CMC): micelle buffering drops the *effective*
  elasticity → `film_stability` 0.60 → worse than baseline. Foam stability peaks near the CMC.
- **Heat to 333 K:** μ 2.0→1.1 mPa·s, σ 21→13 mN/m, faster drainage; `film_stability` rises to
  1.12 (E ∝ RT). **Raise pressure to 3 bar:** ρ_gas doubles (1.75→3.5 kg/m³).
(Total gas is always conserved — the cliff shifts it from dispersed foam gas to retained slug.)

## Inputs that should come from the UPSTREAM unit (not fitted)
`rho_foam`, `rho_liquid`, `rho_p`, particle sizes `d_fine/d_mid/d_crs`, inlet bubble size
`d_b_in`, decant-controlling droplet size `d_sep`, inlet holdups `eps_l_in/eps_g_in/eps_s_in`,
inlet loading `Js_*_in`, and the **operating point** `T_col, P_col` and **surfactant**
`c_surf` (with its isotherm `Gamma_inf/K_ads/cmc`). `mu_cont` is now derived from `T_col`.

## RESOLVED
- **✅ Gas is not vented in this column (fixed).** Gas is no longer a vent sink. `J_g` is split
  into `J_g_foam` (dispersed) + `J_g_slug` (retained), total conserved and carried out the top;
  bursting/coalescence convert foam gas → slug gas, which is only released in the next
  (downstream) unit up. Bursting still coarsens bubbles and dumps coarse particles. Implemented
  in `foam_wash_column_psd.R` section (3).
- **✅ Surfactant → film-elasticity → drainage chain (steps 1–3).** The scalar `film_stability`
  is now derived from surfactant coverage via Langmuir/Szyszkowski/Gibbs–Marangoni; drainage
  runs on a physical `tau_drain(μ,T,mobility)` toward a stability-set equilibrium holdup.
  Temperature and pressure added (`mu(T)`, `σ(T)`, ideal-gas `ρ_gas(P,T)`). Krieger–Dougherty
  border crowding added to settling. All in `derive_state_props()`; baseline unchanged.

## OPEN ITEMS (priority order)
1. **Calibrate the film/kinetics closures** — `K_coal, K_break, K_burst, d_b_burst` and the new
   surfactant constants (`E_stab_ref, Gamma_inf, K_ads, cmc, mu_surf`, disjoining `Π_charge`)
   are placeholder forms; the isotherm and `E_stab_ref` normalization should be pinned to the
   actual surfactant. Highest-value data: measured bubble-size profile (inlet vs overflow),
   surface tension vs dose (σ–c isotherm), and the actual gas split to the next stage.
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
