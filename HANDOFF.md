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
  1/film_stability, surfactant-set) minus breakage (restores toward inlet size) minus a
  **burst sink** (`K_bsink`) that removes the coarsest bubbles — so coarsening self-limits and
  `d_b` stays bounded even in the cliff (it no longer runs to physically absurd sizes).
- **Burst trigger = film rupture.** Bursting fires when the film thins to rupture, expressed as a
  **critical bubble size that is itself set by film physics**:
  `d_b_crit = d_b_burst · film_stability^a_fs · (σ_ref/σ)^a_sig` — weak films (low Gibbs
  elasticity) or high surface tension → smaller `d_b_crit` → burst sooner. The rate is then
  modulated by drainage speed and particle effects (see table below). At baseline `d_b_crit ≈
  3.0 mm` and `d_b` sits just above it, so the column runs **right at its burst threshold** —
  which is why it is sensitive to the operating point.
- **Gas state (conserved):** no sparger and no vent. Coalescence and bursting still happen in
  the column, but the gas is **not released here** — bursting converts dispersed foam gas
  (`J_g_foam`) into large **retained slug gas** (`J_g_slug`) that can only escape in the next
  (downstream) unit up. Total gas `J_g_foam + J_g_slug` is conserved and carried out the top.
  **Bursting still coarsens/removes bubbles and is the single mechanism for coarse dumping**,
  but it no longer removes gas from the column.

#### What affects the burst (verified by sweep)
| Variable | Effect on burst | Why | Channel |
|---|---|---|---|
| Surfactant dose (elasticity, σ) | **Strong** — starving → runaway burst (cliff) | weak films / high σ → small `d_b_crit` | `film_stability`, `σ` |
| Temperature | **Strong, competing** — hot → *less* burst here | σ↓ & elasticity↑ stabilize (win) vs μ↓ destabilizes | `σ(T)`, `E∝RT`, `μ(T)` |
| Viscosity | **Yes** — higher μ → less burst | slower film drainage (Reynolds ∝ 1/μ) | `M_visc=(μ_ref/μ)^n_visc` |
| Local solids **content** | **Yes** — more solids → less burst | particles armor the film (+ raise μ_eff) | `M_armor=1/(1+k_armor·φ_cond)` |
| Local solids **size** (coarse) | **Weak** — coarser → slightly more burst | large particles bridge/rupture thin films | `M_bridge` |
| Pressure | **None in-column** | only sets ρ_gas and (upstream) `d_b_in`; the big expansion effect is **downstream** | — |
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
- **Film thickness (wired into drainage, holdup, AND burst):** Plateau-border radius
  `r_pb = c_pb·d_b_eff·√eps_l`, capillary suction `P_c = σ/r_pb`, DLVO equilibrium film thickness
  `h_eq = λ_D·ln(Π_charge·θ / P_c)` (~48 nm baseline). Couplings:
  - **Drainage rate** ← border permeability: `tau_local = tau_drain·(r_pb_ref/r_pb)^p_perm`
    (p_perm=2). Bigger bubbles → fatter borders → **faster drainage** (old item #4).
  - **Equilibrium holdup** ← disjoining/capillary balance: a thicker equilibrium film (stronger
    disjoining Π_charge·θ, lower σ, longer λ_D) holds a **wetter** foam,
    `eps_l_dry ← eps_l_dry·film_stability^hold_exp·film_wet^q_film`, `film_wet = h_eq/h_eq_ref`.
  - **Burst** ← film rupture: `d_b_crit` shrinks as films thin (see burst table above).
  - **`h_eq(z)` is now robust:** `d_b_eff = min(d_b, d_b_cap)` in the border geometry and the
    burst `d_b`-sink keep `d_b` bounded, so `h_eq` no longer spuriously rises when bubbles run
    away — in the cliff it correctly *drops* (~48→40 nm) as the foam dries.
  All factors normalized to 1.0 at baseline. Verified: starving surfactant or shrinking λ_D
  thins the equilibrium film → drier foam; doubling coalescence fattens borders → `tau` drops.

## Calibrated baseline (current params)
| Quantity | Value |
|---|---|
| Residence | 1.74 h (0.16 pool + 1.58 foam) |
| Top solids content | 6.3% |
| Bubble size | 2.0 → 3.5 mm (burst-limited; `d_b_crit` ≈ 3.0 mm) |
| Gas out top | 100% conserved (top split ~foam 70% / retained slug 30%) |
| Plug-flow height | 91% |
| Retention | fine 92% / mid 64% / coarse 9% |
| Impurity | 100% → 47% |
| Thermo state | T 298 K, P 1.5 bar → μ 2.0 mPa·s, ρ_gas 1.75 kg/m³, σ 21 mN/m |
| Surfactant | c 5 mol/m³ (CMC 6) → coverage 0.91, Gibbs E 99 mN/m → film_stability 1.00 |
| Film thickness | eq. `h_eq` ~42–48 nm (border r_pb ~283 µm; drainage τ ~1800 s) |

### Operational cliff (the key finding)
Film stability controls coalescence, which is make-or-break — and it is now **set by the
surfactant** (coverage → Gibbs elasticity), not a free knob. Verified in the model:
- **Starve surfactant** to `c=1 mol/m³` (below CMC): elasticity collapses → `film_stability`
  0.20, `d_b_crit` shrinks → runaway coalescence, `d_b` grows (now **burst-limited to ~7.5 mm**,
  not an absurd 90 mm), gas goes ~100% to slug, top solids fall to ~3%, plug regime disintegrates.
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
- **✅ Film thickness wired into drainage + holdup (was item #4).** Plateau-border radius
  `r_pb(z)` sets the drainage timescale (permeability ∝ r_pb²: bigger bubbles drain faster);
  the disjoining/capillary equilibrium film thickness `h_eq` sets the equilibrium holdup
  (thicker films → wetter foam). `h_eq(z)`/`r_pb(z)`/`tau_local(z)` emitted as profiles.
- **✅ Film-rupture burst trigger + `h_eq(z)` fixed.** Bursting fires off a film-rupture critical
  size `d_b_crit(film_stability, σ)`, modulated by viscosity / solids content / coarse-particle
  size (see burst table). A burst `d_b`-sink caps runaway coarsening, so `d_b` stays physical and
  `h_eq(z)` no longer spuriously rises in the cliff. Burst dependence on T/P/μ/solids swept and
  documented. Bubble panel now plots `d_b` vs `d_b_crit`.

## OPEN ITEMS (priority order)
1. **Calibrate the film/kinetics closures** — `K_coal, K_break, K_burst, K_bsink, d_b_burst` and
   the surfactant/film constants (`E_stab_ref, Gamma_inf, K_ads, cmc, mu_surf, Π_charge, λ_D,
   c_pb, a_fs, a_sig, n_visc, k_armor, k_bridge`, the `r_pb_ref`/`h_eq_ref` normalizers) are
   placeholder forms. **The net temperature sign of the burst is set by `a_fs/a_sig/dsigma_dT`
   and is not yet anchored** — pin it with data. Highest-value data: measured bubble-size profile
   (inlet vs overflow), surface tension vs dose (σ–c isotherm), and the actual gas split.
2. **Pin `d_sep`** (~45 µm placeholder) — sets the decant timescale. Should be the real fine
   dispersed/emulsion droplet size from upstream.
3. **Couple gas state → holdup:** bursting (foam → slug) should make the collapsing foam wetter
   near the top (currently the gas split and `eps_l` evolve semi-independently).
4. Optional fidelity: full population balance (bimodal distributions, fine-bubble tail) instead
   of a lumped mean diameter; counter-current wash liquid adding to holdup.

## Modeling notes / assumptions
- All closures are PLACEHOLDER forms calibrated to two plant anchors: **total residence
  1.5–2.5 h** and **top solids 3–7%**. Directional trends are trustworthy; absolute numbers
  need the data above.
- Verified end-to-end in R 4.3.3 / deSolve 1.40; single clean solve (zone switch is a smooth
  blend, no discontinuity).
- History of the derivation (why each mechanism was added) is in the git log on this branch.
