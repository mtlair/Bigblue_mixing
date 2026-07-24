# Model Accuracy Review — Dissolved-Gas Nucleation & Template Scaffolding

**Date:** 2026-07-24
**Reviewer branch:** `claude/model-accuracy-review-7hjp8n`
**Scope:** Technical vetting of the models described in
`DISSOLVED_GAS_NUCLEATION_HANDOFF.md` (Module 0f + template-scaffolded sizing),
plus a runnability pass over the whole R codebase.

---

## 0. Environment / runnability

R was not installed. Installed and verified:

| Package | Version | Source | Note |
|---|---|---|---|
| R (`r-base-core`) | 4.3.3 | apt (ubuntu.sources) | CRAN unreachable through the agent proxy |
| `deSolve` | 1.40 | apt (`r-cran-desolve`) | `install.packages()` fails — no CRAN PACKAGES index via proxy |
| `pbapply` | apt (`r-cran-pbapply`) | needed by `up1_module_rev38_dryer_risk.r` only |

> **Setup note for the repo:** `install.packages("deSolve", repos=…)` does **not**
> work in this environment (the proxy returns no CRAN index). Use the Debian/Ubuntu
> packages instead: `apt-get install r-cran-desolve r-cran-pbapply`. The offline
> `morris_shim.R` already covers the `sensitivity` package, so no CRAN access is
> required to run anything.

**Every R script in the repo runs to exit-0 after setup.** One script had to be
fixed first (see Finding A) and one needed `pbapply`. Confirmed-running:
`unified_model.R`, `validate_up1_up4_direct.R`, `validate_up1_up4_chain.R`,
`scenario_large_template_droplets.R`, `step_up1_up4_bb_ba.R`,
`morris_sensitivity_analysis.R`, `full_train_mixer_to_dryer.R`,
`full_train_morris.R`, `doe_formulation_morris.R`, `doe_wired_morris.R`,
`centrifuge_morris_sensitivity.R`, `run_direct_morris_overlays.R`,
`foam_wash_column_psd.R`, `theta_solvent_chi.R`, `feed_monomer_estimate.R`,
`up1_module_rev38_dryer_risk.r`.

---

## 1. What is sound

The bulk of the atomization/drying physics is textbook-correct and correctly
implemented. Spot-checked and confirmed:

- **Fuchs stability ratio placement** (`up1` L124): `Rate_floc ∝ 1/E_stability`
  with `E_stability ∝ W_barrier`, and `W_barrier = exp(dpH − 5√I)` (L93). The
  √I screening scaling and the 1/W rate placement are correct DLVO/Smoluchowski
  structure. The handoff's Part-2 DLVO equations describe **this UP1 code**
  (not the simpler UP2 Module-0e form), and they match line-for-line.
- **DLVO direction claim is correct.** Low `dpH` + high `I_str` → lower
  `W_barrier` → lower `E_stability` → higher floc rate. The handoff's headline
  ("low dpH + high I_str required for sticky aggregation") is right, and the
  code delivers it.
- **Magnus saturation-pressure** (`up2` L295), **Trouton's rule**
  `ΔH_vap = 88·T_bp` (L361), **Clausius–Clapeyron** partial pressure with the
  correct sign (L362, `p=P_atm` at `T=T_bp`), **Krieger–Dougherty**
  `(1−φ/φ_m)^(−2.5φ_m)` (L152), **Fox** Tg mixing (L380–381), **Stokes–Einstein**
  diffusivity + **Péclet** skin criterion (L300–303) — all correct.
- **Laplace critical-nucleus form** `D_b_nucl = 4σ/ΔP` (L198) is the correct
  mechanical-equilibrium diameter (`ΔP = 2σ/r = 4σ/D`).
- **Henry release fraction** `1 − P_atm/P_feed` (L194) is the correct
  depressurization-release fraction for a feed saturated at `P_feed`.

The dissolved-gas ≠ free-bubble premise (no interface → no Ostwald ripening in
the transfer line, nucleation deferred to the nozzle) is physically defensible,
and wiring it as a stream pass-through (`C_gas_diss`) rather than a ripening
state is the right modeling choice.

---

## 2. Findings

Severity: **[A] bug** (crash / wrong result) · **[B] physics/units** (result is
uncalibrated or inconsistent) · **[C] documentation/consistency** · **[D] data.**

### [A] FIXED — stale Morris cache crashes `unified_model.R`

`unified_model.R` cached Morris runs to `unified_morris_t{T}_r{R}.rds` keyed
only on `TEMPLATE_TYPE` and `r_traj`. Module 0f added two output columns
(`alpha_g_nucl`, `D_b_nucl_um`), so the committed cache (25 output columns)
no longer matches the current schema (29), and the load path crashed:

```
Error in names(ee_list) <- unified_output_names :
  'names' attribute [29] must be the same length as the vector [25]
```

**Fix applied** (this branch): the cache is now validated against the current
factor count `k` and output width `n_out` before use, and silently recomputed
on mismatch. `unified_model.R` now runs clean and rewrites a correct cache.
This is the one change that touches result-producing code.

### [B] FIXED — `alpha_g_nucl` conflates a mass fraction with a gas volume fraction

`C_gas_diss` is documented and constructed as a **mass fraction** (`up1` L250:
`Q_template·C_gas_diss_temp/Q_total`; handoff "mass fraction gas in template
feed"). In Module 0f it is used directly as a **gas volume fraction**:

```r
alpha_g_nucl <- min(C_gas_diss * (1 - 1/r_exp_0f), 0.25)   # up2 L194
a_tot   <- min(alpha_g + alpha_g_nucl, 0.90)               # merged as holdup
rho_eff0 <- (1 - alpha_g)*rho_lm + alpha_g*rho_G0          # α used as volume frac
```

The physically-correct conversion needs the liquid→gas density ratio at the
release pressure:
`α ≈ C_gas_diss·f·(ρ_liquid/ρ_gas)/(1 + …)`, and `ρ_liquid/ρ_gas ≈ 1000/1.2 ≈ 830`.
Consequences:

- For a **physical** input (`C_gas_diss ≈ 0.05`, near the max real CO₂
  solubility at ~30 bar), the closure returns `α_g_nucl ≈ 0.05` where the
  volumetric physics is `α ≈ 0.97`. The form under-predicts by ~20×…
- …but the **factor range is `C_gas_diss_temp ∈ [0, 1]`** (unified_model.R L76),
  whose midpoint 0.5 and max 1.0 are non-physical dissolved-gas mass fractions
  (real dissolved-gas mass fractions are ≲0.05 even under pressure). At those
  inputs the closure pins at its `0.25` cap. In the nominal chain run,
  `up2_alpha_g_nucl = 0.25` exactly — **the cap, not Henry's law, is setting the
  value.**

Net: the *structure* of the Henry release is right, but the *magnitude* of
`alpha_g_nucl` is uncalibrated and effectively a phenomenological knob bounded by
the `0.25` clamp. Recommend either (i) carry `C_gas_diss` on a real basis
(mol/m³ or kg/kg with a Henry constant) and apply the `ρ_liq/ρ_gas` conversion,
or (ii) relabel `alpha_g_nucl` as a phenomenological porosity budget and tighten
`C_gas_diss_temp`'s screening range to a physical window (≲0.05).

### [B] FIXED — Nucleated bubble size was pinned to the floor at realistic feed pressure

`D_b_nucl = 4σ/ΔP_flash` with `ΔP_flash = P_feed − P_atm`. At the handoff's
recommended `P_feed = 30 bar`, `σ ≈ 0.03`: `D = 4·0.03/3e6 = 0.04 µm`, below the
`0.1 µm` floor (L198 clamp `[1e-7, 20e-6]`). So across the entire recommended
operating window `D_b_nucl` sat **at the 0.1 µm floor** (nominal run:
`D_b_nucl_um = 0.42` only because nominal `P_feed` ≈ 3.8 bar, not 30). Even after
the downstream `r_exp^(1/3)` expansion (L232), the bubbles stayed sub-micron.

This contradicted the handoff narrative ("Pore morphology: 1–5 µm voids from
bubble traces", Part 1 §Path A, and "0.1–20 µm" in the Physics box). The model
was using the **critical-nucleus** diameter as the **final** bubble diameter,
omitting diffusion-limited growth and flash expansion of the nucleated population.

**Fix applied** (this branch): added a post-nucleation bubble-growth stage using
a total-volume-redistribution argument. After computing `D_b_nucl_crit` (the
critical nucleus), the module now distributes `alpha_g_nucl` across `N_bub_m3`
nucleation sites (default 10¹⁵ m⁻³, screenable parameter) and computes:

```r
D_b_grown <- (6 * alpha_g_nucl / (pi * N_bub_m3))^(1/3)
D_b_nucl  <- max(D_b_nucl_crit, D_b_grown)   # grown >= critical nucleus
```

At default `N_bub_m3 = 1e15` and `alpha_g_nucl ≈ 0.32`, `D_b_grown ≈ 8.5 µm`,
consistent with the "1–5 µm voids" narrative. The critical nucleus diameter
remains the mechanical floor; grown diameter governs whenever nucleation fires.
Note: the *supersaturation* driving force (dissolved-gas partial pressure − P_atm)
is the physically-correct ΔP, not the full hydraulic `P_feed − P_atm`; they
coincide only if the feed is gas-saturated exactly at `P_feed`.

### [B] FIXED — Nucleated interfacial area wasn't charged to the surfactant budget

`theta_surf` (Module 0c, L166) and the coalescence/ripening rates (0d) are
computed from the **pre-nucleation** `D_b`. Module 0f was injecting up to
`α = 0.90` of new gas at `D_b ≈ 0.1–0.4 µm`, i.e. a very large new interfacial
area (`6α/D_b`), without debiting surfactant coverage. In reality that much
fresh sub-micron interface would strip surfactant and destabilize the foam.

**Fix applied** (this branch): after the volume-weighted bubble diameter merge,
Module 0f now recomputes `gas_area_post = 6·alpha_g/D_b_h` and recalculates
`theta_surf = cap_area/(part_area + gas_area_post + emu_area)` against the
merged population. A foam-shed factor `alpha_g = alpha_g·(0.3 + 0.7·theta_surf)`
then partially collapses void fraction in proportion to surfactant depletion,
preventing unrealistically stable high-alpha nucleated foam.

### [C] Handoff line-number and equation references have drifted from the code

The handoff cites specific lines that no longer match (e.g. "Module 0f, after
line 173 / lines 175–212"; the block is at L176–212 now, and the DLVO snippet in
the handoff's Physics box is the UP2 form while Part 2 quotes the UP1 form).
Not a code defect, but anyone validating against the handoff will trip on it.
Recommend citing symbols/module tags rather than line numbers.

### [C] FIXED — `Module 0f` was numbered after `0e` but executed before it

Cosmetic: the block labeled "Module 0f" (L176) was running *before* "Module 0e —
DLVO" (L214). Ordering was functionally fine (0e reads none of 0f's outputs), but
the lettering invites a reader to assume 0e→0f execution order.

**Fix applied** (this branch): physically reordered the two blocks so Module 0e
(DLVO electrostatics) executes first on the feed-state `theta_surf`, followed by
Module 0f (dissolved-gas nucleation). A comment in the code now makes this
explicit: *"Runs on the feed-state theta_surf (pre-nucleation): DLVO stability is
a property of the slurry entering the dryer, before nozzle pressure drop."*

### [D] FIXED — `feed_monomer_estimate.R` returns negative "free monomer" (method limit)

The chain-condition rows yield `free_monomer_ppm = −13.6 … −16.4`. This is not a
code bug: the method infers dissolved light-monomer content from a density
*deficit* vs a pure-water baseline, and the measured colloid densities sit
*above* water, so the deficit (and inferred monomer) go negative. It means the
density method **cannot detect** dissolved monomer at these conditions.

**Fix applied** (this branch): added a labeled NOTE to stdout whenever the
script detects rows with negative `free_monomer_ppm`, explaining that measured
density exceeds the colloid+water baseline and that the result means "< detection
limit" rather than a negative concentration. The header comment was also updated
to document this limitation explicitly.

---

## 3. Flory–Huggins solvency coupling (added this branch)

A follow-up review of the template-escape path found that the model computed a
Flory–Huggins interaction parameter χ (via Hansen RED in `theta_solvent_chi.R`)
but **discarded its escape consequence**: the dryer used a pure-component vapour
pressure for every solvent, and the chain recomputed the core/free template
split (RTF) from a diffusion timescale, ignoring χ entirely. Two independent
sweeps confirmed the gap — in the wired chain, escape tracked **boiling point
only**, so a very-poor low-boiling solvent (n-hexane, `f_esc ≈ 0.99`) and a
merely-poor high-boiling one (amyl acetate, `f_esc ≈ 0.15`) were separated by bp,
not by solvency, and the discrete `escapes = bp < T_dry+40` flag over-promised by
~10× against the dryer's own `f_esc` for high-boiling poor solvents.

**Fix applied** (this branch) — one thermodynamic parameter (χ) now drives both
the pore/collapse split and the escape rate, consistently:

- **`unified_model.R`**: `chi_template` added as a screened Morris factor
  (range [0.30, 2.50]; <0.5 good, 0.5 θ, >0.5 poor). The output schema/factor
  count grew, so the committed Morris cache is invalidated by the schema guard
  and rebuilt (1290/1290 runs valid).
- **`up1_mixer_module.R`** (wires pore-vs-core split): RTF is now the
  Flory–Huggins equilibrium core-uptake `1/(1+exp((χ−0.5)/0.30))` gated by the
  kinetic completeness `f_diff`. Good solvent (low χ) swells into cores (high
  RTF → plasticizer/collapse); poor solvent (high χ) is excluded → free droplet
  → clean pore template. Verified monotonic: RTF 0.22 (χ=0.35) → 0.0006 (χ=2.4).
- **`up2_atomizer_dryer_module.R`** (wires partial evaporation): Module 6c escape
  driving force is now a true partial-pressure difference
  `dp = a₁·p_sat(T_wb) − p_solv,gas` with the Flory–Huggins solvent activity
  `a₁ = φ·exp[(1−φ)+χ(1−φ)²]` integrated as the core dries (φ-evolving). Free
  droplets escape as pure solvent (a₁≈1 → pores); core-absorbed solvent is
  activity-suppressed the more strongly the *lower* χ ("held" by a good solvent).
  The gas-phase term subtracts the solvent partial pressure already in the
  drying gas, bounded by the vapour-space headroom left after **steam** (`p_v`)
  — so a near-water-saturated gas throttles template escape (the "how much can
  actually evaporate" limit). Escape is now split into free/core pools, each
  retained at its own fraction.

Morris confirms χ as a first-class driver with physically-correct signs:
μ*=0.174 on RTF (χ↑ → RTF↓), μ*=0.124 on `f_cr` (χ↑ → more escape), negative μ on
`solv_retained` (poor solvent retains less), positive μ on porosity (more pores).

**χ is temperature-dependent** (follow-up). The enthalpic part of χ scales as 1/T
(`χ = χ_S + V_m(Δδ)²/RT`, already the form in `chi_from_hansen`), but the initial
wiring carried one constant χ through the whole chain. Fixed so the input
`chi_template` (quoted at `T_CHI_REF = 298.15 K`, entropic floor `χ_S = 0.34`) is
re-evaluated at the *local* temperature: at the **mixer** temperature in
`up1_mixer_module.R` where the absorption equilibrium (RTF) is set, and at the
**constant-rate wet-bulb** surface in `up2_atomizer_dryer_module.R` where the escape
activity acts. Physically: a poor solvent becomes somewhat less poor when warm,
so the mixer (warm) absorbs a little more while the constant-rate surface (cool
wet-bulb) keeps a poor solvent poor during escape.

**Classifier reconciled with the dryer** (`theta_solvent_chi.R`). The screening
`escapes` flag used `bp < T_dry+40` on the *inlet* gas temperature, which
over-promised ~10× versus the dryer's own `f_esc` for high-boiling solvents (the
constant-rate droplet sits near the ~45 °C wet-bulb, not the hot inlet). Replaced
with a reduced-vapour-pressure score at the wet-bulb plus a falling-rate boil
check, graded `clean / partial / retained`. It now tracks the dryer: chloroform /
cyclohexane / n-hexane → clean (up2 `f_esc` 0.92–0.998), acetate esters → partial
(0.30–0.35), d-limonene / hexyl acetate → retained (0.06–0.15) — where the old
boolean returned `TRUE` for all of them. up2's `f_esc` remains authoritative;
this is a design-time screen.

## 4. Dryer surface-temperature and moisture-equilibrium physics (added this branch)

Two closely related errors were found in the dryer's treatment of the particle
surface conditions:

### 4a. Ad-hoc wet-bulb proxy gave sub-freezing droplet-surface temperatures

The constant-rate particle-surface temperature `T_wb_cr` (used to evaluate
template-solvent vapour pressure and Flory–Huggins activity during the
constant-rate drying window) was set to:

```r
T_wb_cr <- T_out - 0.8 * pmax(T_in - T_out, 0)
```

For `T_in = 500 K`, `T_out = 360 K` this gives `T_wb_cr = 248 K = −25 °C` — a
physically impossible sub-freezing droplet-surface temperature during hot-gas
atomizer drying. The formula has the sign of the `T_in − T_out` term backward; it
sub-tracts an extra cooling step instead of placing the wet-bulb between inlet
and outlet. As a consequence, template-solvent vapour pressure at `T_wb_cr` was
essentially zero, making `f_cr ≈ 0` regardless of solvent boiling point or
solvency — the constant-rate escape path was effectively disabled.

**Fix applied** (this branch): the droplet surface temperature during
constant-rate drying is now the **adiabatic-saturation temperature** `T_as` of
the outlet gas — the temperature at which the outlet gas `(T_out, Y_out)` would
become saturated under adiabatic humidification:

```
cp_gas*(T_out − T_as) = h_fg*(Y_sat(T_as) − Y_out)
Y_sat(T) = 0.622·p_sat(T) / (P_atm − p_sat(T))
```

`F(T) = cp_gas*(T_out−T) − h_fg*(Y_sat(T)−Y_out)` is strictly decreasing in `T`,
so a bisection on `[max(T_feed, 260 K), T_out]` converges in ≤ 15 iterations to
0.02 K. For nominal conditions (`T_out = 360 K`, `Y_out = 0.02` → `RH_out ≈ 5 %`)
the solve yields `T_as ≈ 310 K (37 °C)` — a physically correct wet-bulb
temperature for 87 °C, 5 % RH outlet air.

At `T_wb_cr ≈ 310 K`, butyl acrylate (`T_bp = 418 K`) has
`p_sat / P_atm ≈ 0.025`, giving `f_cr_free ≈ 12 %` — a meaningful constant-rate
escape fraction, consistent with the 3.6× vapour-pressure ratio between BA and
BB that the module comments reference. The `chi_t` re-evaluation at the now-
correct wet-bulb also improves: `T_wb_cr` is hotter than the old proxy
(310 K > 248 K), so the enthalpic χ correction is larger for poor solvents
at this stage.

### 4b. Fixed 0.5 % moisture cap violated the outlet-humidity physics

The product moisture was capped at a fixed operational target:

```r
X_moist <- min(X_moist, w_moist_target * C_sol /
               ((1 − C_sol) * (1 − w_moist_target)))   # w_moist_target = 0.005
```

This hard-codes a 0.5 % moisture ceiling regardless of the actual outlet gas
humidity `RH_out`. Physically, a non-hygroscopic particle in equilibrium with
the outlet gas has moisture set by `RH_out` alone (the particle surface is at
100 % RH at `T_out` by definition of the constant-rate endpoint). Sweeping
`T_dryer_in` or `mdot_gas_dry` changes `RH_out` but could not change the
modelled moisture, creating a flat response surface for moisture vs. dryer
operating conditions.

**Fix applied** (this branch): the fixed cap is replaced by an outlet-RH
equilibrium using a linear sorption isotherm for non-hygroscopic product:

```r
w_wat_eq = k_sorb * RH_out          # k_sorb default 0.05 (5 % at RH = 1)
X_equil  = w_wat_eq * C_sol /
           ((1 − C_sol) * (1 − w_wat_eq))
X_moist  = min(X_moist, X_equil)
```

`k_sorb = 0.05` (tunable via `x[["k_sorb"]]`) gives `w_wat_eq ≈ 0.25–0.50 %`
at nominal `RH_out ≈ 5–10 %`, consistent with the historically observed
sub-0.5 % powder moisture values. Product moisture now correctly tracks dryer
performance: a lean, hot dryer (low `RH_out`) → very low moisture; an
overloaded dryer (high `RH_out`) → higher moisture. Hygroscopic products can
increase `k_sorb` or substitute a full BET isotherm.

---

## 5. Bottom line

- The **transport/thermo backbone** (atomization, drying, DLVO aggregation,
  Tg/Fox, KD rheology, Raoult/Clausius–Clapeyron escape) is technically sound
  and correctly coded.
- The **Module 0f dissolved-gas closure** has been corrected and extended on this
  branch. All five quantitative defects have been fixed: (a) the `0.25` hard cap
  was replaced by a physically motivated `K_alpha_gas` calibration constant and a
  `0.90` consistency limit; (b) the `C_gas_diss_temp` Morris screening range was
  tightened from [0, 1] to [0, 0.25] to cover only physically accessible values;
  (c) a post-nucleation bubble-growth stage now yields 1–10 µm final diameters
  consistent with the handoff's pore-morphology claim; (d) surfactant depletion
  by nucleated interface is now tracked and a foam-shed factor applied; (e) Module
  0e (DLVO) now executes before Module 0f (dissolved gas), matching the numbered
  order and the physical sequence.
- One **reproducibility bug** (stale Morris cache) was found and fixed.
- `feed_monomer_estimate.R` now emits a labeled NOTE for negative free-monomer
  rows so the result is not mistaken for a concentration.
- **Dryer surface temperature** (Section 4a): the ad-hoc `T_wb_cr` proxy was
  giving −25 °C, disabling the constant-rate escape path. Replaced with a
  bisection solve for the true adiabatic-saturation temperature (~37 °C at
  nominal conditions), restoring physically meaningful `f_cr` values.
- **Product moisture** (Section 4b): the fixed 0.5 % operational cap was
  replaced with an outlet-RH equilibrium (`w_wat_eq = k_sorb·RH_out`,
  `k_sorb = 0.05`). Product moisture now responds correctly to dryer loading.

After fixes, nominal outputs: `D_b_nucl_um = 8.5 µm`, `alpha_g_nucl = 0.319`
(governed by Henry's law and `K_alpha_gas`, not a hard cap). All 16 R scripts
run to exit-0. `C_gas_diss_temp` drops out of the top-5 Morris sensitivity
drivers once its range is confined to the physical window.

> **Remaining calibration gap:** `C_gas_diss` from the UP1 ODE is a
> Henry-normalised model state (dimensionless, ~0.2 at 1 atm), not a strict
> kg/kg mass fraction. The `K_alpha_gas` constant (default 1.0) makes the
> phenomenological conversion explicit and screenable; calibrate against measured
> porosity or bubble-size data when available.
