# DoE proposal — wired chain UP1→UP2→UP3→UP4

**Date:** 2026-07-23 · Source: `doe_wired_morris.R` → `unified_output/doe_wired_morris_indices.csv`,
`doe_proposal.txt`. Morris elementary-effects screen (41 factors, r=30, 1260 runs)
on the **fully wired** chain (foam wash + centrifuge enabled, size-template
morphology path on), so the dryer sees the concentrated ~40 % cake.

## Top drivers per product CQA (model ranking)

| CQA | Top controllable drivers | Direction | SEM anchor |
|-----|--------------------------|-----------|------------|
| **d50** particle size | v_tip, C_solid_mass, ionic_strength, (n_flow, C_surfactant, ΔpH) | ↑v_tip→↑ρ but ↓d50; ↑solids→↓d50 | — (laser PSD) |
| **porosity** | C_solid_mass ↑, template_dose ↑, T_dryer_in ↑, Q_template ↓ | more solids/template/heat → more porosity | **~0.37, robust** |
| **sphericity** | T_dryer_in ↓, C_monomer ↑, T_bp_solv ↑, σ ↑ | hotter inlet → less round | **~0.91, robust** |
| **tapped density** | v_tip ↑, ionic_strength ↓, C_solid_mass ↓, ΔpH ↑ | — | — |
| **product Tg** | Tg_polymer ↑ (dominant), C_plasticizer ↓, C_monomer ↓ | plasticizer/monomer soften | — |
| **residual moisture** | Q_template ↓, C_solid_mass ↑, T_dryer_in ↓ | — | — |
| **skin/fusion** | mu_L ↑, T_mix ↓, T_dryer_in ↑, v_tip ↓ | more viscous/hotter → more skin | **DEFERRED** (not repeatable) |

## Recommended DoE — split by controllability

The Morris "top-4 of ≥2 CQAs" set mixes three kinds of factor. A process DoE
should move only what the operator sets; formulation and intrinsic-material
factors belong to separate studies.

### A. Process setpoints — the primary DoE (run these)
Bounds narrowed to the **real operating window** from the validation data
(`cond_process.csv` / `cond_up1234.csv`), not the wide screening range.

| factor | unit | DoE low | DoE high | nominal | CQAs it drives |
|--------|------|---------|----------|---------|----------------|
| **v_tip** (UP1 tip speed) | m/s | 6 | 20 | 16 | d50, tapped density, skin |
| **C_solid_mass** (feed solids) | wt/wt | 0.15 | 0.30 | 0.25 | d50, porosity, tapped density, moisture |
| **T_dryer_in** (UP4 inlet) | K | 410 | 435 | 425 | porosity, sphericity, moisture |
| **ALR** (atomizing air-liquid ratio) | – | **0.1** | 1.8 | 0.6 | span/breadth, droplet size |

→ A **2⁴ full factorial (16 runs) + 3 center points**. Responses: d50 & span
(laser), **porosity (SEM, robust)**, **sphericity (SEM, robust)**, tapped density,
residual moisture. Resolves all main effects + 2-factor interactions (v_tip and
C_solid_mass both show Morris σ ≳ μ*, i.e. real interactions); add a face-centred
axial set → CCD if porosity/d50 curvature needs a quadratic surface.

**Changes from the raw Morris set:**
- **ALR replaces P_atom_air**, and is the right choice because it is dimensionless
  and **atomizer- and scale-agnostic** (air/liquid mass ratio) — the effect
  transfers across nozzles and throughputs, unlike an absolute pressure. P_atom_air
  is not independently settable and co-varies with ALR anyway.
- **ALR range extended DOWN to 0.1** (from the measured 0.9–1.8) to probe the
  low-atomization region where droplet size and span change fastest. 0.05 can be
  added as an extra low axial point. **Caveat:** the atomization closure was
  calibrated at ALR 0.9–1.8, so model predictions below ~0.9 are extrapolation —
  which is exactly why it is worth an experimental point rather than a model claim.
  ALR mainly moves the distribution **span** (model rank #3) — median d50 is
  dominated by v_tip/solids once the feed is a ~40 % paste.
- **Q_template dropped** — no data and not a current operational variable. Its CQA
  roles (porosity, moisture) are still carried by C_solid_mass and T_dryer_in.

**Why this set is a solid basis (not just model ranking):** every factor is (a) a
real operator setpoint, (b) top-ranked for ≥1 CQA in the wired screen, and (c)
bounded by *measured* operating data — v_tip 4.4–19.3 m/s and C_solid_mass 13–25 %
(`cond_process.csv`), T_dryer_in ≈143 °C and ALR 0.9–1.8 (`cond_up1234.csv`). The
chain that ranks them is validated to product d50 RMS 1.10. So the effects are
anchored to hardware you can set and to data, not to free model parameters.

### B. Formulation / surface-chemistry DoE — the higher-value study
Focused Morris (r=40, wired chain, process held at nominal) over the
monomer / plasticizer / binder / solvent design space
(`doe_formulation_morris.R` → `doe_formulation_proposal.txt`). This is where the
**morphology and structure** CQAs live, and one factor dominates:

**C_plasticizer is the master formulation lever — top-3 for 8 of 10 CQAs.**
It sets the softness/free-volume axis and pulls almost everything:

| CQA | dominant formulation driver(s) | direction |
|-----|-------------------------------|-----------|
| product Tg | **C_plasticizer**, C_monomer | both soften (↓Tg) |
| porosity | template_dose (solvent amount) ↑, **C_plasticizer** ↓ | more solvent → more pores |
| sphericity | **C_plasticizer** ↑ | plasticized → rounder |
| shell permeability | k_perm_plast ↑, C_binder ↓, **C_plasticizer** ↑ | plasticizer opens the shell |
| pore size | **C_plasticizer** ↓, C_monomer ↓ | |
| retained solvent | **T_bp_solv** ↑ | higher-bp solvent stays in |
| micro-explosion risk | **T_bp_solv** ↓, C_binder ↑ | high-bp solvent → less burst |
| cake strength | **C_plasticizer** ↓, C_monomer ↓, C_binder ↓ | softeners weaken the cake |
| skin/fusion | **C_plasticizer**, C_monomer, C_binder ↑ | all promote fusion |

**Recommended formulation DoE — dosable knobs (proposed ranges):**

| factor | role | DoE low | DoE high | nominal | drives |
|--------|------|---------|----------|---------|--------|
| **C_plasticizer** | plasticizer conc. | 0.005 | 0.04 | 0.0225 | Tg, sphericity, porosity, pore, strength, permeability, skin (8/10) |
| **C_monomer** | monomer conc. | 0.001 | 0.03 | 0.0155 | Tg, pore size, strength, skin (4/10) |
| **C_binder** | binder conc. | 0.01 | 0.15 | 0.08 | permeability, burst, strength, skin (4/10) |
| **T_bp_solv** | solvent choice (bp) | 300 K | 360 K | 330 | solvent retention, burst, d50 (4/10) |
| **template_dose** | solvent amount | 0.02 | 1.1 | 0.56 | porosity, d50 (3/10) |

→ Two clean axes: **plasticizer+monomer** = softness/Tg/strength; **binder** =
permeability/burst; **solvent amount + type** = porosity/retention/burst. A
**2⁵⁻¹ resolution-V (16 runs) + centers** resolves these with C_plasticizer as the
expected dominant main effect. Hold the process box (section A) at center.
`k_perm_mono/plast/bind` rank high but are each additive's intrinsic shell-
permeation property — they tell you *which additive chemistry* the response keys
on (plasticizer permeation dominates), not a knob to set; screen candidate
chemistries against them rather than dialing them.

### C. Intrinsic/material — characterize, don't set
`n_flow`, `mu_L` (slurry rheology), `k_perm_plast`, `T_bp_solv`, `σ`, `ionic_strength`,
`ΔpH`. These rank high (n_flow is the top d50 driver, mu_L the top skin driver) but
are **material/serum properties**, not direct setpoints — they move through solids,
formulation and additive chemistry. Measure them per run (viscometry, ζ-potential)
as covariates so the process-DoE effects are not confounded by drifting rheology.

## Calibration status of the responses
- **Porosity (~0.37) and sphericity (~0.91): robust SEM anchors** — usable as DoE
  responses and to calibrate `phi_porosity_z` / `Omega_struct_z` (note the
  definitional gap: SEM surface-porosity vs model bulk `phi_porosity_z`).
- **Skin/fusion: deferred.** SEM skin is not repeatable per sample yet
  (view-to-view scatter > between-sample spread). Before using it as a DoE
  response, collect **several high-mag views per sample at a fixed magnification**
  and average. Until then, `theta_skin_z` stays model-only.

## Caveats
- Rankings are **model-based**; the wired chain matches product d50 to RMS 1.10 but
  UP3 viscosity and cake solids are calibrated (not first-principles), so treat the
  DoE as hypothesis-generating — confirm effects experimentally.
- The wired chain requires `size_template = 1`; the muted path overpredicts at
  40 % cake solids.
- Bounds in section A are the realistic operating window; widen only within safe
  process limits.
