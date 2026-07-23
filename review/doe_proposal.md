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
| **Q_template** (template/colloid ratio) | – | 0.1 | 2.0 | 1.05 | porosity, moisture, tapped density |
| **P_atom_air** (atomizing pressure) | psig | 3 | 11 | 7 | d50 (dominant atomizer knob) |

→ A **2^5 fractional factorial (resolution V, 16 runs) + 3 center points**, or a
face-centred CCD if curvature on porosity/d50 is expected. Responses: d50 & span
(laser), **porosity (SEM, robust)**, **sphericity (SEM, robust)**, tapped density,
residual moisture. This resolves the main effects and 2-factor interactions the
Morris σ/μ* flags as non-linear (v_tip, C_solid_mass, ionic_strength all show
σ ≳ μ*).

### B. Formulation factors — a separate formulation DoE
`C_monomer`, `C_plasticizer`, `C_binder`, `template_dose`, `Tg_polymer`.
These drive **Tg** (Tg_polymer dominant; plasticizer/monomer soften) and modulate
porosity/sphericity. Hold the process box at its center while varying these.

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
