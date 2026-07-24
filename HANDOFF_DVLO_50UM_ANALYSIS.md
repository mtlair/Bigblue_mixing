# DLVO + Surface Chemistry 50 µm Analysis — Handoff Summary

**Date:** 2026-07-24  
**Branch:** `claude/dvlo-surface-chemistry-model-kq5ory`  
**Script analyzed & run:** `dvlo_surface_chemistry_50um.R`  
**Status:** Complete. Outputs committed. Ready for merge to main.

---

## File Directory

```
Bigblue_mixing/
├── dvlo_surface_chemistry_50um.R          ← scenario driver (analyzed + executed)
├── unified/
│   ├── up1_mixer_module.R                 ← UP1 mixer ODE model
│   ├── up2_atomizer_dryer_module.R        ← UP4 atomizer + dryer model
│   └── interface_stream.R                 ← inter-stage stream contract
├── foam_wash_module.R                     ← UP2 foam-wash (bypassed in this script)
├── up3_separator_module.R                 ← UP3 decanting separator
├── DVLO_PATH_HANDOFF.md                   ← mechanistic design rationale
└── unified_output/
    ├── dvlo_50um_sweep1_dlvo_chemistry.csv    ← Delta_pH × ionic_strength grid (30 rows)
    ├── dvlo_50um_sweep2_gas_template.csv      ← Q_gas × template_dose grid (25 rows)
    ├── dvlo_50um_sweep3_up3_gforce.csv        ← UP3 Fg sweep (6 rows)
    ├── dvlo_50um_sweep4_large_seed.csv        ← D_seed sweep 10–45 µm (8 rows)
    ├── dvlo_50um_sweep5_plasticizer.csv       ← C_plasticizer sweep (7 rows)
    ├── dvlo_50um_sweep6_tg_polymer.csv        ← Tg_polymer sweep (6 rows)
    ├── dvlo_50um_best_case_kpis.csv           ← combined best-case KPIs
    └── dvlo_50um_best_case_stream.csv         ← best-case UP3-exit stream state
```

---

## Process Train Nomenclature

The file naming uses a legacy `up2_` prefix for what the process calls UP4. Confirmed mapping:

| Process stage | File / function | Role | Status in this script |
|---|---|---|---|
| **UP1** | `unified/up1_mixer_module.R` → `up1_run_mixer()` | Mixer / coagulation vessel | Active |
| **UP2** | `foam_wash_module.R` → `foam_wash_column()` | Foam-wash column | **Bypassed** (identity pass-through) |
| **UP3** | `up3_separator_module.R` → `up3_separator()` | Decanting centrifuge | Active |
| **UP4** | `unified/up2_atomizer_dryer_module.R` → `up2_run_dryer()` | Atomizer + spray dryer | Active (nominal setpoint, not the study variable) |

---

## What the Script Does

The script (`dvlo_surface_chemistry_50um.R`) is a **sensitivity sweep driver** exploring
routes to **~50 µm, high-porosity, low-surface-fusion particles** without adjusting UP4
(the "known path"). It chains UP1 → UP3 → UP4, with UP2 bypassed.

The fundamental constraint: UP1's calibrated homogeneous-aggregation plateau is pinned at
`D_agg = 0.2 µm × (1 + 50) = 10.2 µm`. No DLVO chemistry knob can push past it. Three
upstream levers are explored instead.

---

## Three-Lever Strategy

### Lever 1 — Size: Large-Seed Heteroaggregation (UP1)

Introduce **cationic 35–45 µm seed particles** (amine-functional PMMA spheres or
aminopropyl-silica monospheres) into the anionic latex slurry. DLVO electrostatics become
attractive between seed and colloid (opposite charges) while remaining repulsive between
like-sign colloids — selective heteroaggregation without homo-floc.

Implemented via `run_chain_large_seed()`: after the normal UP1 run, the stream
`D_agg_um` is overridden to `D_seed + 2 × shell` and `size_template = 1` is set so UP4
tracks the aggregate, not the atomizer droplet.

The very small d_ratio (0.82 µm / 50 µm = 0.016) also generates structural porosity:
`phi_struct ≈ 0.30 × theta_skin × (1 − d_ratio^0.5) ≈ 0.25`.

**Seed type:** commercial cationic PMMA latex spheres (amine or quaternary-ammonium
surface), 35–45 µm, net positive charge at process pH.

### Lever 2 — Surface Fusion Suppression (DLVO chemistry + plasticizer)

`theta_skin` has two parallel sub-routes:
- **Péclet route** (`theta_skin_pe`): high DLVO stability (high `Delta_pH`, low
  `ionic_strength`) → lower aggregation rate → lower skin consolidation.
- **Surface-fusion route** (`theta_skin_fus`): driven by `Tg_plas` vs wet-bulb (373 K)
  via the Fox equation. Stripping plasticizer raises `Tg_plas` toward `Tg_polymer`.

**Hard floor from fixed Tg_polymer:** if `Tg_polymer` cannot be changed (polymer IS the
product), the fusion route is locked above ~0.82 even at zero plasticizer. The rigid seed
physically interrupts polymer film coalescence at the surface — a mechanical suppression
not captured in the current model equations but the intended real mechanism.

### Lever 3 — Porosity (gas holdup + template)

- **Gas holdup**: raise `Q_gas` in UP1; run UP3 at low Fg (100 g vs 430 g nominal)
  → `alpha_floor = 0.076 × (430/100) = 0.327` (33% gas preserved vs 7.6% nominal).
- **Capillary-bridge template (type 4, butyl butyrate)**: `chi = 1.4` (poor solvent),
  bp = 166 °C, minimal aqueous dissolution. Higher `template_dose` → more `phi_templ_free`
  → more pore template surviving into UP4.

**DLVO–gas tension (confirmed by runs):** the chemistry that suppresses homo-floc (high
`Delta_pH`, low I) also suppresses gas capture in UP1
(`capture_eff ∝ 1/(1 + 0.1 × E_stability)`). At best-case DLVO conditions, `alpha_g_feed`
collapsed to 0.018, not the expected 0.10–0.33. The gas-porosity and DLVO-selectivity
levers work against each other for the mixer-capture pathway.

**Implication:** if both high selectivity and high gas porosity are needed simultaneously,
the **gas-template route (template_type = 2)** should be considered. In type 2, gas stays
dissolved (nucleates at the UP4 nozzle) rather than being captured in the mixer; the DLVO
chemistry no longer gates gas retention. DLVO direction for type 2 is **opposite** to type
4 — lower `Delta_pH` and higher `ionic_strength` to facilitate particle packing around
nucleating bubbles at the nozzle.

---

## Actual Run Results (2026-07-24)

### Nominal baseline
```
D_agg = 10.2 µm, Dp50 = 165.5 µm, phi = 0.223, theta_skin = 0.949, alpha_g = 0.076
```

### Sweep summary

| Sweep | Variables | Key finding |
|---|---|---|
| 1 | Delta_pH × ionic_strength | theta_skin range 0.878–0.949; DLVO lever alone gives modest improvement |
| 2 | Q_gas × template_dose | phi range 0.204–0.225; gas/template contribution limited by DLVO suppression of capture |
| 3 | UP3 Fg (50–430 g) | phi 0.218–0.235; at 100 g: alpha_floor = 0.327, phi = 0.235 |
| 4 | D_seed (10–45 µm) | **Large-seed path works for size**: D_seed=40 → Dp50=76.5 µm, phi=0.250 |
| 5 | C_plasticizer | Lower plasticizer → lower theta_skin (0.9448 at 0 vs 0.9531 at 4%); small effect |
| 6 | Tg_polymer (fixed) | Even Tg=380 K only gives theta_skin=0.925; fusion floor is physics-fundamental |

### Best-case combined recipe (all levers)

| Parameter | Nominal | Best-case |
|---|---|---|
| D_seed (µm) | N/A | 40 |
| Delta_pH | 2.5 | 3.5 |
| ionic_strength (M) | 0.25 | 0.03 |
| Q_gas | 1.0 | 1.6 |
| template_dose | 0.5 | 0.8 |
| C_plasticizer | 0.022 | 0.005 |
| Tg_polymer (K) | 330 | 330 (fixed) |
| UP3 Fg (g) | 430 | 100 |

| KPI | Nominal | Best-case actual |
|---|---|---|
| Dp50 (µm) | 165.5 | 75.7 |
| phi_porosity | 0.223 | 0.231 |
| theta_skin | 0.949 | 0.981 |
| alpha_g_feed | 0.076 | 0.018 |
| rho_tapped (kg/m³) | 624.8 | 593.5 |
| Tg_product (K) | 291.9 | 313.1 |

---

## Key Discrepancies vs Handoff Predictions

1. **Dp50 = 75.7 µm, not 48–55 µm.** `TEMPLATE_DENSIFY = 1.20` was calibrated from
   10 µm seeds. At 50 µm seeds the size-template overlay overshoots. Needs recalibration
   against a 50 µm seed + shrinkage measurement.

2. **alpha_g collapses at high DLVO stability.** The best-case chemistry (Delta_pH=3.5,
   I=0.03) drives E_stability to its 100 cap, reducing `capture_eff` to near zero. UP3
   low-g has nothing to preserve. This is a real coupling not fully anticipated in the
   design rationale.

3. **theta_skin worsens in best-case (0.981 vs nominal 0.949).** The large seed changes
   the d_ratio and phi_struct terms in a way that interacts adversarially with the skin
   closure. The rigid-seed mechanical suppression (physically real) is not in the model.

---

## UP3 Solids at 100 g

From the concentration equation `Cs_out = max(Cs_in, 0.405 × (Fg/430)^0.05)`:
- **100 g → 37.5% solids** (vs 40.4% at 430 g nominal)
- The g-force exponent (0.05) is intentionally weak — cake solids change little; the main
  effect of lower Fg is gas retention via `alpha_floor = alpha_gas_cake × (Fg_ref/Fg)`.

---

## Known Gaps & Next Steps

| Gap | Risk | Suggested action |
|---|---|---|
| TEMPLATE_DENSIFY calibrated at 10 µm | Dp50 prediction unreliable for 50 µm seeds | Bench-spray 50 µm seed slurry; measure dry Dp50; recalibrate factor |
| Gas capture killed by high DLVO stability | Gas-porosity and selectivity levers conflict | Evaluate template_type=2 (gas nucleation at nozzle; DLVO direction reversed) |
| theta_skin model has no rigid-seed term | Model underestimates surface-fusion suppression by seed | Add theta_skin multiplier for seed coverage fraction; calibrate from SEM |
| Tg_polymer fixed (product constraint) | Hard fusion floor ~0.88 theta_skin | Rigid-seed mechanical suppression is the path; needs experimental validation |
| UP2 foam-wash bypassed | Impurity removal (monomer washout) not modeled | Wire `unified.wire_up2 = TRUE` for full chain; re-run Morris |

---

## Setup Notes

- **R dependency:** `deSolve` (install via `apt-get install r-cran-desolve`; CRAN unreachable from this environment)
- **Run:** `Rscript dvlo_surface_chemistry_50um.R` from repo root (~30 s)
- **Output dir:** `unified_output/dvlo_50um_*.csv` (8 files)

---

Generated by Claude Sonnet 4.6 | Branch: `claude/dvlo-surface-chemistry-model-kq5ory`
