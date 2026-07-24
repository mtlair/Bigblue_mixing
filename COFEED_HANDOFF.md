# Co-Feed Template-in-Water Mechanism: Summary Handoff

**Date:** 2026-07-23  
**Branch:** `claude/main-handoff-review-0qo61a`  
**Status:** Ready for implementation

---

## Executive Summary

A new **co-feed template-in-water mechanism** has been developed for atomizer-dryer pore formation. Instead of pre-emulsifying the template solvent into the polymer slurry, the template is injected as a fine aqueous dispersion directly into UP1 (the mixer) alongside the colloid feed. This approach:

- **Decouples template droplet size from slurry homogenization** → fine control over pore size
- **Maintains type-4 (capillary_bridge) pore formation** → clean spherical pores with capillary wetting
- **Minimizes core absorption** → RTF_mixer ≈ 0.5%, cores stay glassy
- **Enables chemistry-driven candidate screening** → poor solvents (RED > 1.1) are the sweet spot

**Recommended templates:** Butyl acetate (BA) or butyl butyrate (BB)  
**Target porosity:** 15–25% (fill = 0.10–0.25 pendular)  
**Expected pore size:** 0.3–0.5 μm (set by dispersion D_template)

---

## Chemistry Foundation: Hansen Solubility Parameters

### Polymer Hansen Sphere (Pinned)
| δD | δP | δH | r₀ | Hildebrand equiv. |
|----|----|----|-----|-------------------|
| 17.00 | 12.10 | 10.20 MPa^0.5 | 8.0 | 23.2 MPa^0.5 |

Source: User-supplied from characterization; cross-validated against xlsx-fitted sphere2 (closest fit).

### Template Screening: RED & Water Dispersibility

**Red = Ra / r₀** where Ra² = 4(ΔδD)² + (ΔδP)² + (ΔδH)²

Solvency bands:
- RED < 0.9: good solvent (theta point RED = 1.0)
- RED 0.9–1.1: theta (marginal)
- RED 1.1–1.5: poor (pore templates, type 4)
- RED > 1.5: very poor

Water dispersibility criterion: `water_sol < 5 g/100mL` OR `logP > 1.5`  
(Template must form stable droplets in aqueous slurry, not dissolve into continuous phase.)

### Recommended Candidates (Type 4, Capillary Bridge)

| template | RED | bp (°C) | logP | water_sol (g/100mL) | RTF (dryer) | pore_yield |
|---|---|---|---|---|---|---|
| butyl acetate | 1.20 | 126 | 1.78 | 0.7 | 0.12 | 88% |
| **butyl butyrate** | **1.33** | **166** | **2.83** | **0.04** | **0.035** | **96%** |

**Why BB over BA:**
- Lower RTF (dryer model): 0.035 vs 0.12 → less core absorption in atomizer dryer
- Better water immiscibility: 0.04 g/100mL vs 0.7 → droplet stability
- Higher bp: 166°C vs 126°C → slower evaporation (more control in dryer, needs inlet T ≥ 140°C)
- Pore yield: 96% vs 88% (minimal softening of polymer cores)

---

## Co-Feed Mechanism: UP1 to Spray Dryer

### Stage 1: UP1 Mixing (τ = 15–25 min, T = 20–50°C)

**Input state:**
- Colloid feed: Q_colloid = 3.3 mL/min (or scaled)
- Template feed: Q_template = 0.3–0.5 mL/min (10–15% of colloid flow)
- Template dispersion: Fine aqueous emulsion, D_template = 0.3–0.5 μm

**UP1 kinetics (type 4 capillary_bridge):**

Residual Template Fraction in mixer:
```
tau_D = r_particle² / D_polymer_eff  (diffusion timescale, ~30–120 min for glassy cores)
f_diff = 1 - exp(-3.0 * tau_eff / tau_D)  (fractional penetration into cores)
RTF_mixer = 0.05 * f_diff * fill  (type 4 capping, max 0.15)
```

At nominal conditions (D_particle ≈ 1 μm, tau ≈ 20 min):
- **RTF_mixer ≈ 0.005–0.010** (0.5–1% of template diffuses into cores)
- **99% stays free** as discrete liquid inclusions in interstitial water

**Key sensitivity results:**
- Droplet size D_template: **No effect on RTF_mixer** (diffusion is INTO polymer cores, not template droplets)
- Residence time: **Logarithmic response** (8× tau → only 2× RTF rise)
- Temperature: **Weak effect** (~50% RTF rise over 70°C span)

**UP1 exit state:**
- Polymer cores: Glassy, ~0.5% core saturation with template
- Interstitial space: 99% of template still as free liquid droplets
- Porosity retained: ~36% void fraction (full retention, Rumpf capillary physics)
- Bond strength (pendular): ~0.09–0.10 (moderate bonds from capillary wetting)

### Stage 2: Spray Dryer (Inlet T ≈ 140–160°C)

**Particle entry state (t = 0, inlet):**
- Template droplets: Encapsulated as liquid inclusions inside consolidating polymer shell
- Core state: Glassy (Tg ≈ 70°C >> T_inlet, no softening from template)
- Water: Rapid evaporation begins

**Drying sequence (t = 0–100 ms):**

1. **Shell consolidation (0–10 ms):** Polymer bonds form around template droplets; shell permeability begins to drop
2. **Template evaporation (10–50 ms):** Template bp exceeded → liquid → vapor phase transition inside pore space
   - BA (bp 126°C): Actively evaporating at inlet T
   - BB (bp 166°C): At threshold, slower evaporation (requires adequate dwell time)
3. **Vapor transport (10–100 ms):** Template vapor + water vapor diffuse through semi-permeable shell
   - **Critical:** Cores stay glassy → shell stays permeable → vapor escapes cleanly
   - Capillary network remains open (type 4 physics)
4. **Final cooling:** Vapor fully escaped, void remains open → pore retained

**Why type-4 succeeds (vs type-3 collapse):**
- Type-3 (good/theta solvent): RTF_mixer ≈ 0.20–0.40 → cores soften → Tg drops → shell collapses as polymer flows → pores flatten
- Type-4 (BA/BB poor solvent): RTF_mixer ≈ 0.005 → cores stay glassy → shell rigid → pores survive → open porosity

**Exit state (t > 100 ms):**
- Pore size: ~0.3–0.5 μm (set by D_template input to UP1)
- Porosity: φ_pore ≈ (1 − RTF_mixer) × fill ≈ 0.99 × fill
  - At fill = 0.20 → ~20% open porosity
  - At fill = 0.25 → ~25% open porosity
- Pore morphology: Spherical/ellipsoidal voids, capillary-connected (Rumpf wetting bridges)

---

## Implementation Roadmap

### Step 1: Prepare Template Dispersion
- **Solvent:** BA or BB (recommend BB for robustness)
- **Surfactant:** ~0.5–1% (e.g., Tween 80, SDS, or lecithin) to stabilize droplets
- **Homogenization:** Microfluidizer or high-shear mixer → D_template ≈ 0.3–0.5 μm
- **Storage:** Stable at room temperature (BB: very water-immiscible, minimal phase separation)

### Step 2: Co-Feed Setup
- Inject template dispersion into UP1 inlet, alongside colloid feed
- Flow rate: Q_template ≈ 10–15% of Q_colloid (for fill ≈ 0.15–0.20)
- Residence time: τ ≈ 15–25 min (no urgency to minimize; logarithmic benefit)
- Temperature: 20–50°C (weak effect on RTF_mixer; keep convenient)

### Step 3: Spray Dryer Conditions
- **Inlet temperature:** ≥140°C (BA: clear margin; BB: at threshold)
- **Outlet temperature:** ≥90°C (for complete water removal)
- **Dwell time:** ≥100 ms (allow template evaporation to complete)
- **Product collection:** Porous particles with engineered porosity (~15–25%)

### Step 4: Validation
- **SEM morphology:** Spherical pores, capillary-connected network
- **Porosity measurement:** BET (gas adsorption) or mercury porosimetry
- **Expected pore size distribution:** Narrow, centered at 0.3–0.5 μm
- **Batch consistency:** Pore size reproducible across atomizer runs

---

## Code Integration

### Files Updated/Created

1. **theta_solvent_chi.R** (primary module)
   - Hansen HSP calculation: `RED_from_hansen(dD, dP, dH, poly)`
   - Hildebrand chi (informational): `chi_from_hansen(dD, dP, dH, Vm, poly, T)`
   - Water dispersibility filter: `is_water_dispersible(water_sol, logP)`
   - Chemistry-to-regime predictor: `template_from_chemistry(template, poly, T_dry_C, db)`
   - Database: 29 species (polymer + 28 candidates), HSP + physical properties
   - Demo: Full candidate screen, type classification, volatility/water-solubility filtering
   - Added: butyl_butyrate, octyl_acetate (but octyl needs T≥170°C)

2. **up1_module_rev38_dryer_risk.r** (mixer kinetics)
   - Unchanged; existing kinetics already support type-4 capillary_bridge
   - RTF_mixer calculation uses LDF diffusion model (f_diff)
   - Type gating: rigid (1), gas (2), surface_weld (3), capillary_bridge (4)

3. **Sensitivity analyses** (scratchpad, for reference)
   - `ba_bb_cofeed.R`: Dose-sweep at nominal conditions
   - `ba_bb_sensitivity.R`: Parameter sensitivity (D_template, tau, T_system)
   - `pore_formation_mechanism.txt`: Detailed mechanistic walkthrough

### Integration Points

**To use co-feed in the unified chain:**

```r
source("theta_solvent_chi.R")

# Step 1: Chemistry-driven template selection
bb_props <- template_from_chemistry("butyl_butyrate", T_dry_C = 140)
# Returns: RED, chi, solvency, water_dispersible, escapes, regime, template_type=4,
#          RTF (dryer), rho_template, dD, dP, dH, logP, water_sol

# Step 2: Feed into UP1 with type=4 classification
up1_params <- list(
  template_type = bb_props$template_type,  # 4 (capillary_bridge)
  D_template = 0.5,                         # μm (dispersion droplet size)
  Q_template = 0.5,                         # mL/min
  C_temp_mass = 0.05,                       # 5% BB in aqueous dispersion
  template_dose = 0.20,                     # fill = 0.20 (pendular optimum)
  tau = 20, T_system = 70,                  # mixer conditions
  # ... other colloid, particle, surfactant params ...
)

# Step 3: Run UP1
up1_result <- run_single_simulation(up1_params)
RTF_mixer <- up1_result$Residual_Template_Fraction  # ~0.005–0.010

# Step 4: Spray dryer (existing UP4 module)
# Pore formation driven by (1 - RTF_mixer) ≈ 0.99 → 99% free template → pores
```

---

## Key Discoveries & Tradeoffs

| Finding | Implication |
|---|---|
| **RTF_mixer << D_template effect** | Droplet size doesn't gate RTF; use any D ∈ [0.3–0.5 μm] for consistency |
| **Type-4 chemistry dominates kinetics** | RED > 1.1 (poor solvent) suppresses core penetration more than tau/T |
| **Capillary bonds are weak (~0.09)** | Particles held by interfacial tension, not strong adhesion; acceptable for porous products |
| **BB pore yield > BA (96% vs 88%)** | Firmer poor solvent (RED 1.33 vs 1.20) → less core softening in dryer |
| **Water-immiscibility is critical** | BA water_sol = 0.7 g/100mL acceptable; BB at 0.04 g/100mL is more robust |
| **Porosity scales with fill** | φ_pore ≈ 0.99 × fill; for 20% porosity, use fill ≈ 0.20 (pendular) |

---

## Comparison to Alternatives

| Mechanism | Pore Size | Porosity | Droplet Control | Ease |
|---|---|---|---|---|
| **Co-feed (new)** | 0.3–0.5 μm | 15–25% | Excellent (D_template) | Moderate (dispersion prep) |
| Pre-emulsified slurry | 0.3–0.5 μm | 15–25% | Good (homogenization energy) | Easy (existing slurry) |
| Gas templating (type 2) | 1–10 μm | 20–40% | Moderate (bubble size) | Simple (co-inject gas) |
| Rigid filler (type 1) | ~10 μm | 10–20% | Poor (seed deposition) | Simple (add seeds) |
| Leachable filler (PMMA removal) | 1–10 μm | 20–40% | Excellent | Complex (2-step, leaching) |

---

## References & Data Files

**Documentation:**
- `review/theta_solvent_chi.md` — Hansen HSP fundamentals, candidate screen, regime map
- `COFEED_HANDOFF.md` (this file) — Co-feed mechanism, implementation roadmap

**Code:**
- `theta_solvent_chi.R` — Chemistry predictor, 29-species database
- `up1_module_rev38_dryer_risk.r` — Mixer kinetics (no changes needed for type-4)
- `scratchpad/ba_bb_cofeed.R` — Dose-sweep simulation
- `scratchpad/ba_bb_sensitivity.R` — Parameter sensitivity analysis

**Data files (xlsx, uploaded):**
- `references/HSP_Calculations_Ruben-Manuel.xlsx` — Abbott 1269-chemical database
- `references/NXTLEVEL Green solvents Solubility Parameters Jan 2023.xlsx` — Bio-based candidates

**Output files (generated):**
- `scratchpad/ba_bb_cofeed_sweep.csv` — Full dose-response table (BA, BB)
- `scratchpad/pore_formation_mechanism.txt` — Detailed stage-by-stage walkthrough

---

## Next Steps (Post-Handoff)

1. **Experimental validation:**
   - Prepare BB dispersion (0.5% Tween 80, homogenize to D ≈ 0.5 μm)
   - Run atomizer-dryer trial at inlet T = 150°C, fill = 0.20
   - Measure pore size (SEM, BET), porosity, and mechanical properties
   - Compare to pre-emulsified slurry baseline

2. **Scale-up:**
   - Adjust Q_template flow rate for desired fill (proportional to Q_colloid)
   - Confirm dispersion stability over production run (~1–2 hr)
   - Validate consistency batch-to-batch

3. **Optimization:**
   - Fine-tune D_template for target pore size (TEM SEM correlation)
   - Explore ternary dispersions (e.g., BA + BB blend) for pore-size distribution control
   - Thermal cycling tests (freeze–thaw stability of dispersion)

4. **Documentation:**
   - Update process specifications with co-feed SOP
   - Archive dispersion prep protocol (supplier, homogenizer settings, stability data)
   - Record SEM/BET results in data/ directory for audit trail

---

## Contact & Questions

- **Chemistry module:** `theta_solvent_chi.R` — Hansen HSP, RED, RTF (dryer model)
- **Kinetics module:** `up1_module_rev38_dryer_risk.r` — UP1 RTF_mixer, f_diff, type-4 gating
- **Implementation:** Co-feed dispersion prep, atomizer-dryer conditions, validation plan

*Handoff ready for experimental validation. All simulation results are reproducible via R scripts.*
