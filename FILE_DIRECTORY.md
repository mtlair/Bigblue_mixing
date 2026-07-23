# Bigblue_mixing Repository: Co-Feed Mechanism Implementation

**Last Updated:** 2026-07-23  
**Status:** Merged to main, ready for experimental validation  
**Current Branch:** `claude/main-handoff-review-0qo61a` (working branch, synced with main)

---

## Primary Chemistry Module

### `theta_solvent_chi.R` (NEW/UPDATED)
**Purpose:** Chemistry-driven template predictor using Hansen Solubility Parameters

**Key Functions:**
- `RED_from_hansen(dD, dP, dH, poly)` — Hansen relative energy distance
- `chi_from_hansen(dD, dP, dH, Vm, poly, T_K)` — Flory-Huggins interaction parameter (informational)
- `is_water_dispersible(water_sol, logP)` — Water dispersibility filter (< 5 g/100mL OR logP > 1.5)
- `template_from_chemistry(template, T_dry_C, db)` — Regime predictor (type 1/2/3/4)
- `screen_hsp_database(xlsx_path, RED_min, RED_max)` — Abbott 1269-chemical screen

**Database:** 29 species (polymer + 28 candidates)
- Original set: water, styrene, MMA, toluene, hexane, cyclohexane, acetone, ethanol, DBP, paraffin, CO2, d-limonene, butyl acetate, chloroform, methylene chloride, TEC, IPM
- Expanded set: ethyl butyrate, isoamyl acetate, amyl acetate, hexyl acetate, dibutyl sebacate, chlorobenzene, xylene, ethyl levulinate, butyl levulinate
- **NEW co-feed candidates:** butyl_butyrate, octyl_acetate

**Output:** template_type (1–4), RED, chi, RTF, solvency, water_dispersible, escapes, regime, mechanism

**Run demo:** `Rscript theta_solvent_chi.R`

---

## Mixer Kinetics Module

### `up1_module_rev38_dryer_risk.r` (EXISTING, NO CHANGES)
**Purpose:** UP1 spray-dryer mixing kinetics and template diffusion into polymer cores

**Key Output Parameters:**
- `RTF_mixer` (Residual Template Fraction in mixer): ~0.5–1% for type-4 capillary_bridge
- `Bond_Strength` (capillary wetting, Rumpf model): ~0.09–0.10 at pendular fill
- `Retained_Porosity`: ~0.36 (full void retention for type-4)
- `Template_Target_Score`: Bond × Porosity / EPS_RCP
- `Swelling_Softness_exit`: Polymer network openness (fed to UP2 dryer)

**Type-4 Gating:** RTF_mixer = 0.05 × f_diff × fill (capped 0–0.15)

**Demo sweep:** Lines 679–765 (pendular dose window, template_type sweep)

---

## Documentation & Handoff

### `COFEED_HANDOFF.md` (NEW, PRIMARY)
**Purpose:** Complete summary of co-feed mechanism for implementation & validation

**Sections:**
1. Executive summary
2. Hansen HSP chemistry foundation
3. Co-feed mechanism (UP1 → dryer stages)
4. Implementation roadmap (dispersion prep, flow setup, validation)
5. Code integration guide
6. Key discoveries & tradeoffs
7. Comparison to alternative templating
8. References & next steps

**Table of Contents:**
- Recommended templates: Butyl butyrate (BB, RED=1.33) vs butyl acetate (BA, RED=1.20)
- Stage-by-stage mechanism: UP1 mixing → shell consolidation → template evaporation → pore formation
- Practical setup: Inlet T ≥ 140°C, fill 0.15–0.25 pendular, D_template 0.3–0.5 μm

---

### `review/theta_solvent_chi.md` (EXISTING, UPDATED)
**Purpose:** Technical documentation of Hansen HSP method, Hildebrand fallacies, regime map

**Sections:**
- Method: Hansen Solubility Sphere, RED calculation, theta point definition
- Why Hansen (not single Hildebrand): DBP, MMA, toluene comparison
- Water-dispersibility filter
- Regime map (good/theta/poor/gas → types 3/NA/4/2)
- Candidate screen (27 species): type 4 pore templates, plasticizers, rigid fillers, non-viable
- What's needed: r₀ calibration, absolute pore size/porosity calibration

---

## Sensitivity & Validation

### `scratchpad/ba_bb_cofeed.R` (NEW)
**Purpose:** Dose-sweep simulation at nominal UP1 conditions (D_template=0.5 μm)

**Results:** BA and BB produce identical RTF_mixer (~0.005–0.025 across dose range)
- Pendular optimum: dose ≈ 0.795 (high, near funicular boundary)
- RTF_mixer at fill=0.20: ~0.0076
- Interpretation: 99.2% of template stays free (99.2% → pores)

**Output:** `scratchpad/ba_bb_cofeed_sweep.csv` (40 dose points × 2 templates)

---

### `scratchpad/ba_bb_sensitivity.R` (NEW)
**Purpose:** Parameter sensitivity analysis (D_template, tau, T_system)

**Key Findings:**
1. **Droplet size has NO effect on RTF_mixer** (diffusion INTO polymer cores, not template droplets)
2. **Residence time:** Logarithmic response (8× tau → only 2× RTF)
3. **Temperature:** Weak effect (~50% RTF rise over 70°C span)

**Interpretation:**
- Type-4 chemistry (RED > 1.1) dominates RTF suppression
- Physical droplet size ≠ kinetic factor
- Use any D_template ∈ [0.3–0.5 μm]; chemistry gates the RTF, not droplet size

---

### `scratchpad/pore_formation_mechanism.txt` (NEW)
**Purpose:** Detailed mechanistic walkthrough of template → pore transformation in spray dryer

**Stages:**
1. Inlet conditions (free template, glassy cores, RTF_mixer ≈ 0.5%)
2. Particle atomisation & early drying (0–10 ms): encapsulation, shell consolidation
3. Core softening check (glassy cores remain, RTF_mixer ≠ spray-dryer RTF)
4. Template evaporation (10–50 ms): liquid → vapor inside pore space
5. Vapor transport & water co-evaporation (10–100 ms): escape through semi-permeable shell
6. Final structure: spherical pores, capillary-connected network

**Type-3 vs Type-4 Comparison:**
- Type-3 (good/theta, e.g., MMA RED=0.95): RTF_mixer 0.20–0.40 → cores soften → shell collapses → pores flatten
- Type-4 (poor, e.g., BA/BB RED=1.2–1.3): RTF_mixer 0.005–0.010 → cores stay glassy → shell rigid → pores survive

---

## Reference Data Files

### `references/HSP_Calculations_Ruben-Manuel.xlsx`
**Source:** Steven Abbott HSP Database (1269 chemicals)

**Sheets:**
- Chemicals: CAS, name, dD, dP, dH, Vm
- Pre-fitted polymer spheres: 1-sphere and 2-sphere (Sphere2 closest to user-supplied values)

**Usage:** `screen_hsp_database()` in theta_solvent_chi.R queries this for candidates in RED ranges

---

### `references/NXTLEVEL Green solvents Solubility Parameters Jan 2023.xlsx`
**Source:** NxtLevel bio-based green solvent database

**Proprietary Solvents:** Ethyl levulinate (EtLev), Butyl levulinate (ButLev), EtLPK, EtLGK, EHL
- HSP provided; bp/water_sol/logP for most unknown (needs supplier confirmation)
- EtLev RED=0.75 (good, non-volatile, plasticizer, not pore template)
- ButLev RED=0.71 (good, non-volatile, bio-based plasticizer option)
- EtLPK/EHL: RED 1.12–1.15 (pore-template range if bp < T_dry+40)

---

## Output & Analysis

### `scratchpad/ba_bb_cofeed_sweep.csv`
**Columns:** template, dose, Bond_Strength, Retained_Porosity, Template_Target_Score, Residual_Template_Fraction, Blended_Porosity

**Rows:** 40 dose levels × 2 templates (BA, BB)

**Key metric:** RTF_mixer vs fill (dose); used to confirm type-4 porosity retention across pendular window

---

## Integration Checkpoints

### Code Flow
```
theta_solvent_chi.R
  ├─ polymer Hansen sphere (dD=17, dP=12.1, dH=10.2, r0=8)
  ├─ candidate database (29 species)
  ├─ Hansen RED calculator
  ├─ water dispersibility filter
  └─ template_from_chemistry() returns type/RTF/regime
       ↓
up1_module_rev38_dryer_risk.r
  ├─ type=4 (capillary_bridge) gating
  ├─ RTF_mixer = LDF diffusion into polymer cores
  ├─ Bond_Strength = Rumpf capillary adhesion
  ├─ Retained_Porosity = pendular dose window
  └─ output: Bond, Porosity, RTF_mixer, Swelling_Softness_exit
       ↓
spray_dryer (UP2, existing)
  ├─ inlet T ≥ 140°C
  ├─ template evaporates (bp-controlled)
  ├─ vapor escapes (glassy cores keep shell permeable)
  └─ final pores: φ_pore ≈ (1-RTF_mixer) × fill ≈ 20% at fill=0.20
```

---

## Validation Plan

### Experimental (Post-Handoff)

1. **Dispersion prep:**
   - Butyl butyrate: 5% w/w in water + 0.5% Tween 80 (or equivalent surfactant)
   - Homogenize to D ≈ 0.5 μm (microfluidizer or high-shear mixer)
   - Verify stability over 2 hr

2. **UP1 co-feed trial:**
   - Q_colloid = 3.3 mL/min
   - Q_template = 0.3–0.5 mL/min (fill ≈ 0.15–0.20)
   - tau ≈ 15–25 min
   - Measure viscosity, colloid stability, RTF_mixer (optional ODE simulation)

3. **Spray dryer:**
   - Inlet T = 150°C (>BB bp threshold 166°C barely met; may need 160°C for robust evaporation)
   - Outlet T ≥ 90°C
   - Dwell time ≥ 100 ms
   - Collect product

4. **Pore characterization:**
   - **SEM:** Morphology, pore size distribution, surface connectivity
   - **BET:** Specific surface area, porosity (15–25% expected)
   - **Mercury porosimetry:** Pore size, capillary throat distribution
   - **Mechanical:** Tensile, flexural (check capillary-bonded network strength)

5. **Batch consistency:**
   - Repeat 2–3 trials
   - Monitor pore size (SEM) and porosity (BET) reproducibility

---

## Quick Reference

**Co-Feed at a Glance:**
- **Mechanism:** Template-in-water dispersion co-fed into UP1 mixer
- **Chemistry:** Butyl butyrate (RED=1.33, bp=166°C, water_sol=0.04 g/100mL)
- **UP1 outcome:** RTF_mixer ≈ 0.5% → 99% free template
- **Dryer outcome:** Template evaporates → spherical pores (0.3–0.5 μm)
- **Expected porosity:** 15–25% (pendular fill 0.15–0.25)
- **Target inlet T:** ≥ 140°C (preferably 150–160°C for BB)

**Files to Run:**
- `Rscript theta_solvent_chi.R` — Verify chemistry predictor & candidates
- `Rscript scratchpad/ba_bb_sensitivity.R` — Confirm UP1 parameter sensitivities
- `up1_module_rev38_dryer_risk.r` (within R session) — Run full UP1 simulation with type=4

**Documentation:**
- `COFEED_HANDOFF.md` — Implementation guide & mechanism explanation
- `review/theta_solvent_chi.md` — Hansen HSP fundamentals
- `scratchpad/pore_formation_mechanism.txt` — Detailed physics walkthrough

---

## Branch & Repository Status

**Working Branch:** `claude/main-handoff-review-0qo61a`  
**Main Branch:** Synced, includes COFEED_HANDOFF.md + theta_solvent_chi.R updates  
**Remote:** Pushed to origin; no uncommitted changes

**Commits (this session):**
1. Add butyl butyrate and octyl acetate to HSP database
2. Fix escapes logic (bp-based, not volatile flag)
3. Add co-feed handoff documentation

---

*Prepared by: Claude Sonnet 4.6 | Session: https://claude.ai/code/session_01NEqQDtgpcYNVBFZmLwPGS6*
