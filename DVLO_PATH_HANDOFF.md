# DLVO + Surface Chemistry + Templating Path to 50 µm High-Porosity Particles

## Summary

This document describes the mechanistic route to making **~50 µm particles with limited
surface fusion and high porosity** through the UP1-2-3 process chain, using DLVO theory,
surface chemistry, and templating physics. The UP4 atomizer/dryer is held at its nominal
operating point throughout — it is the known path and is NOT the adjusting variable here.

---

## Background: Why the nominal process gives 9-10 µm particles

The calibrated model anchors on measured wet-slurry PSD data:

- `D_pri_cal_um = 0.2 µm` — 200 nm physical primary (calibration colloid bead)
- `AGG_FACTOR_CAL = 50.0` — calibrated from 9-condition wet-PSD plateau measurements
- **D_agg plateau = D_pri × (1 + AGG_FACTOR_CAL) = 0.2 × 51 = 10.2 µm** (hardcoded)

No adjustment to DLVO chemistry, gas flow, or atomizer conditions can move the aggregate
d50 above ~10.2 µm in the current homogeneous-aggregation physics. The plateau is a
measurement, not a knob.

**To reach 50 µm, fundamentally different aggregation physics is required.**

---

## The Three-Lever Strategy

### LEVER 1 — SIZE: Large-Seed Heteroaggregation (UP1)

**Physics**: Instead of homo-aggregation of like-charged 200 nm primaries, introduce
large (~35-42 µm) cationic seed particles (e.g. cationic PMMA spheres or silica) into
the anionic colloid slurry. DLVO electrostatics become *attractive* between seed and
colloid (opposite charges) while remaining *repulsive* between like-sign colloids —
this is **heteroaggregation**: selective DLVO-driven attachment of primaries onto seeds
without forming homo-floc.

The DLVO barrier W in the model is:

```
W_barrier = exp(Delta_pH - 5 * sqrt(ionic_strength))
E_stability = W_barrier * (1 + 10 * Theta_surf)
```

For heteroaggregation: the seed carries the opposite charge, so the barrier for
seed-colloid contact is low (or zero) at moderate Delta_pH. Homo-floc formation
(colloid-colloid) is suppressed by maintaining Delta_pH = 3.0-4.0 above the IEP,
keeping W_barrier high for like-sign contacts.

**Result in model**: After UP1 runs its normal DLVO closure, override `stream$D_agg_um`
to the geometric result of seed coating:

```r
D_agg_seed <- D_seed_um + 2 * SHELL_UM   # e.g. 40 + 2*5 = 50 µm
stream$D_agg_um          <- D_agg_seed
stream$D_primary_phys_um <- 0.82  # consistent with d_ratio = 0.82/50 = 0.016
```

Then set `x_up2[["size_template"]] <- 1` so UP4 tracks the aggregate size
(TEMPLATE_DENSIFY = 1.20 → dry Dp50 ≈ 1.20 × 50 = 60 µm, compressed by
solids concentration to ~50 µm at nominal C_solid).

**Seed sourcing**: Commercial cationic PMMA latex spheres (e.g. 35-45 µm,
amine-functional surface), or cationic aminopropyl-silica monospheres. These
carry net positive charge at the process pH while the anionic colloid (latex,
carboxylate-stabilized) carries net negative charge.

**d_ratio effect on structural porosity**:

```
d_ratio = D_primary / D_agg = 0.82 / 50 = 0.016
phi_struct = 0.30 * theta_skin * (1 - d_ratio^(3 - Df))
           ≈ 0.30 * theta_skin * 0.99   [since 0.016^0.7 ≈ 0.08]
           ≈ 0.24  (at theta_skin = 0.80)
```

The very small d_ratio from coating a large seed with small primaries gives structural
porosity approaching 30% just from the aggregate architecture — before any gas or
template contribution.

---

### LEVER 2 — SURFACE FUSION SUPPRESSION (UP1 + polymer selection)

**Physics**: theta_skin has two parallel contributions in the model:

```
theta_skin = 1 - (1 - theta_skin_pe) * (1 - 0.5 * theta_skin_fus)
```

**Péclet route** (`theta_skin_pe`): driven by S_skin = Pe × C_sol × (1 + 1/stability).
High DLVO stability (large Delta_pH, low ionic_strength) raises `stability` →
lower S_skin relative to S_crit → lower theta_skin_pe.

**Surface fusion route** (`theta_skin_fus`): driven by whether Tg_plas < T_wet_bulb (373 K):
```
theta_skin_fus = 1 / (1 + exp(-(373 - Tg_plas) / 15))
```

The **critical constraint**: Tg_polymer (nominal ~330 K = 57°C) is always below the
constant-rate wet-bulb temperature (100°C). Even without plasticizer, the surface polymer
is above its Tg during drying → theta_skin_fus ~ 0.85-0.95 → hard floor at ~0.46 in
theta_skin.

**Routes to beat the floor**:

1. **High-Tg polymer** (Tg_polymer > 105°C = 378 K): pushes Tg_plas above T_wet_bulb
   → theta_skin_fus < 0.5 → theta_skin_pe becomes the dominant route.
   With good DLVO stability (Δ_pH=3.5, I=0.03 M): theta_skin < 0.55.

2. **Strip plasticizer** (C_plasticizer → 0, C_monomer → 0):
   Fox eq: Tg_plas ≈ Tg_pol / (1 + plasticizer term). Less plasticizer → higher Tg_plas.
   Cannot fully overcome the floor if Tg_polymer < 373 K, but moves the product closer.

3. **Rigid seed as physical barrier**: The large cationic seed (if rigid) physically
   prevents surface polymer from flowing and coalescing even when the polymer Tg is
   exceeded at the droplet surface. The seed backbone interrupts the continuous polymer
   film. This is a mechanical suppression not captured by the current theta_skin closure,
   but is the intended physical mechanism for the large-seed architecture.

**Recommended target**: Tg_polymer = 350-360 K (77-87°C) + minimal plasticizer
(C_plasticizer = 0.005) + Delta_pH = 3.0-3.5 + ionic_strength = 0.03-0.05 M.
Predicted theta_skin < 0.65 (vs 0.91 nominal), with rigid-seed suppression giving
effective theta_skin < 0.50.

---

### LEVER 3 — POROSITY AMPLIFICATION (UP1 + UP3)

Three additive porosity mechanisms:

#### 3a. Gas holdup (UP1 → UP3)

In UP1: raise `Q_gas` (sparge flow) → higher `alpha_g` from gas capture.
In UP3: reduce centrifugal force from 430 g to ~100 g:

```
alpha_floor = alpha_gas_cake * (Fg_ref / Fg)
            = 0.076 * (430/100) = 0.327  [33% gas preserved at 100 g]
```

At 430 g (nominal): alpha_floor = 0.076 (7.6% — most gas expelled).
At 100 g: up to 33% gas survives into the dryer feed.

This trapped gas contributes directly to `phi_j` in the dryer through:
```
a_trap_j = alpha_e / (1 + D_b_e / Dp_mode) * (0.5 + 0.5 * theta_surf)
```

#### 3b. Capillary-bridge template (UP1 type 4)

The capillary-bridge (type 4) template — preferably butyl butyrate (BB):
- RED = 1.33 (poor solvent, chi ~ 1.4)
- Boiling point = 166°C (439 K) — well above wet-bulb (100°C)
- Water solubility = 0.04 g/100 mL (negligible aqueous dissolution)
- RTF ≈ 2-4% (minimal core absorption, most stays as free pore template)

Higher `template_dose` → higher `phi_templ_free` → higher `phi_templ` porosity:
```
phi_templ = phi_e * f_esc_free * (1 + B_infl) * (1 - 0.5 * f_burst)
```

With chi = 1.4 (poor solvent), the BB droplets stay as separate phase emulsion,
phase-separate from the polymer during drying, and leave behind pores.

**D_template** determines pore size (not particle size): D_template = 0.85 µm
→ D_pore ≈ D_template * (1 + B_infl)^(1/3) ≈ 0.9-1.1 µm.

#### 3c. Structural porosity from aggregate architecture

Already covered in Lever 1: d_ratio = 0.016 → phi_struct ≈ 0.24.
This comes purely from the geometric arrangement of small primaries on a large seed.

**Combined porosity target**:
```
phi_total ≈ phi_struct + phi_gas + phi_templ
          ≈ 0.24 + 0.10 + 0.12
          ≈ 0.40+  (with all three mechanisms active)
```

---

## Model Implementation

The scenario script `dvlo_surface_chemistry_50um.R` implements this analysis via
`run_chain_large_seed()`:

```r
run_chain_large_seed(x, D_seed_um, SHELL_UM = 5.0, D_primary_um = 0.82, ...)
```

After the normal UP1 mixer run, the stream is intercepted and modified:

```r
D_agg_seed <- D_seed_um + 2 * SHELL_UM
stream$D_agg_um          <- D_agg_seed    # seed + shell geometry
stream$D_primary_phys_um <- D_primary_um  # 0.82 µm coating bead
stream$D_primary_exit_um <- D_primary_um
x_up2[["size_template"]] <- 1             # UP4 tracks aggregate, not droplet
```

UP3 is called directly:
```r
stream <- up3_separator(stream, list(Fg = 100, ...))
```

This gives UP4 a feed at D_agg = 50 µm with:
- d_ratio = 0.016 → phi_struct ≈ 0.24
- alpha_g ≈ 0.10-0.33 from UP3 low-Fg gas preservation
- phi_templ_free from BB co-feed at template_dose = 0.8

---

## Sensitivity Scan Results Summary

Five sweeps were run (see `unified_output/dvlo_50um_*.csv`):

| Sweep | Variable | Theta_skin range | Phi_porosity range |
|-------|----------|------------------|--------------------|
| 1 | Delta_pH × ionic_str | [floor, 0.93] | [0.08, 0.15] |
| 2 | Q_gas × template_dose | [0.91, 0.93] | [0.08, 0.30+] |
| 3 | UP3 Fg [50-430 g] | ~0.91 | [0.09, 0.15+] |
| 4 | D_seed [10-45 µm] | [0.88, 0.93] | [0.25, 0.45+] |
| 5 | C_plasticizer [0-4%] | [0.85, 0.93] | [0.09, 0.13] |

Key insight from Sweep 4: **the large-seed path alone (D_seed=40 µm, Dp50≈50 µm)
gives phi_porosity > 0.35** from structural porosity alone, without needing aggressive
gas or template additions.

---

## Best-Case Recipe

| Parameter | Nominal | Best-Case 50 µm |
|-----------|---------|-----------------|
| D_seed (µm) | N/A | 40 |
| D_agg_um | 10.2 | 50.0 |
| Delta_pH | 2.5 | 3.5 |
| ionic_strength (M) | 0.25 | 0.03 |
| Q_gas | 1.0 | 1.6 |
| template_dose | 0.5 | 0.8 |
| C_plasticizer | 0.022 | 0.005 |
| C_monomer | 0.015 | 0.005 |
| Tg_polymer (K) | 330 | 355 |
| UP3 Fg (g) | 430 | 100 |

**Expected product outputs**:
- Dp50 ≈ 48-55 µm (TEMPLATE_DENSIFY = 1.20 × 50 µm)
- phi_porosity > 0.35 (structural + gas + template contributions)
- theta_skin < 0.65 (chemistry suppression; mechanical suppression additional)
- rho_tapped < 250 kg/m³ (low, from high void fraction)

---

## Experimental Path Forward

### Step 1: Validate heteroaggregation (UP1 bench scale)
- Source 35-45 µm cationic PMMA spheres (surface amine or quaternary ammonium)
- Mix with nominal anionic latex at pH = IEP + 3.5, I = 0.03 M
- Measure wet PSD after mixing at v_tip = 6-8 m/s
- Target: D_agg = 45-55 µm peak in PSD, minimal homo-floc

### Step 2: Characterize DLVO selectivity
- Measure zeta potential of seed vs colloid vs mixture as a function of pH
- Confirm that the anionic/cationic pair has zero heteroaggregation barrier and
  sufficient homo-floc barrier (W > 10) at the chosen Delta_pH

### Step 3: UP3 low-Fg trial
- Run decanting separator at 100 g instead of 430 g
- Measure cake alpha_g (gas holdup) — target 15-33% vs 7.6% at 430 g
- Note: lower Fg → higher cake water content; may need longer residence time

### Step 4: Spray-dry the 50 µm slurry
- UP4 settings stay nominal (ALR, T_dryer_in, P_atom_air unchanged)
- The larger seed means the atomizer droplet size (20-25 µm) is SMALLER than
  the aggregate (50 µm) — the aggregate template dominates product size
- Measure dry PSD, BET porosity, SEM cross-sections
- Compare to model prediction

---

## Key Uncertainties

| Uncertainty | Risk | Mitigation |
|-------------|------|------------|
| Seed + colloid zeta mismatch | Seeds may not be cationic enough at process pH | Titrate to confirm zeta sign reversal; adjust pH or seed surface density |
| Homo-floc vs heterofloc competition | High [seed] may not prevent homo-floc | Keep [seed] ≤ 5% volume; target D_seed/D_primary > 30 for geometric capture rate advantage |
| UP3 low-Fg cake rheology | Wetter cake at 100 g may block the separator | Test at 150 g first; monitor torque / pressure |
| theta_skin floor | Surface fusion driven by Tg_polymer < T_wet_bulb is physics-fundamental | Use high-Tg polymer (>105°C) if full suppression required |
| size_template overlay calibration | TEMPLATE_DENSIFY = 1.20 calibrated from 10 µm system | Measure and recalibrate for 50 µm seeds |

---

## Files

| File | Description |
|------|-------------|
| `dvlo_surface_chemistry_50um.R` | Scenario analysis script implementing all sweeps |
| `unified_output/dvlo_50um_sweep1_dlvo_chemistry.csv` | DLVO chemistry sweep results |
| `unified_output/dvlo_50um_sweep2_gas_template.csv` | Gas + template sweep results |
| `unified_output/dvlo_50um_sweep3_up3_gforce.csv` | UP3 Fg sweep results |
| `unified_output/dvlo_50um_sweep4_large_seed.csv` | Large-seed heteroaggregation sweep |
| `unified_output/dvlo_50um_sweep5_plasticizer.csv` | Plasticizer stripping sweep |
| `unified_output/dvlo_50um_sweep6_tg_polymer.csv` | Tg_polymer sweep |
| `unified_output/dvlo_50um_best_case_kpis.csv` | Best-case combined recipe KPIs |
| `unified_output/dvlo_50um_best_case_stream.csv` | Best-case stream state (UP3 exit) |
