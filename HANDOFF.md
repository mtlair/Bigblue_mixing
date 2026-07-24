# Bigblue Mixing — Session Handoff

**Date:** 2026-07-23  
**Branch:** `claude/psd-viscosity-data-review-khqg39`  
**Repo:** `mtlair/Bigblue_mixing`

---

## What was done this session

### 1. theta_skin closure fix (`unified/up2_atomizer_dryer_module.R`)

The original surface_fusion-formation closure had the wrong sign for plasticizer/monomer
effects (more solvent → less surface_fusion, which is backwards). Fixed by adding a
**surface-fusion route** (Fox Tg depression at constant-rate wet-bulb 100 °C) and
correcting the Péclet route so softness *lowers* the aggregation threshold:

```r
# Péclet route  — softness LOWERS resistance threshold
S_crit_pe     <- 500 / (1 + 0.05 * Softness)
theta_skin_pe <- S_skin / (S_skin + S_crit_pe)

# Surface-fusion route — Fox equation + wet-bulb surface temp
T_surface_cr  <- 373.15   # K (constant-rate period)
inv_Tg_plas   <- (1 - phi_solvent) / Tg_polymer + phi_solvent / Tg_solv
theta_skin_fus <- 1 / (1 + exp(-(T_surface_cr - (1/inv_Tg_plas)) / 15))

# Parallel routes, series non-surface_fusion probabilities
theta_skin <- 1 - (1 - theta_skin_pe) * (1 - 0.5 * theta_skin_fus)
```

After the fix C_monomer/C_plasticizer drive MORE surface_fusion (ranks 14/16, positive μ*),
up from ranks 16–17 with the wrong direction. All other factor rankings unchanged.

### 2. Morris sensitivity screen re-run

Deleted stale cache `unified_output/unified_morris_t4_r30.rds`, re-ran
`unified_model.R` (41 factors, r = 30 trajectories → 1260 evaluations).
1260/1260 valid. Top-5 chain drivers per key output:

| Output | Drivers (rank order) |
|--------|----------------------|
| D_particle | P_atom_air[up2], tau_mix[up1], C_solid_mass[up1], sigma[up2], T_mix[up1] |
| phi_porosity | C_solid_mass[up1], template_dose[up1], Q_template[up1], T_dryer_in[up2], Q_colloid[up1] |
| rho_tapped | P_atom_air[up2], tau_mix[up1], template_dose[up1], C_solid_mass[up1], sigma[up2] |
| Tg_product | Tg_polymer[up2], C_solid_mass[up1], C_plasticizer[up1], Q_template[up1], Q_colloid[up1] |
| X_moisture | C_solid_mass[up1], Q_template[up1], Q_colloid[up1], k_perm_bind[up2], C_monomer[up1] |

Morris output in `unified_output/` — see `morris_run.log` for validation counts.

### 3. Full-chain data processing (`data/up1_2_3_4_visc_sem.xlsx`)

Parsed 4 complete UP1→UP2→UP3→UP4 chain conditions (`up3_1`–`up3_4`) and the
new `data/up1_2_3_solids_rho.xlsx`. Derived quantities written to
`data/cond_up1234.csv`. Key findings:

**Solid particle density: 1730–1742 kg/m³**  
Inferred from feed density measurement (`1/ρ_feed = s/ρ_solid + (1−s)/ρ_water`).
Higher than pure polymer (~1200) because particles are mineral-loaded. Slightly
lower than pure mineral skeleton due to residual monomer — consistent across all
xlsx files.

**UP2 water routing confirmed by mass balance**  
UP3 mass balance closes against UP1 output alone (249.84 ≈ 152.8 + 97.2 lb/hr).
The UP2 wash water (49 lb/hr, column `up2_h2o_lbhr`) exits UP2 as a separate
effluent — it does NOT pass through the separator. UP3 therefore receives
the UP1 foam unchanged in total flow; it only concentrates from 25 % → 40 % solid.

**UP3 cake gas holdup: 7.6 %** (for up3_1, the only measured condition)  
Using ρ_solid = 1730 kg/m³: gas-free cake density at 40.7 % solid = 1209 kg/m³.
Measured = 1117.3 kg/m³ → gas holdup = (1209−1117)/(1209−1) = 7.6 %. Real
residual trapped gas from UP1 sparging surviving centrifugation at 430 g.

**UP1 exit density (two estimates)**  

| Method | alpha_g at exit | ρ_exit (kg/m³) | Basis |
|--------|----------------|----------------|-------|
| All sparge gas entrained (upper bound) | 30–31 % | 770–790 | ideal-gas expansion of up1_scfh |
| Back-calc from UP3 7 % holdup (Boyle) | 4.4–4.9 % | 1063–1070 | compressed 7 % atm → operating P |

Model should use the back-calc estimate; the "all-entrained" value assumes no gas
escapes the headspace before product exits.

**UP4 feed rate**  

| Cond | up4_feed (mass balance, lb/hr) | up4_feed (heat duty, lb/hr) | Dryer efficiency |
|------|-------------------------------|----------------------------|-----------------|
| up3_1 | 152.8 | 234.9 | 65 % |
| up3_2 | 123.3 | 236.4 | 52 % |
| up3_3 | 129.7 | 212.8 | 61 % |
| up3_4 | 135.2 | 195.7 | 69 % |

Mass-balance estimate is primary. Heat-duty estimate is secondary (variable dryer
thermal efficiency 52–69 %, probably because `up4_dry_air_scfm` was operator-fixed
at a constant high setting regardless of throughput).

---

## File directory

```
Bigblue_mixing/
├── DATA_REVIEW.md          ← primary narrative; Parts 1–9 complete
├── HANDOFF.md              ← this file
│
├── data/
│   ├── visc.xlsx               original UP1/UP4 direct-path data (9 conditions)
│   ├── up1_2_3_4_visc_sem.xlsx full-chain data (4 conditions, PSD + process vars)
│   ├── up1_2_3_solids_rho.xlsx historical run set with UP3 density + reject data
│   ├── cond_up1234.csv         ← parsed + derived full-chain table (primary output)
│   ├── up3_viscometry.csv      post-UP3 flow curves (4 conditions)
│   └── cond_process.csv        direct-path 9-condition process table
│
├── unified/
│   ├── up1_mixer_module.R      UP1 gassed mixer model
│   ├── up2_atomizer_dryer_module.R  ← theta_skin fixed here (two-route closure)
│   └── interface_stream.R      stream state definitions (SKIN_SURVIVAL = 0.30)
│
├── unified_model.R             top-level chain runner + Morris screen
├── validate_up1_up4_direct.R   validation harness for 9 direct-path conditions
│
├── unified_output/
│   ├── unified_morris_t4_r30.rds    Morris cache (regenerated after surface_fusion fix)
│   ├── unified_morris_indices.csv   41-factor μ* and σ rankings
│   ├── unified_morris_{process,surface,polymer}_variables.png
│   ├── nominal_chain_summary.txt    nominal-point chain outputs
│   ├── nominal_chain_outputs.csv
│   ├── morris_run.log               1260/1260 valid
│   ├── up1_up4_direct_validation.csv
│   └── direct_morris_overlays.csv
│
└── sem/                        SEM images (sample numbers 114641–114644)
```

---

## What is NOT yet done (next steps)

1. **UP1→UP4 full-chain validation harness** — the 4 conditions in
   `cond_up1234.csv` need a `validate_up1_up4_chain.R` script analogous to
   `validate_up1_up4_direct.R`, routing UP1→UP2→UP3→UP4 with the correct
   post-UP3 feed state (40 % solid, ρ ≈ 1117 kg/m³, η ≈ 12–27 Pa·s).

2. **UP2 foam-wash model** — currently a pass-through. Mass balance is closed but
   the impurity-removal function of UP2 (monomer washout, template cleaning) is
   not modelled.

3. **UP3 separator model** — currently a pass-through. Could be modelled as a
   solid-liquid separator with a G-force–dependent moisture closure (weir height
   also affects cake wetness; `up3_weir` = hi/low is a control variable).

4. **Reject data for up3_2–up3_4** — only up3_1 has `up3_reject_solid%` and
   `up3_reject_lbhr` measured. The other three are null. Reject loss is small
   (~0.28 % solid) so neglecting it introduces < 0.5 % error in UP4 feed.

5. **SEM image analysis** — images in `sem/` have not been processed. Sample
   numbers 114641–114644 map to up3_1–up3_4 (see `sem_sample_no_filename` sheet).
   Filename convention: `SEM_X <sample_no>_<a/b/null><magnification><k/null>x0`.

6. **DoE proposal** — Morris screen identifies top factors per output property.
   Next step is to translate factor rankings into an experimental design with
   bounds derived from the calibrated model.

---

## Calibration status

| Metric | Calibration ratio | Status |
|--------|------------------|--------|
| Wet d50 (UP1-control 8) | 1.07 | ✓ well-calibrated |
| Viscosity (UP1-control 8) | 2.79 | acceptable |
| Dry d50, muted (8) | 1.68 | ✓ |
| Dry d50, tmpl=1 (8) | 1.24 | ✓ |

Surface_fusion closure fix did not change the PSD calibration metrics (theta_skin feeds
morphology, not particle size). The full-chain calibration (through UP3) is the
remaining open item.
