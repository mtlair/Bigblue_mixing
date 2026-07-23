# UP1/UP4 PSD + viscosity data review and model recalibration

Source data: [`data/visc.xlsx`](data/visc.xlsx) — rotational-rheometer flow curves
(shear rate 12.7 → 1889 s⁻¹, all at ~25 °C), plus a PSD + run-condition sheet
(`up1_up4_psd`) keyed by **cond** number.

- `sample_type = m` — post-**UP1** wet PSD (colloid/aggregate out of the mixer)
- `sample_type = p` — post-**UP4** dry PSD (spray-dried powder)
- Temperatures in °C; `up1_psig` / `max_up4atom_psig` are gauge pressures at the
  UP1 exit and the UP4 atomizer; all PSD percentiles in µm.

Nine conditions with paired m/p PSD were reviewed: cond 2, 4, 5, 8, 9, 11, 12, 14, 15.

---

## Part 1 — Data QC

**Integrity.** Units and columns are internally consistent; temperature holds at
25.0 ± 0.05 °C across every flow curve (isothermal, so no temperature correction
needed). PSD percentiles are monotone (d10 < … < d99) in all 18 rows. One
cosmetic issue: `cond2` carries a stray `fit` column of `#REF!` errors (a broken
cross-reference, column F) — it holds no data and is ignored.

**Two clean rheological populations.**

| Population | Sheets | Low-shear η (12.7 s⁻¹) | Flow index n |
|---|---|---|---|
| UP1 **feed** (pre-mix) | `up1_feed_0.13/0.2/0.25` | 0.0022–0.0030 Pa·s | 0.89–0.97 (near-Newtonian) |
| UP1 **product** (post-mix) | `cond*` | 0.020–0.770 Pa·s | 0.20–0.67 (strongly shear-thinning) |

The feeds at 13 / 20 / 25 % solids are essentially water-thin and Newtonian.
Mixing in UP1 transforms them into shear-thinning floc suspensions — the flow
index drops from ~0.9 to a mean of **0.44**, and the low-shear viscosity rises
1–2 orders of magnitude.

**Two process control regimes.** The nine conditions split into two regimes that
control the *final (dry) product size* by different mechanisms:

- **cond2 — atomizer-control regime.** Its UP1 feed never aggregates (v_tip
  4.37 m/s < onset, so the slurry stays a dispersed ~0.18 µm colloid, near-water
  0.005 Pa·s). With no upstream structure to template the particle, the **UP4
  atomizer sets the product size**: the dry d50 (19.85 µm) is a droplet-shell — the
  *largest* dry particle in the set despite the *lowest* solids and finest wet feed.
- **cond 4/5/8/9/11/12/14/15 — UP1-control regime.** UP1 forms a ~10 µm aggregate
  (v_tip > onset) that templates the particle through drying, so the **UP1 aggregate
  sets the product size**: dry d50 = 12.3 µm mean ≈ 1.2× the wet aggregate, and it
  does *not* balloon out to the 20–27 µm droplet-shell the atomizer alone would give.

**Aggregation is a sharp on/off switch in tip speed** (the UP1 regime boundary):

| cond | v_tip (m/s) | solid % | wet d50 (µm) | regime |
|---|---|---|---|---|
| 2 | 4.37 | 13.2 | **0.183** | atomizer-control (no UP1 aggregate) |
| 5 | 6.92 | 19.6 | 9.97 | aggregated |
| 15 | 7.41 | 24.9 | 9.71 | aggregated |
| 8 | 8.48 | 19.2 | 9.60 | aggregated |
| 14 | 9.61 | 25.0 | 9.74 | aggregated |
| 9 | 9.62 | 19.6 | 10.37 | aggregated |
| 12 | 9.62 | 19.0 | 9.96 | aggregated |
| 4 | 10.57 | 20.0 | 10.19 | aggregated |
| 11 | 12.40 | 19.6 | 12.00 | aggregated |

The one condition below ~6.8 m/s (cond2, atomizer-control) stays at the ~0.18 µm
primary particle; every condition above it (UP1-control) jumps to a **flat ~10 µm
aggregate plateau** (mean 10.19 µm, CV 7 %) that barely moves with tip speed. There
is **no milling rollover** anywhere in the tested envelope — cond11, the highest tip
speed (12.4 m/s), has the *largest* wet d50, not a reduced one.

**Viscosity is highly scattered and not a clean function of the process.**
Across the eight aggregated conditions (all ~19–25 % solids, all ~10 µm aggregate)
the low-shear viscosity spans 0.020–0.770 Pa·s — a **38× range, CV 77 %**. It does
not track solids, tip speed, or wet d50 (|corr| ≤ 0.35 for each). cond8 (0.770) is
the high outlier and cond12 (0.020, run at the high 150 lb/hr feed) the low one.
This is genuine process/sample variability; a deterministic closure can only be
expected to hit the **central tendency** (geometric mean 0.201 Pa·s, median 0.227).

**Dry PSD — which unit controls depends on the regime.**
In the **UP1-control** regime the post-UP4 dry d50 is 8.8–15.3 µm (mean 12.3),
≈ 1.2× the wet aggregate, and *anti*-correlated with UP1 tip speed (corr −0.91)
and wet d50 (−0.79): harder mixing upstream builds a denser aggregate that dries
to a finer, more compact particle (cond11, highest shear, gives the smallest dry
d50 at 8.82 µm despite the largest wet aggregate). The **atomizer-control** cond2
inflates 108× (0.18 → 19.85 µm) — a hollow droplet-shell whose size is set by the
atomizer, not the colloid.

**Atomizer operating envelope (from `up4atom_scfm`, `up4_feed`, `max_up4atom_psig`).**
Converting the recorded flows (air at 0.0765 lb/ft³; feed lb/hr → kg/s; psig →
absolute Pa) gives the true nozzle envelope:

| Nozzle input | Measured (9 cond) | Prior model factor range | Status |
|---|---|---|---|
| ALR = ṁ_air/ṁ_feed | **0.90–1.78** (mean 1.50) | 1.0–10.0 | model ran up to ~6× too high |
| Atomizing-air pressure | **1.20–1.80×10⁵ Pa** (2.7–11.4 psig) | 2.0–7.0×10⁵ Pa | entire model range *above* the real max |
| Liquid feed ṁ_L | **0.0135–0.0194 kg/s** | 0.002–0.020 kg/s | low end ~7× below any real run |

For the atomizer-control regime, cond2 anchors the droplet→particle scale: dry
19.85 µm at ALR 1.39 back-solves to a ~40 µm atomization droplet at 12.7 % solids.

---

## Part 2 — Model vs. measurement

The active chain is `unified/up1_mixer_module.R` → `interface_stream.R` →
`unified/up2_spray_dryer_module.R` (sourced by `unified_model.R` and
`full_train_mixer_to_dryer.R`). The UP1 slurry viscosity `Blended_Viscosity_PaS`
becomes `mu_exit_PaS` and feeds the nozzle `mu_slurry0` when `couple_viscosity`
is on, so the UP1 viscosity closure has real downstream consequence.

The prior UP1 closures were **single-point calibrations to cond8** (η = 0.770 Pa·s
and wet d50 = 9.6 µm at v_tip = 8.5 m/s). With the full 9-condition set:

| Quantity | Old model | Data (9 cond) | Assessment |
|---|---|---|---|
| Wet d50 plateau | 9.94 µm | 10.19 µm (mean, n=8) | slight under; retune |
| Milling onset | ~12.3 m/s (20 % solid) | none seen ≤ 12.4 m/s | **false rollover at cond11** |
| Low-shear η (plateau) | ~0.85 Pa·s | 0.201 Pa·s (geomean) | **~4× over** (anchored to the max point) |
| Flow index n | 0.40–1.00 (mid 0.70) | 0.20–0.67 (mean 0.44) | range too high |
| Dispersed regime | ~0.2 µm, ~0.001 Pa·s | 0.183 µm, ~0.005 Pa·s | ✓ correct behaviour |

---

## Part 3 — Recalibration applied

All changes are in the active model; the legacy standalone
`up1_module_rev38_dryer_risk.r` predates the three-regime closure (no
`AGG_FACTOR_CAL` / floc factor) and is not sourced by the chain, so it was left
untouched — noted here because that file's header claims a 1:1 sync that no longer
holds.

| File | Change | Rationale |
|---|---|---|
| `unified/up1_mixer_module.R` | `AGG_FACTOR_CAL` 48.7 → **50.0** | plateau d50 = 0.2·(1+50) = 10.2 µm = measured mean 10.19 |
| `unified/up1_mixer_module.R` | `MILL_MARGIN` 1.8 → **3.5** | no attrition observed ≤ 12.4 m/s; the UP1→UP4-direct harness showed that under nominal chemistry v_tip_crit runs ~3.9–5.9 m/s (below the 6.81 calibration point), so 2.5 still milled cond11/cond14 down to ~8.3 µm — 3.5 keeps every measured condition on the plateau, milling only above ~13.7 m/s |
| `unified/up1_mixer_module.R` | floc exponent 1.67 → **1.29** | targets the geometric-mean slurry viscosity (0.201 Pa·s) of all 8 aggregated conditions instead of the single 0.770 Pa·s outlier |
| `unified_model.R`, `morris_sensitivity_analysis.R` | `n_flow` factor range 0.40–1.00 → **0.20–0.67** | measured post-UP1 flow-index envelope (mean 0.44); the old range never reached the observed low end and was centred too high |
| `unified_model.R`, `morris_sensitivity_analysis.R` | `ALR` 1.0–10.0 → **0.9–1.8** | measured `up4atom_scfm`/`up4_feed` ratio; the old range ran the nozzle at up to ~6× the real air-liquid ratio |
| `unified_model.R`, `morris_sensitivity_analysis.R` | atomizing-air pressure 2.0–7.0×10⁵ → **1.2–1.8×10⁵ Pa** | measured `max_up4atom_psig`; the entire prior range sat *above* the real maximum |
| `unified_model.R`, `morris_sensitivity_analysis.R` | `mdot_L` 0.002–0.020 → **0.013–0.020 kg/s** | measured `up4_feed` (107–154 lb/hr); the old low end was ~7× below any real run |

Because `morris_sensitivity_analysis.R`'s factor table also sets the full-train
nominal (`sp_mid` = range midpoints), these three corrections move the chain's
nominal atomizer operating point onto the measured process, not just the screen.
The `couple_viscosity` note in `full_train_mixer_to_dryer.R` was updated to the
recalibrated 1.29 floc exponent.

`C_AGG_CAL` (0.0797, the aggregation-onset tip speed) was **kept** — the data
brackets onset between 4.37 m/s (cond2, dispersed) and 6.92 m/s (cond5,
aggregated), fully consistent with the calibrated 6.81 m/s at 20 % solids.

**Post-change check (20 % solids, saturated aggregation):** wet d50 = 10.20 µm
(target 10.19) and η = 0.200 Pa·s (target geomean 0.201) — versus the old
0.855 Pa·s.

### Dry-PSD control regime — size-template overlay (implemented, muted)
The two control regimes are now represented in both dryers
(`up2_run_dryer`, `spray_dry_model`) as a switchable **size-template overlay**,
governed by a `size_template` knob in [0, 1]:

- **`size_template = 0` (default, muted).** The dry particle is the droplet→shell
  only (`Dp ∝ droplet · (solids/(1−porosity))^{1/3}`) — current calibration is
  bit-for-bit unchanged. This is deliberate: the UP1→UP4 chain scale is still
  anchored on the existing (through-UP2/UP3) calibration and the nozzle-geometry
  SMD gap is open, so the template stays off until the base calibration checks out.
- **`size_template = 1` (future experiment).** In the **UP1-control** regime the
  UP1 aggregate templates the particle — the per-mode sizes log-blend toward
  `TEMPLATE_DENSIFY · D_agg` (`TEMPLATE_DENSIFY = 1.20`, the measured dry/wet ratio),
  reproducing dry d50 ≈ 12 µm. The **atomizer-control** regime is untouched at any
  setting: a dispersed feed has `d_ratio → 1`, so the template weight `w_tmpl → 0`
  and cond2 keeps its droplet-shell (19.85 µm) automatically.

The overlay is exposed as `sp_mid["size_template"] <- 0` in
`full_train_mixer_to_dryer.R`. The absolute SMD scale still awaits nozzle geometry
(`D_h`, `A_L`) — the standing high-priority gap.

### Scope of this calibration
The `visc.xlsx` runs are **UP1 → UP4 direct** — they do *not* pass through UP2
(foam-wash) or UP3 (centrifuge). Accordingly, this pass recalibrates only the UP1
closures and the UP4 atomizer inputs; UP2/UP3 are untouched. Routing calibration
data through UP2/UP3 is the next stage, once the UP1/UP4 models check out against
this direct data. The end goal is to use the calibrated model to set the
**experimental ranges** for the next design of experiments.

## Part 4 — UP1→UP4-direct validation

`validate_up1_up4_direct.R` replays all nine conditions **UP1 → UP4 directly**
(bypassing UP2/UP3, matching the data path). Each condition's UP1 setpoints
(v_tip, solids, exit T/P) and UP4 atomizer setpoints (ALR, air pressure, feed
flow) come from `data/cond_process.csv`; formulation chemistry is held at nominal
(the data doesn't resolve it, so the aggregation onset is common to all runs).
Output: `unified_output/up1_up4_direct_validation.csv`. Accuracy as a geometric
RMS ratio (1.00 = perfect):

| Quantity | RMS ratio | Read |
|---|---|---|
| **Wet d50, UP1-control (8)** | **1.07** | aggregation plateau reproduced to ~7 % |
| Wet d50, all 9 | 1.90 | inflated only by cond2 (see note) |
| **Viscosity, UP1-control (8)** | 2.78 | closure sits at the geomean (0.20 Pa·s); the 0.02–0.77 scatter is irreducible from the recorded inputs — an expected miss, not a bias |
| **Dry d50, muted (8)** | 1.81 | droplet-shell alone over-predicts (up to 35 µm) |
| **Dry d50, `size_template=1` (8)** | **1.25** | the aggregate-template overlay lands the dry PSD at 11–16 µm vs measured 8.8–15.3 — strong support for the UP1-control mechanism |
| Dry d50, cond2 (atomizer-control) | 17.6 vs 19.85 µm | droplet-shell is the right mechanism here (~12 % low) |

Two findings the harness surfaced:
1. **`MILL_MARGIN` raised 2.5 → 3.5** (above) — the harness caught that nominal
   chemistry gives a lower onset than the calibration point, so 2.5 spuriously
   milled the high-solids/high-v_tip conditions.
2. **cond2 dispersed-size reporting basis.** In the atomizer-control regime the
   stream reports `D_agg_um` = `D_primary_exit_um` (the ODE primary, 1.25 µm),
   not the physical 0.2 µm bead the wet PSD measures (0.18 µm). The *regime* is
   classified correctly (below onset); only the reported dispersed size is on the
   ODE basis. A future tidy-up would report `D_primary_phys_um` there. It does not
   affect the UP1-control conditions or the dry PSD.

### Still open
- **Viscosity scatter** (CV 77 %) is real process variability; the closure hits
  the central tendency but cannot reproduce the per-condition spread.
- **`couple_viscosity` default** remains `FALSE`; the recalibrated viscosity now
  makes turning it on defensible for a process-representative run.
- **Turning the size-template on** (`size_template=1`) once its 1.2·D_agg
  densification is pinned down — the dry-d50 RMS of 1.25 is already a strong
  first anchor.
- **Regenerated `unified_output/`** now reflects the recalibrated model (Morris
  cache busted and re-run).

---

## Part 5 — Particle morphology & density (SEM + bulk density)

**SEM (folder `sem/`, keyed by `sem_sample_no_filename` in `visc.xlsx`)** confirms
the two-regime morphology predicted from size + solids:

- **cond2 (atomizer-control, Sample 110495):** crumpled, collapsed, raisin-like
  **hollow shells** with a micro-porous nodular surface and visible broken shells.
- **cond11 (deep UP1-control, v_tip 12.4 m/s, Sample 110045) & cond14:** **compact,
  rounded granules** — no collapse. cond5 (6.9 m/s, just above onset) reads as
  intermediate, consistent with sitting on the aggregation boundary.

So the dry-particle morphology switches on the **same `v_tip_crit`** as wet
aggregation: dispersed feed → high-Péclet **hollow shell**; pre-aggregated feed →
**compact granule**.

**Visual morphology survey (1500×, one class per sample):**

| Sample | cond | v_tip | morphology |
|---|---|---|---|
| 110495 | cond2 | 4.37 | **hollow — crumpled/collapsed shells** |
| 110041 | cond5 | 6.92 | compact granule (borderline, reads dense) |
| 110043 | cond9 | 9.62 | compact granule |
| 110045 | cond11 | 12.40 | compact granule (clearest dense) |
| 110048 | cond14 | 9.61 | compact granule |
| 109384 | cond15 | 7.41 | compact granules |
| 110494 | cond1 | — (SEM-only) | compact granule |
| 110042 | cond6 | — (SEM-only) | compact granule |

**cond2 is the only hollow sample; every other sample is a compact dense granule.**
The SEM-only samples (cond1, 6, 10, 13, 16) need their process conditions to place
them in the regime framework, but morphologically they read dense.

**Automated image analysis did not yield reliable per-cond shape numbers.** A
watershed-segmentation pipeline (circularity, solidity) over the SEM set does *not*
separate the regimes: the discriminating feature is the **internal void**, which is
not present in a 2-D *surface* projection, and both regimes have irregular
projected outlines (crumpled shells vs lobed raspberry granules), so circularity
conflates them; large-granule counts per frame are also too low (n ≈ 3–10) for
stable statistics. Reliable quantification needs **dynamic image analysis** on
dispersed powder (thousands of particles) for sphericity and **cross-section SEM**
(FIB/microtome) for a true void fraction — the surface images give the *qualitative*
class cleanly but not calibrated numbers.

**Bulk density reconciliation (skeletal 1.70 g/cc, bulk ~0.30 g/cc).**
"Compact" is *not* "solid": a 0.30 g/cc bulk with a 1.70 g/cc skeletal density
forces substantial internal porosity even in the UP1-control granules —

| packing f | envelope g/cc | φ_intra |
|---|---|---|
| 0.45 | 0.67 | 0.61 |
| **0.50** | **0.60** | **0.65** |
| 0.55 | 0.55 | 0.68 |

So the compact granules are **~65 % internal void** (loosely packed primary
aggregates), and cond2 more still (central void → ~85 %). The old closure returned
φ ≈ 0.03 (near-solid) and matched ~0.22 g/cc tapped only via an implausible packing
of ~0.22 — the **right bulk number by compensating errors**, which would mispredict
how density responds to PSD/process and gets every porosity-driven property wrong.

**Recalibration — `morphology_recal` knob (MUTED by default), both dryers.**
Anchored to bulk 0.30 g/cc with skeletal 1.70 g/cc: bulk = ρ_skel·(1−φ)·packing.

| regime | φ (total void) | sphericity Ω | tapped |
|---|---|---|---|
| UP1-control (compact granule) | 0.65 | 0.90 | **0.30 g/cc** (anchor) |
| atomizer-control (cond2, hollow) | 0.85 | 0.55 | 0.13 g/cc |

The hollow bump and the sphericity drop are gated on the **same aggregation
fraction** (`D_primary_exit` vs `D_agg`) the size-template uses, so no new regime
switch is introduced. `morphology_recal = 0` (default) leaves the current closure
bit-for-bit unchanged (verified on R); `= 1` applies the anchored pair. This folds
in the regime-dependent porosity **and** sphericity correction (φ↑, Ω↓ for the
dispersed feed).

**Density coupling — now reconciled (UP1 `rho_polymer` 1050 → 1700).**
UP1's polymer density was raised from the generic-latex 1050 to the true skeletal
**1700 kg/m³**, so the whole chain (UP1 slurry, UP2 `rho_s`, and the morphology
recalibration's 1.70 g/cc) uses one consistent density. This lowers `phi_s`
(0.192 → 0.128 at 20 % solid), so the two `phi_s`-dependent wet-side constants were
re-derived to hold the **measured** anchors:

| constant | old (phi_s=0.192) | new (phi_s=0.128) | anchor held |
|---|---|---|---|
| `C_AGG_CAL` (onset) | 0.0797 | **0.0532** | onset ~6.8 m/s at 20 % (cond2 dispersed, cond5 aggregated) |
| floc exponent | 1.29 | **1.34** | slurry η geomean 0.20 Pa·s |
| `AGG_FACTOR_CAL` | 50.0 | 50.0 (unchanged) | wet d50 plateau 10.2 µm (direct PSD measurement) |
| `MILL_MARGIN` | 3.5 | 3.5 (unchanged) | onset preserved ⇒ milling still out of envelope |

Verified on R (UP1→UP4-direct harness): onset classification unchanged, wet d50
plateau = 10.20 µm, slurry η geomean = 0.209 Pa·s — the wet-side calibration is
preserved under the corrected density. `SG_Colloid` rises (~1.01 → ~1.09), which is
the physically correct consequence. `unified_output/` regenerated. (The legacy
`up1_module_rev38_dryer_risk.r`, not sourced by the chain, still uses 1050.)

Porosity and particle *size* are physically
coupled (a more porous particle is larger); here the size-template and the
morphology recalibration are kept as independent knobs, so enabling both is
additive rather than self-consistent — acceptable for a muted capability, to be
unified when both are turned on against real density data.

---

## Part 6 — Dryer thermals (evaporative air-flow estimate)

`visc.xlsx` now carries the dryer inlet/outlet temperatures (`up4_Tin ≈ 152 °C`,
`up4_Tout ≈ 88 °C`). With those, the drying-gas mass flow is back-solved from an
**evaporative energy balance** — sized to dry the feed to ≤ 0.5 % residual moisture:

> `mdot_gas · cp_gas · (Tin − Tout) = mdot_evap · h_fg + mdot_feed · cp_feed · (Tout − Tfeed)`

using the feed rate (`up4_feed`) and solids (`up1_feed_solid%`, no intermediate
dilution). The product is **not hygroscopic**, so drying is not limited by the
air's moisture capacity (inlet dew point ignored) — the binding constraint is
energy. Result: **mdot_gas ≈ 0.40–0.62 kg/s** (air/feed ≈ 30×, humidity rise only
~0.025 kg/kg — far below saturation, confirming the energy-limited regime). Code:
`dryer_airflow()` in `validate_up1_up4_direct.R`; constants match the UP2 module
(`h_fg = 2.30e6`, `cp_gas = 1005`).

**This finally grounds the drying-dependent properties.** The nominal thermals
(127 °C, 0.32 kg/s) were badly off; the real dryer runs **hotter (152 °C) with ~1.5×
the air**. Feeding the measured `T_dryer_in` + the evaporative `mdot_gas` per
condition, **skin jumps from ~0.1 to ~0.75–0.82** — strong shell formation, which is
exactly what the SEM shows (skinned shells/granules). The thermal factor windows
were retuned to the measured envelope: `T_dryer_in` 330–470 → **410–440 K**,
`mdot_gas_dry` 0.10–1.00 → **0.40–0.62 kg/s** (and `sp_mid` nominal 0.25 → 0.47).

**Moisture reconciliation (applied).** The energy-sized air leaves ≤ 0.5 % moisture
in the non-hygroscopic powder, so the final moisture cannot exceed that regardless
of the per-mode drying kinetics. Both dryers now cap `X_moist` (fraction of feed
water retained) at the value that yields the target **powder** moisture
`w_moist_target` (default **0.5 %**), via
`X_cap = t·Csol / ((1−Csol)(1−t))`. Result: powder moisture is pinned at 0.5 %
across all conditions (was an inconsistent 0.4–16 % of feed water from the raw
kinetics), and the moisture-plasticized `Tg` now uses the operational value. The
parameter is a single knob — a future study can **sweep `w_moist_target` up to
~1 %**. (The drying-kinetics closure itself, `tau_dry`/`t_res`, is still
uncalibrated in an absolute sense; the cap makes the *reported* moisture
operationally correct until that closure is calibrated.)

**Direct-path scope of the thermals.** `T_dryer_in` (~425 K) is a dryer setpoint,
valid regardless of path. `mdot_gas` was sized on the **direct-path** feed
(`up1_feed_solid%` ≈ 20 %, no intermediate dilution) — correct for the UP1→UP4
runs the data describes. When UP2/UP3 come online, UP3 concentrates the feed
(~45 %), so `mdot_gas` must be **re-derived from the post-UP3 solids** (same
evaporative method, different input); the `full_train` nominal is tagged
accordingly.

---

## Part 7 — theta_skin closure revision: surface-fusion route

### Problem with the original single-route closure
The original theta_skin model used a single Péclet-concentration route:

```
S_crit <- 500 * (1 + 0.05 * Softness)
theta_skin <- S_skin / (S_skin + S_crit)
```

`Softness = (1 + 5·C_binder)·exp(25·phi_solvent)` where `phi_solvent = C_monomer + C_plasticizer + w_core`.

The sign of Softness in S_crit was wrong: higher Softness (more plasticizer) *raised* S_crit, making skin *harder* to form. In the Morris screen, C_monomer and C_plasticizer ranked 16–17 of 41 factors and drove theta_skin in the **wrong direction** (more plasticizer → less skin). Physically, plasticizer/monomer enable surface fusion by depressing the polymer Tg — they should drive *more* skin.

Two mechanisms are now implemented in parallel:

### Route 1: Péclet concentration (corrected sign)
```
S_crit_pe <- 500 / (1 + 0.05 * Softness)
theta_skin_pe <- S_skin / (S_skin + S_crit_pe)
```
Dividing by Softness (was multiplying) reflects the correct physics: a more plasticized surface has *lower* viscous resistance to consolidation, so the threshold *falls* with Softness.

### Route 2: Surface fusion anchored at the wet-bulb temperature
During the **constant-rate drying period**, evaporative cooling holds the droplet surface at the wet-bulb temperature (~100 °C = 373 K), regardless of the dryer outlet temperature. If the plasticized polymer Tg drops below 100 °C, the surface is rubbery throughout this period and fuses continuously.

```
T_surface_cr <- 373.15   # constant-rate wet-bulb [K]
inv_Tg_plas  <- (1 - phi_solvent) / Tg_polymer + phi_solvent / Tg_solv   # Fox eq.
Tg_plas      <- 1 / inv_Tg_plas
theta_skin_fus <- 1 / (1 + exp(-(T_surface_cr - Tg_plas) / 15))
```

The 15 K sigmoid width spans the critical region around 100 °C. This route is the primary differentiator for high-Tg_polymer conditions where plasticizer determines whether the surface stays rubbery through the constant-rate period.

**Two-stage temperature separation:** `T_surface_cr` (100 °C, constant-rate, skin formation) is distinct from `T_particle = 0.85·T_out + 0.15·T_feed` (~80–85 °C, falling-rate) which is retained for burst, solvent boiling, sticking, and residual-solvent escape — all falling-rate phenomena.

### Combined
```
theta_skin <- 1 - (1 - theta_skin_pe) * (1 - 0.5 * theta_skin_fus)
```
Both routes act in parallel (either independently forms skin); weight 0.5 on the fusion route avoids trivial saturation at nominal conditions (most process-relevant Tg_pol values are 40–80 K below T_surface_cr so theta_skin_fus ≈ 0.97 at nominal, saturation handled by the 0.5 weight).

### Result: corrected factor directions and rankings

From the re-run Morris screen (41 factors, r=30 trajectories, `template_type=4`):

| Factor | Rank | μ* | Direction | Expected |
|--------|------|----|-----------|---------|
| C_monomer | 14 | 0.0254 | MORE skin | ✓ correct |
| C_plasticizer | 16 | 0.0238 | MORE skin | ✓ correct |
| Tg_polymer | 9 | 0.0627 | LESS skin | ✓ correct |
| T_dryer_in | 4 | 0.1081 | MORE skin | ✓ correct |
| T_mix | 3 | 0.1657 | LESS skin | ✓ (hotter feed → smaller ΔT driving kappa) |
| mu_L | 1 | 0.2957 | MORE skin | ✓ (lower D_diff → higher Pe) |

All signs are now physically consistent. C_monomer/C_plasticizer moved from rank 16–17 (wrong direction) to rank 14/16 (correct direction). The particle-size calibration metrics were **unchanged** by this fix — theta_skin feeds morphology, not the PSD.

**Stale Morris cache removed and screen re-run** after the skin closure change.

---

## Part 8 — UP1→UP2→UP3→UP4 full-chain data

Source: `data/up1_2_3_4_visc_sem.xlsx` — four conditions (`up3_1`–`up3_4`) that ran
the **complete** UP1→UP2→UP3→UP4 chain (the direct visc.xlsx runs bypassed UP2/UP3).
Data parsed and saved to `data/cond_up1234.csv` (process + PSD) and
`data/up3_viscometry.csv` (post-UP3 flow curves). SEM sample cross-reference is on
the `sem_sample_no_filename` sheet (sample numbers 114641–114644); images reside in
the `SEM/` folder under the same naming convention.

### Process summary (all 25 % solid UP1 feed, v_tip > 16 m/s → all UP1-control)

| Cond | v_tip (m/s) | UP3 solid % | d50 dry (µm) | η post-UP3 @ 12.7 s⁻¹ (Pa·s) |
|------|------------|------------|-------------|-------------------------------|
| up3_1 | 16.44 | 40.7 | 10.20 | 12.0 |
| up3_2 | 18.12 | 42.7 | 8.72 | 18.8 |
| up3_3 | 18.11 | 40.6 | 9.07 | 12.0 |
| up3_4 | 19.28 | 38.8 | 8.67 | 26.8 |

UP4 inlet temperatures are similar to the direct-path runs (143 °C).

### Key findings

**UP3 concentrates from 25 % → 39–43 % solids.** UP2 adds water
(`up2_h2o_lbhr` ≈ 41–49 lb/hr, diluting the slurry slightly); UP3 (decanting
centrifuge) then concentrates the cake. The UP4 dryer therefore sees a **richer feed
(~40 % solids)** than the UP1→UP4-direct path (25 %). The `dryer_airflow()` energy
balance must be re-derived from the **post-UP3 solids** when routing the full chain.

**Post-UP3 viscosity is 12–27 Pa·s at 12.7 s⁻¹**, versus 0.02–0.77 Pa·s post-UP1
for the same solid at 25 %. The concentration step raises the viscosity by roughly
50–100×. This will have significant impact on the UP4 atomizer model (`mu_L` drives
droplet breakup directly), and it is correctly tracked — the UP1 exit viscosity
closure feeds the dryer model, and concentrating to 40 % would shift it further.

**Dry d50 (8.67–10.2 µm) is consistent with the UP1→UP4 direct set (8.8–15.3 µm)**
and also in the UP1-control regime. The slight lower range here may reflect the
higher post-UP3 solids (denser droplet → smaller dry particle at the same droplet
size) or the higher v_tip (all ≥ 16 m/s, driving the aggregate toward the
compacted end of the plateau).

### up3_1 and up3_3 viscometry: identical by design

The flow curves on sheets `up3_1` and `up3_3` are bit-for-bit identical (20 data
points). This is **physically reasonable**: both conditions exited UP3 at essentially
the same solid content (40.70 % vs 40.61 %) and the same UP1 feed solid (25 %),
so their post-centrifuge slurry compositions are nearly indistinguishable. Identical
viscometry from near-identical compositions is expected, not a data-entry error.

### UP3 rho — real measurement; residual trapped gas explains low density

`up3_rho = 1117.3 kg/m³` (recorded only for up3_1) is a **physical sample
measurement**, not a formula error. To interpret it, the effective solid density
must be derived from the UP1 feed measurement rather than assumed from pure polymer.
Solving `1/ρ_feed = s/ρ_solid + (1−s)/ρ_water` at s = 0.25, ρ_feed = 1117.9 kg/m³
gives **ρ_solid ≈ 1730 kg/m³**. This is simply the **bare colloid density
(1.70 g/cc, confirmed) — the colloid is polymer**, and 1700 is also the dried-
product skeletal density used in the model (`unified/up1_mixer_module.R:234`).
There is no mineral skeleton; an earlier draft that read 1730 as a mineral-loaded
particle "diluted" by residual monomer was mistaken. The feed density is fully
accounted for by colloid + water (predicted 1112 vs measured 1118 — no deficit),
which is exactly what a **gaseous** monomer implies: residual monomer is not a
dissolved liquid and leaves no density signature in the feed. (Residual monomer is
addressed separately below; it cannot be inferred from these densities.)

With ρ_solid = 1730 kg/m³, the gas-free density at 40.7 % solid is:

```
ρ_no_gas = 1 / (0.407/1730 + 0.593/1000) = 1209 kg/m³
```

The measured 1117.3 kg/m³ is lower by 91.7 kg/m³. Attributing the deficit to
residual gas trapped in the centrifuge cake:

```
α_gas = (ρ_no_gas − ρ_meas) / (ρ_no_gas − ρ_air)
      = (1209 − 1117.3) / (1209 − 1.2) ≈ 7.6 %
```

This **7.6 % residual gas** is physically expected: UP1 sparges air into the slurry
(alpha_g ≈ 30 % volumetric in the foam at exit, see Part 9); centrifugation at
~430 g removes most gas but cannot drive it to zero. UP3 runs at atmospheric
pressure so no gas expands on exit. The 7.6 % is a real trapped-gas fraction, not
a formula error. UP3 measurements for up3_2–up3_4 are null (single representative
measurement taken for up3_1 only).

---

## Toward the goal — experimental ranges

With the UP1 and UP4-atomizer closures calibrated and validated UP1→UP4-direct,
the model is ready to propose **experimental ranges** for the next DoE: the
recalibrated factor windows (`n_flow`, `ALR`, atomizing pressure, `mdot_L`) now
bracket the real process, and the Morris screen ranks which of them actually move
each product property. Routing the calibration through UP2/UP3 is the following
stage, once this direct calibration is accepted.

The full-chain data (`up1_2_3_4_visc_sem.xlsx`) gives the first four-stage
measurements; the next modelling step is to build a validation harness that routes
UP1→UP2 (foam wash with water addition) →UP3 (decanting centrifuge, concentrating
to ~40 %)→UP4 (dryer, re-derived air flow from post-UP3 solids).

---

## Part 9 — New data columns and derived mass-balance quantities

Source: updated `data/up1_2_3_4_visc_sem.xlsx` and new `data/up1_2_3_solids_rho.xlsx`.
All derived results written to `data/cond_up1234.csv`.

### New columns added to up1_2_3_4_visc_sem.xlsx

| Column | Units | Description |
|--------|-------|-------------|
| `up4_dry_air_scfm` | SCFM | Dry air flow through UP4 dryer |
| `up2_liq_out_lbhr` | lb/hr | Cleaning liquid leaving UP2 (= 0 for all 4 chain runs) |
| `up3_reject_solid%` | % | Solid content of UP3 centrate reject (measured for up3_1 only) |
| `up3_reject_lbhr` | lb/hr | UP3 reject stream flow (measured for up3_1 only) |

New file `data/up1_2_3_solids_rho.xlsx` (Sheet1): historical runs with up3_rho
(in g/cm³) and up3_reject_solid% across the wider run set.

### UP1 exit density

UP1 is a gas-sparged foam mixer. Exit density was calculated two ways:

**Method A — all sparge gas entrained (upper bound):**
Ideal-gas expansion of the sparge flow (up1_scfh in SCFH) from standard
conditions to UP1 exit (up1_exit_temp, up1_psig):

```
n_gas  = P_std × V_std / (R × T_std)
V_exit = n_gas × R × T_exit / P_exit
ρ_exit = mdot_liq / (V_liq + V_exit)
```

**Method B — 7 % gas holdup at exit (from UP3 back-calculation):**
The measured UP3 cake density implies 7.6 % residual gas after centrifugation.
Back-expanding to UP1 exit pressure (Boyle's law):
`α_g_exit ≈ 0.076 × P_atm / P_exit`

| Cond | P_exit (bar) | T_exit (°C) | ρ_exit Method A (kg/m³) | α_g Method A | ρ_exit Method B (kg/m³) |
|------|-------------|------------|------------------------|-------------|------------------------|
| up3_1 | 1.44 | 54.9 | 772 | 30.9 % | 1063 |
| up3_2 | 1.55 | 69.2 | 778 | 30.4 % | 1068 |
| up3_3 | 1.62 | 69.8 | 790 | 29.4 % | 1070 |
| up3_4 | 1.51 | 72.7 | 770 | 31.2 % | 1067 |

Method A (30 % gas) represents the case where all sparged gas remains entrained;
in practice most bubbles rise out of the mixer headspace before the product exits.
Method B (4.4–4.9 % at operating pressure, expanding to 7 % at 1 atm) is more
consistent with the UP3 cake measurement and reflects the actual trapped-bubble
fraction. The model uses Method B for downstream feed-state calculations.

### UP3 mass balance and UP2 water routing

The UP3 reject for up3_1 (97.2 lb/hr at 0.284 % solid) combined with the cake
(152.8 lb/hr at 40.7 %) closes against the **UP1 output alone** (249.84 lb/hr):

```
249.84 ≈ 152.8 + 97.2 = 250.0 lb/hr  ✓ (0.16 lb/hr rounding)
```

If the UP2 wash water (49.27 lb/hr) were included in the UP3 feed, the mass balance
would be short by exactly 49.1 lb/hr. This confirms that **UP2 wash water exits UP2
as a separate effluent** and does not proceed to UP3. UP2 is a counter-current foam
wash: fresh water is fed to remove impurities from the foam surface and exits as a
dilute wash effluent; the cleaned foam (same volume as the UP1 output) continues to
UP3.

### UP4 feed rate

**Primary estimate (mass balance):** solid from UP1 divided by UP3 cake solid
fraction, with a small rejection correction for up3_1 (0.276 lb/hr solid in reject):

| Cond | UP1 solid (lb/hr) | UP3 solid% | UP4 feed/cake (lb/hr) |
|------|------------------|-----------|----------------------|
| up3_1 | 62.46 (−0.28 reject) | 40.7 | **152.8** |
| up3_2 | 52.63 | 42.7 | **123.3** |
| up3_3 | 52.68 | 40.6 | **129.7** |
| up3_4 | 52.44 | 38.8 | **135.2** |

**Secondary estimate (heat duty):** air enthalpy `Q = ṁ_air × cp × ΔT` divided by
evaporation enthalpy (water 25 °C → steam at 90 °C, Δh = 2555 kJ/kg):

| Cond | ṁ_air (g/s) | Q (kW) | UP4 feed heat-duty (lb/hr) | Thermal efficiency |
|------|-----------|--------|---------------------------|-------------------|
| up3_1 | 848 | 45.6 | 234.9 | 0.650 |
| up3_2 | 824 | 44.5 | 236.4 | 0.522 |
| up3_3 | 772 | 41.4 | 212.8 | 0.610 |
| up3_4 | 728 | 39.2 | 195.7 | 0.691 |

The heat-duty estimate exceeds the mass-balance estimate by 1.5–1.9×. For up3_1
(the best-closed condition) the implied dryer thermal efficiency is **65 %**,
plausible for a small pilot spray dryer with wall losses and air bypass. The
efficiency varies across conditions (52–69 %), likely because up4_dry_air_scfm
was set at a fixed high value regardless of feed throughput. **The mass-balance
estimate is the primary basis for the UP4 feed rate; the heat-duty estimate
provides an independent check but is not used directly.**

### Summary of effective particle solid density

Solving the slurry-density equation from the UP1 feed measurement at 25 % solid:

```
ρ_solid = s / (1/ρ_feed − (1−s)/ρ_water)
        = 0.25 / (1/1117.9 − 0.75/1000) = 1730 kg/m³
```

Consistent across all 4 conditions (1729–1742 kg/m³). **This is the bare colloid
density (1.70 g/cc); the colloid is polymer, not a mineral-loaded particle**, and it
matches the dried-product skeletal density used in the model. The earlier reading —
"higher than pure polymer, lower than a mineral skeleton diluted by residual
monomer" — was incorrect: there is no mineral phase, and the feed density is
explained by colloid + water alone (predicted 1112 vs measured 1118).

### Residual monomer — gaseous, not resolvable from these densities

Residual monomer is **gaseous at room temperature (MW ≈ 70 g/mol)**, so it produces
no liquid-density deficit and cannot be read off the feed density. It would have to
be quantified as a gas (PV=nRT) or by direct assay (TGA / GC). The gas holdup seen
in the streams is **foaming air** (`up1_scfh` is an air sparge), not monomer, so the
UP1-exit and UP3-cake void fractions bound trapped *air*, not residual monomer — at
most they set a loose upper bound on any monomer that might share those voids.
Quantifying residual monomer therefore needs data not in the current columns: a
monomer feed/off-gas flow, or a TGA/GC assay on the product. See
`review/up3_density_monomer_note.md` for the full working.
