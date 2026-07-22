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

**Open coupling to flag:** the recalibration uses skeletal **1.70 g/cc for the dry
powder**, but the **UP1 slurry still uses `rho_polymer = 1050`**. If 1.70 is the
true polymer density, UP1 `phi_s` (and therefore the viscosity / aggregation-onset
calibration, which assumed 1050) should be revisited — a distinct pass, since it
shifts the wet-side anchors. Also, porosity and particle *size* are physically
coupled (a more porous particle is larger); here the size-template and the
morphology recalibration are kept as independent knobs, so enabling both is
additive rather than self-consistent — acceptable for a muted capability, to be
unified when both are turned on against real density data.

---

## Toward the goal — experimental ranges

With the UP1 and UP4-atomizer closures calibrated and validated UP1→UP4-direct,
the model is ready to propose **experimental ranges** for the next DoE: the
recalibrated factor windows (`n_flow`, `ALR`, atomizing pressure, `mdot_L`) now
bracket the real process, and the Morris screen ranks which of them actually move
each product property. Routing the calibration through UP2/UP3 is the following
stage, once this direct calibration is accepted.
