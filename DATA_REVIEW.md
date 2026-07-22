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

**Aggregation is a sharp on/off switch in tip speed.**

| cond | v_tip (m/s) | solid % | wet d50 (µm) | regime |
|---|---|---|---|---|
| 2 | 4.37 | 13.2 | **0.183** | dispersed (primary bead) |
| 5 | 6.92 | 19.6 | 9.97 | aggregated |
| 15 | 7.41 | 24.9 | 9.71 | aggregated |
| 8 | 8.48 | 19.2 | 9.60 | aggregated |
| 14 | 9.61 | 25.0 | 9.74 | aggregated |
| 9 | 9.62 | 19.6 | 10.37 | aggregated |
| 12 | 9.62 | 19.0 | 9.96 | aggregated |
| 4 | 10.57 | 20.0 | 10.19 | aggregated |
| 11 | 12.40 | 19.6 | 12.00 | aggregated |

The one condition below ~6.8 m/s (cond2) stays at the ~0.18 µm primary particle;
every condition above it jumps to a **flat ~10 µm aggregate plateau**
(mean 10.19 µm, CV 7 %) that barely moves with tip speed. There is **no milling
rollover** anywhere in the tested envelope — cond11, the highest tip speed
(12.4 m/s), has the *largest* wet d50, not a reduced one.

**Viscosity is highly scattered and not a clean function of the process.**
Across the eight aggregated conditions (all ~19–25 % solids, all ~10 µm aggregate)
the low-shear viscosity spans 0.020–0.770 Pa·s — a **38× range, CV 77 %**. It does
not track solids, tip speed, or wet d50 (|corr| ≤ 0.35 for each). cond8 (0.770) is
the high outlier and cond12 (0.020, run at the high 150 lb/hr feed) the low one.
This is genuine process/sample variability; a deterministic closure can only be
expected to hit the **central tendency** (geometric mean 0.201 Pa·s, median 0.227).

**Dry PSD is set by atomization, not by the wet aggregate.**
Post-UP4 dry d50 is 8.8–15.3 µm for the aggregated feeds (dry/wet ratio ≈ 0.7–1.6,
mean ~1.2) and 19.85 µm for the dispersed cond2 — a 108× jump over its 0.18 µm wet
size, because each atomized droplet dries to a particle built from many colloidal
beads. Backing out droplet size from dry d50 and the solids loading gives a
physically sensible airblast SMD band of **15–27 µm**. Within the aggregated set
the dry d50 is *anti*-correlated with UP1 tip speed (corr −0.91) and wet d50
(−0.79): harder mixing upstream yields a finer dry powder. No single measured
atomizer variable (`max_up4atom_psig`, `up4atom_scfm`, `up4_feed`) explains it on
its own, so the dry PSD is genuinely multi-factor.

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
| `unified/up1_mixer_module.R` | `MILL_MARGIN` 1.8 → **2.5** | no attrition observed ≤ 12.4 m/s; pushes v_tip_mill to ~17 m/s, out of the validated range, removing the spurious cond11 rollover |
| `unified/up1_mixer_module.R` | floc exponent 1.67 → **1.29** | targets the geometric-mean slurry viscosity (0.201 Pa·s) of all 8 aggregated conditions instead of the single 0.770 Pa·s outlier |
| `unified_model.R`, `morris_sensitivity_analysis.R` | `n_flow` factor range 0.40–1.00 → **0.20–0.67** | measured post-UP1 flow-index envelope (mean 0.44); the old range never reached the observed low end and was centred too high |

`C_AGG_CAL` (0.0797, the aggregation-onset tip speed) was **kept** — the data
brackets onset between 4.37 m/s (cond2, dispersed) and 6.92 m/s (cond5,
aggregated), fully consistent with the calibrated 6.81 m/s at 20 % solids.

**Post-change check (20 % solids, saturated aggregation):** wet d50 = 10.20 µm
(target 10.19) and η = 0.200 Pa·s (target geomean 0.201) — versus the old
0.855 Pa·s.

### Not changed — flagged for follow-up
- **Dry-PSD (UP4) recalibration** was not attempted: the target band is well
  characterised (dry/wet ≈ 1.2, implied droplet SMD 15–27 µm, anti-correlated
  with UP1 shear), but it needs each condition's atomizer settings mapped to the
  dryer's `ALR` / `P_atom_air` inputs before `Dp50` can be tuned against it.
- **Viscosity scatter** (CV 77 %) is real process variability; the closure now
  hits the central tendency but cannot reproduce the per-condition spread from
  the recorded factors alone.
- **`couple_viscosity` default** remains `FALSE`; the recalibrated viscosity now
  makes turning it on defensible for a process-representative run.
