# Handoff — Full Mixer→Dryer Train + Morris Screen

**Branch (all work here):** `claude/mixer-workflow-summary-setup-vigdp5`
**Repo:** `mtlair/Bigblue_mixing`   **Remote HEAD:** `bda4b77`
**Runtime:** R 4.3.3 + `deSolve` 1.40 (installed). Everything else base R.

---

## What is built and working

A four-stage process train wired through one shared `stream` object, plus a
full-train Morris sensitivity screen:

```
UP1 mixer ──▶ UP2 foam-wash (491-line ODE) ──▶ UP3 separator ──▶ reslurry ──▶ UP4 dryer
 up1_mixer_module   foam_wash_column_psd     centrifuge_...      (dilution)   morris_sens...
```

- `Rscript full_train_mixer_to_dryer.R` — nominal end-to-end run + a wash sweep.
- `Rscript full_train_morris.R` — **54-factor Morris screen** over the 5 particle
  features. r=12 → 660 runs, ~66 s, **2.7% NA**. `R_TRAJ=20` for a finer screen.

**Screened features (after UP4):** porosity, sphericity, size, skin, tapped density.

**Merge rules that make the screen valid:**
1. **UP1 output bounds downstream** — every downstream input the stream sets is
   dropped as a factor (not double-counted). Includes the dryer's `T_feed`,
   `phi_emulsion`, `D_template`, now bound from the mixer stream.
2. **Equipment/geometry fixed** — UP2 column + UP3 bowl geometry held at nominal.
3. **UP2 enters via coalescence knobs only** (`K_coal/K_break/K_burst/d_b_burst/
   frac_gas_coarse_ref`).
4. **NA-guarded** — each eval is `tryCatch` + finiteness-checked → all-NA row, not
   an abort. Solver chatter suppressed. Widest ranges narrowed to keep NA low.

Current merged dictionary: **54 factors** (UP1 24 + UP2 5 + UP3 11 + UP4 14).
Ranges are in `full_train_output/full_train_factor_ranges.csv`.

---

## Key findings this thread

- **Drivers of the 5 features:** UP4 thermals dominate (`T_system`, `mdot_L`,
  `mdot_gas_dry`, `T_bp_solv`); after binding the template path upstream, **UP1
  `template_dose`, `C_solid_mass`, `T_mix`** rose into the top drivers; **`ALR`**
  (widened to 0.01–2) became the #1 driver of particle size.
- **Low tip-speed puzzle (resolved):** `v_tip` 2→22 m/s strongly moves the *mixer*
  state (WetSkin 0.003→0.093, alpha_g 0.115→0.053, born bubble 110→2.9 mm) but the
  5 final features stay flat. The outputs `v_tip` moves — gas holdup, wet-skin seed,
  bubble size — are exactly the ones the downstream **buffers neutralize**: UP2
  washes/coalesces gas, UP3 separates gas and **resets solids to a fixed 30%**, and
  the dryer is insensitive to feed bubble size. This is the strongest argument for
  the pending redesign (below): untie the fixed reslurry and tie bubble/film to
  local viscosity + surfactant so the mixer's hydrodynamics reach UP4.
- **`T_bp_solv`** = template-solvent boiling point [K, 300–360]. Sets the dryer
  thermal window via three gates: nozzle flash (`S_flashn`, T_feed vs bp),
  boil-vs-retain (`S_bs`, T_particle vs bp), and Clausius–Clapeyron burst
  overpressure (`dP_vap`). Below bp → solvent retained (plasticizer/VOC); near bp →
  clean pore templating; well above → over-pressurization/micro-explosion.
- **Monomer / plasticizer / binder are independent formulation inputs** (no
  correlation among them in the model). Their coupling is via shared downstream
  targets in UP4:
  - `phi_solvent = C_mono + C_plas + w_core` → `Softness = (1+5·C_bind)·exp(25·phi_solvent)`
  - `Perm_shell = exp(k_pm·C_mono + k_pp·(C_plas+w_core)) / (1 + k_pb·C_bind)`
    (monomer/plasticizer open the shell; binder plugs pores)
  - `mu_serum = mu_L·(1+10·C_bind)` (binder thickens serum)
  - Fox `w_res` residual-solvent → depresses product Tg (monomer/plasticizer).

---

## Design work — status

**DONE (this session):**
- **UP2 wash-liquid feed** — `foam_wash_column_ode()` now takes `wash_ratio`
  (wash liquid / foam, drives the column's internal impurity wash) and
  `wash_carry` (overhead fraction leaving WITH the foam into UP3). The carried
  wash dilutes the stream solids/additives handed to UP3 (`dil_w = 1 +
  wash_carry·wash_ratio`).
- **UP2 film/bubble tied to local viscosity + surfactant** — the column
  continuous-phase (serum) viscosity is now `1e-3·(1+10·C_binder)` with Andrade
  scaling to the column T (feeds settling, drainage, bubble burst); film
  elasticity / born-bubble size already keyed to the stream surfactant. Nominal
  σ shifted 36→32 mN/m and bubble size grew accordingly.
- **UP3→UP4 untie** — `centrifuge_to_spray()` now takes `reslurry_add` (added
  liquid per unit cake, default 0 = spray the cake as-is) instead of a fixed 0.30
  target. The dryer feed solids track the separator cake (`C_solid = Cs_cake/(1+
  reslurry_add)`), and the entrained gas is carried from UP3. Effect: **UP3
  `S_base` (dewatering floor) and `delta_rpm` rose into the top drivers** — UP3
  output now bounds UP4. Screen still 54 factors, 2.7% NA.

**DONE (mechanism E — colloid milling/destabilization):**
- **v_tip effect on primary particle size** — implemented colloid milling coupling:
  - Higher v_tip (mixer shear) reduces primary particle size `a_prim` (linear scaling)
  - Higher v_tip increases destabilization factor for `d_ratio` 
  - Final particle size reduced via improved packing (smaller primaries → smaller voids)
  - Morris ranking: v_tip now #4 driver for particle size (mu* = 57.4, was ~1 before)
  - Reaches top-3 for porosity, top-7 for sphericity, top-5 for skin
  - ALR still dominant but gap narrowed from 356x to 5.5x

**STILL PENDING:**
- **Mechanism A (airblast microbubble)** — couple v_tip-set microbubble size to
  bi-fluid nozzle SMD. Microbubbles (D_b_m from UP2) already flow to dryer via
  centrifuge handoff; dryer already uses them in g_in term. Could strengthen coupling
  by modulating SMD_ab directly or feed-side viscosity feedback.
- **Two-pump pressure structure** — pump 1 (UP1+UP2+UP3, one back-pressure) and
  pump 2 (UP4). Currently pressures are per-stage/stream-carried; not yet
  collapsed to two pump levels.
- **UP3 wash as a distinct sprayed feed** — today `reslurry_add` is a lumped
  liquid add. The spec wants a UP3 wash sprayed on the cake, part to the light
  phase, layered on top of the g-force clearing (which already exists). 
- **Refine v_tip scaling exponents** — current exponents (1.0 for a_prim scaling,
  1.0 for mill_size_factor) may be adjustable. Diagnostic sweeps on v_tip at fixed
  thermal/composition could reveal optimal coupling strength.

Original spec (for reference):

- **UP2 and UP3 each get their own wash-liquid feed input.**
- **UP2 wash:** co-current / counter-current flow in the column; the wash liquid
  **leaves with the foam into the UP3 inlet** (it becomes part of the stream to UP3).
- **UP3 wash:** sprayed, disperses into the cake, **part exits with the light phase**.
  The **wash itself washes little** (sprayed on a viscous foam/slurry) — model it as a
  minor dilution/entrainment term. But **do not conflate the weak wash with the unit's
  clearing power**: UP3 concentrates solids in the liquid feed and its **impurity
  removal is g-force-driven (centrifugal), which can be large**. The unit as a whole is
  a strong impurity remover even though the wash spray is not.
  - Already in the model (`centrifuge_morris_sensitivity.R`): separation scales with
    g-force — `omega = rpm·2π/60`, centrifugal settling `v_c ∝ omega²`, recovery
    `1 − exp(−|v_c|·t_clarify/dist)`, `g_eff = r_pool·omega²`. Higher rpm → faster
    settling → sharper clarification cut (`Cut_Size_d50c_um`), more solids into the
    cake (`Product_Solids_MassFrac`), more fines/dissolved rejected to the centrate
    (`Yield_Loss_Centrate`, `Centrate_Conc_kg_m3`). So `rpm` is the impurity-clearing
    lever; the new wash feed rides on top as a small dilution/light-phase term.
- **Post-separator feed to UP4 must be tied to UP3 output** (solids + liquid + some
  entrained gas in foam) — **replace the fixed 30% reslurry** target with the
  separator's actual output.
- **Temperature** is set by the mixer upstream and carried — **not** an independent
  per-stage input.
- **Pressure** is system back-pressure from **two pumps**: **pump 1 pushes UP1+UP2+UP3**
  (one shared pressure level), **pump 2 pushes UP4** (a second level). So pressure is
  two pump levels, not per-stage-independent.
- **Bubble size and film elasticity must tie directly to local viscosity and
  surfactant properties.** (Much of this already exists in the 491-line column ODE —
  `derive_state_props`: Szyszkowski σ / Gibbs elasticity from Langmuir surfactant,
  Hinze born-bubble size, drainage τ from μ(T). The wrapper currently ties σ/CMC to
  the stream surfactant but keeps column viscosity at default — wire it to the stream.)

### The two scope decisions (were open; resolved by the spec above)

1. **What is an independent (screened) input vs. inherited/tied.**
   Resolution from the spec: **temperature is inherited** from the mixer (not a
   factor); **pressure collapses to two pump levels** (pump-1 for UP1–3, pump-2 for
   UP4) rather than per-stage; the **two wash-liquid feeds are new inputs**; and the
   **UP3→UP4 feed is tied to UP3 output** (drop the fixed 30% reslurry). Remaining
   sub-choice: screen the two wash feeds + two pump pressures as Morris factors, or
   hold them as fixed setpoints.

2. **How much physics for the UP3 wash and the UP3→UP4 tie.**
   Resolution from the spec: **UP2 wash = real co/counter-current transport** that
   joins the foam into UP3; **UP3 wash = minimal** (light dilution, part to light
   phase) — but **UP3 as a unit is a strong impurity remover via g-force**, which the
   model already captures (rpm → omega² → settling → clarification cut + solids
   concentration). So the wash is a small term layered on top of the existing
   centrifugal clearing, not the impurity mechanism. The **reslurry is replaced** by a
   UP3-output-tied feed carrying solids + liquid + entrained gas.

---

## File directory (branch `claude/mixer-workflow-summary-setup-vigdp5`)

| File | Role |
|---|---|
| `full_train_mixer_to_dryer.R` | 4-stage train driver; `run_full_train()`; UP2 `up2="ode"`/`"algebraic"`; dryer T_feed/phi_emulsion/D_template bound from stream |
| `full_train_morris.R` | **Full-train Morris screen** (54 factors, NA-guarded, 5 features) |
| `unified/up1_mixer_module.R` | UP1 mixer physics (deSolve ODE) |
| `unified/interface_stream.R` | `stream` object + `stream_from_up1()` + placeholder intermediate stages |
| `unified/up2_spray_dryer_module.R` | UP2-as-dryer module (used by `unified_model.R`, not the full train) |
| `unified/STAGE_CONTRACT.md` | Stream field dictionary + stage contract |
| `foam_wash_column_psd.R` | **UP2 foam-wash 491-line ODE** (bimodal bubbles, Langmuir/Gibbs film, T/P, conserved gas). Plots muted for non-interactive runs |
| `foam_wash_module.R` | UP2 algebraic closure (`up2="algebraic"` path) |
| `centrifuge_morris_sensitivity.R` | UP3 separator/centrifuge model (`unified_centrifuge_model`) |
| `morris_sensitivity_analysis.R` | UP4 spray-dryer model (`spray_dry_model`) |
| `unified_model.R` | 2-module unified chain (UP1→dryer) + its own Morris; source of the mixer factor dict |
| `morris_shim.R` | Offline `sensitivity::morris` fallback |
| `up1_module_rev38_dryer_risk.r` | Original standalone UP1 mixer (source of truth) |
| `GAS_TEMPLATE_ACCOUNTING.md` | Open design note: three-stream gas/template split |
| `WORKSPACE_SUMMARY.md` | Earlier workspace summary |
| `HANDOFF.md` | This file |
| `full_train_output/` | `full_train_morris_indices.csv`, `full_train_na_by_factor.csv`, `full_train_factor_ranges.csv`, `full_train_morris_features.png` (`.rds` cache gitignored) |
| `unified_output/` | `unified_model.R` outputs |

**Not in this branch** (elsewhere): `foam-wash-column-review-pgfauf` and
`column-coalescence-bubble-burst-inrt41` hold earlier foam-wash variants;
`polymer-centrifuge-review-teykp0` and `r-script-missing-factors-jmlr8b` hold
centrifuge / mixer-rev history.

## Branches

- **`claude/mixer-workflow-summary-setup-vigdp5`** — active, all current work.
- `main` — base.
- Feature history: `claude/up1-mixer-workflow-model-11tcoc`,
  `claude/unified-model-two-modules-hvv4qd`, `claude/polymer-centrifuge-review-teykp0`,
  `claude/column-coalescence-bubble-burst-inrt41`, `claude/foam-wash-column-review-pgfauf`,
  `claude/r-script-missing-factors-jmlr8b`, and others (see `git branch -r`).
