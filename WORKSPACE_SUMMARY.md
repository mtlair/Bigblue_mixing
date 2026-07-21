# Mixer → Dryer Workflow — Workspace Summary

**Branch:** `claude/mixer-workflow-summary-setup-vigdp5`
**Purpose of this pass:** consolidate the unified mixer→dryer model files into one
runnable workspace, install R, verify the code runs end-to-end without errors,
and record what is done vs. still open.

Source material pulled together here:
- Unified chain + modules — from `claude/up1-mixer-workflow-model-11tcoc`
  (`unified_model.R`, `unified/`, `full_train_mixer_to_dryer.R`, helpers).
- `GAS_TEMPLATE_ACCOUNTING.md` — from `claude/mixer-workflow-summary-0kvecp`.

---

## How to run (workspace is runnable directly)

R 4.3.3 + `deSolve` 1.40 are installed. From the repo root:

```sh
Rscript unified_model.R              # UP1 mixer -> 2 placeholder stages -> UP2 dryer
Rscript full_train_mixer_to_dryer.R  # UP1 -> foam-wash -> separator -> reslurry -> dryer
```

`unified_model.R` writes to `unified_output/`:
- `nominal_chain_summary.txt` / `nominal_chain_outputs.csv` — nominal end-to-end
  run for all four templating strategies (rigid / gas / surface_weld / capillary_bridge).
- `unified_morris_indices.csv` — Morris μ*/σ indices, 41 factors × 35 outputs.
- `unified_morris_{process,surface,polymer}_variables.png` — μ*–σ lens plots.

Only hard dependency is `deSolve` (UP1 ODE integration); everything else is base R.
The Morris design is a base-R OAT implementation, so the CRAN `sensitivity`
package is not required (`morris_shim.R` is an offline fallback if it were).

---

## Error check — unified model: **clean, no errors**

| Check | Result |
|---|---|
| Parse all four files (`unified_model.R`, `unified/*.R`) | OK |
| Module sourcing + name tables (`up1_output_names`=13, `up2_output_names`=19) | OK |
| Nominal run, all 4 templating strategies | all outputs finite |
| Full Morris screen (41 factors, 30 trajectories → 1260 chain runs) | **1260/1260 valid**, ~36 s |
| `full_train_mixer_to_dryer.R` (4-stage train + wash sensitivity) | runs, EXIT 0 |

The only `NA` that appears is `D_b_m` in the verbose stream print — this is
**intentional**. Per the stream contract, `D_b_m` (feed bubble diameter) is `NA`
at mixer exit until a transfer-line/pump stage sets it; the dryer detects the
`NA` and falls back to the `D_b` screening factor. No numeric output is `NA`.

Cross-file symbol resolution verified: `UP1_EPS_RCP` (used in
`interface_stream.R`) and `SKIN_SURVIVAL` (used in `up2_spray_dryer_module.R`)
are defined in `up1_mixer_module.R` / `interface_stream.R` respectively, both
sourced before use in `unified_model.R`.

---

## Completed work

### Unified chain architecture (`unified_model.R` + `unified/`)
- **Single stream interface** (`unified/interface_stream.R`): a named-list
  `stream` object carries composition, physical state, particulate/gas/template
  state and structure history from the mixer exit into the dryer feed. Factors
  the two standalone screens used to duplicate (solids, additives, surface
  chemistry, template) are now defined **once, at the mixer**; the dryer sees
  their transformed values (dilution, gas holdup, wet skin, residual template).
- **UP1 mixer** (`unified/up1_mixer_module.R`): rev38/v52 two-compartment
  (cavern A / dead-wall B) ODE model, extracted verbatim from
  `up1_module_rev38_dryer_risk.r`. C-∞ smoothed algebra, surface chemistry
  (CMC cap, HLB, DLVO), templating closure (rigid/gas/surface-weld/capillary),
  rev38 collapse-risk outputs (`Swelling_Softness_exit`, `Residual_Template_Fraction`).
- **UP2 spray dryer** (`unified/up2_spray_dryer_module.R`): effervescent–airblast
  atomization + particle formation, refactored to take the `stream` as feed
  instead of independent knobs. Shell-permeability closure and immiscible
  template-solvent emulsion (pore templating, balloon inflation, micro-explosion)
  are wired in.
- **Cross-module couplings:** mixer wet-skin seeds the dryer skin state
  (`SKIN_SURVIVAL`); core-absorbed template (`w_core`, from UP1 RTF) feeds the
  dryer softness / permeability / Fox-Tg / burst inventory.
- **Unified factor dictionary:** 41 factors (24 mixer-side + 17 dryer-side),
  name collisions resolved (`P_mix`/`T_mix`/`tau_mix` vs `P_atom_air`/`T_dryer_in`).
- **Contract doc** (`unified/STAGE_CONTRACT.md`): the stream field dictionary,
  interface constants, stage contracts and the checklist for adding a stage.

### Full four-stage train (`full_train_mixer_to_dryer.R`)
- UP1 mixer → UP2 pressurized foam-wash column → UP3 separator/centrifuge →
  reslurry → UP4 nozzle dryer, wired through the same stream contract.
- Includes a wash-efficiency sensitivity sweep (more foam washed out → less gas
  into UP3/UP4 → denser, less porous powder).
- `foam_wash_column()` is a **deliberately thin placeholder** — it removes/
  compresses entrained gas and sets discharge pressure; everything else is
  pass-through, ready to swap for the real column model.

---

## Open items

### 1. Two intermediate stages are identity placeholders
In `unified_model.R`, `intermediate_stage_1` (transfer line / feed pump) and
`intermediate_stage_2` (hold tank / pre-heater / degasser) return the stream
unchanged. Per `STAGE_CONTRACT.md`:
- **Stage 1** should set `P_Pa` (pump discharge) and, most importantly, `D_b_m`
  (line-shear bubble size) — setting it retires the dryer's `D_b` screening factor.
- **Stage 2** should update `T_K` (pre-heat), `alpha_g` (venting) and ripen
  `D_b_m` / `D_template_um`. **Watch the overlap:** UP2 Module 0d already models a
  pressurized hold + bubble coarsening at the nozzle (`t_hold`, `k_rip0`,
  `k_emu0`) — if stage 2 becomes a real hold tank, move that residence physics
  upstream and shrink `t_hold`, or the hold is double-counted.

### 2. Gas & template three-stream split (`GAS_TEMPLATE_ACCOUNTING.md`) — not implemented
Template gas and foam gas are currently folded into one `alpha_g` holdup. They
must be tracked as three distinct populations:
- **Trapped template gas** — *inside* the solid particles, sized by `D_template`;
  survives to the dryer as intra-particle porosity (needs a dedicated field, a
  gas analogue of `w_core`/`RTF`, kept out of `alpha_g`).
- **Interstitial foam holdup** — *between* particles; this is what `alpha_g`
  represents and what the foam-wash column removes.
- **Slug escape** — large slugs venting overhead; not carried downstream.

Requirement: close the gas mass balance `G_in = G_trapped + G_holdup + G_slug`,
ensure UP2 washes **only** the interstitial stream (never the trapped fraction),
and route the trapped-template stream to `D_template`-sized intra-particle
porosity in UP4.

### 3. Real foam-wash column model
`foam_wash_column()` is a stub. The real pressurized foam-wash column model is
being built separately; drop it in by replacing the one function (it already
conforms to `f(stream, pars) -> stream`).

### 4. Interface assumptions to revisit when real stages land
`RHO_TEMPLATE_LIQ = 750`, `SKIN_SURVIVAL = 0.30`, `UP1_EPS_RCP = 0.36` are
assumptions, not derived (see `STAGE_CONTRACT.md`).

### 5. Emulsion constants still fixed (candidate future factors)
`rho_solv` (750), `h_fg_solv` (350 kJ/kg), `k_emu0` (2e-20), and a `D_pore_min`
capillary-collapse threshold are hard-coded in UP2 constants.

---

## File map

| File | Role |
|---|---|
| `unified_model.R` | Unified chain driver: nominal runs + Morris screen + plots |
| `unified/up1_mixer_module.R` | UP1 mixer physics (deSolve ODE) |
| `unified/up2_spray_dryer_module.R` | UP2 spray dryer + particle formation |
| `unified/interface_stream.R` | Stream object + placeholder intermediate stages |
| `unified/STAGE_CONTRACT.md` | Stream/stage contract reference |
| `full_train_mixer_to_dryer.R` | 4-stage train (adds foam-wash + separator) |
| `foam_wash_module.R` | Foam-wash column (placeholder) + build notes |
| `centrifuge_morris_sensitivity.R` | UP3 separator/centrifuge model |
| `morris_sensitivity_analysis.R` | Standalone UP4 dryer screen (`spray_dry_model`) |
| `up1_module_rev38_dryer_risk.r` | Original standalone UP1 mixer (source of truth) |
| `morris_shim.R` | Offline `sensitivity::morris` fallback |
| `GAS_TEMPLATE_ACCOUNTING.md` | Open design item: three-stream gas/template split |
| `FOAM_WASH_BUILD.md` | Foam-wash column build notes |
