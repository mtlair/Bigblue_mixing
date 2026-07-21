# Handoff — Polymer process Morris sensitivity model

**Repo:** `mtlair/bigblue_mixing`
**Working branch:** `claude/polymer-centrifuge-review-teykp0` (all work below is pushed here)
**Date:** 2026-07-21

---

## 1. The flowsheet being modeled

```
UP1 mixer ─► pressurized foam-wash column ─► decanter centrifuge ─► reslurry ─► spray dryer
 (pulled)        (stub — you're building)      (evolved model)      (dilution)   (evolved model)
```

- **UP1 mixer** = "2nd step prior to the centrifuge" (gassed/templated mixing vessel, ODE-based).
- **Foam-wash column** = "1st step prior to the centrifuge" — a pressurized column that washes
  the foam out of UP1. **You are building this model.** In the chain it is currently a thin
  placeholder that only strips gas + sets discharge pressure/bubble size.
- Composition, gas holdup, surface chemistry, template state, T, P, and **floc strength** are
  defined **once** at the mixer and flow downstream through a shared `stream` interface.

---

## 2. File locations (all in repo root unless noted)

### Active model — the full train
| File | Role |
|---|---|
| `full_train_mixer_to_dryer.R` | **Top-level assembly.** Wires mixer → foam-wash → centrifuge → reslurry → dryer. Contains `foam_wash_column()` (placeholder), `stream_to_centrifuge()` adapter, `run_full_train()`. Run: `Rscript full_train_mixer_to_dryer.R` |
| `centrifuge_morris_sensitivity.R` | **Evolved decanter-centrifuge model.** `unified_centrifuge_model(run)` (27 outputs), `centrifuge_to_spray(run, target_solid_mass)` handoff, 25-factor Morris. Calibrated to plant (exit density 0.5–1.4 g/cc, cake ~40–57%, floc 2000 Pa). |
| `morris_sensitivity_analysis.R` | **Evolved spray-dryer model.** `spray_dry_model(x)`, 31 factors. Internal-mix bi-fluid nozzle, melt sticky-point (`T_sticky_K`), MW effect, surface-enrichment skin. **This is the dryer the train ends on** (not the pulled `up2`). |
| `chained_centrifuge_spray.R` | Earlier 2-unit chain (centrifuge → reslurry → dryer). Superseded by `full_train_*` but kept for the driver-decomposition / ΔT-vs-feed analyses. |

### Pulled upstream code — `unified/` directory + roots
Pulled from branch `claude/unified-model-two-modules-hvv4qd` (the follow-up to
`claude/r-script-missing-factors-jmlr8b`). This is a **parallel framework** with its own dryer;
the train reuses only its **mixer + stream interface**.

| File | Role |
|---|---|
| `unified/up1_mixer_module.R` | **UP1 mixer.** `up1_run_mixer(pars, equipment)`, `up1_default_equipment()`. deSolve ODE. |
| `unified/interface_stream.R` | **Shared stream contract.** `stream_from_up1()`, the two placeholder stages, `print_stream()`. |
| `unified/STAGE_CONTRACT.md` | **Field dictionary** for the stream — read this before building a stage. |
| `unified/up2_spray_dryer_module.R` | Their spray dryer (NOT used by the train; we end on the evolved one). |
| `unified_model.R` | Their standalone chain (mixer → 2 placeholders → their dryer) + own Morris. The train sources only its **definitions** (nominal_x, `up1_pars_from_x`, equipment) — not its run. |
| `up1_module_rev38_dryer_risk.r` | Original 915-line monolithic mixer (reference). |
| `morris_shim.R` | Offline base-R Morris fallback (CRAN `sensitivity` pkg is blocked in this sandbox). |

---

## 3. How the units connect (the interfaces)

- **mixer → foam-wash → centrifuge:** `stream_from_up1()` builds the stream; `foam_wash_column()`
  transforms it; `stream_to_centrifuge()` maps it onto the centrifuge input vector. It overwrites
  only the ~11 inputs the upstream physically sets (solids, gas holdup, additives, surface
  chemistry, T, P, bubble size) and leaves the centrifuge's operating knobs alone.
- **floc coupling:** mixer `Bond_Strength` drives centrifuge `floc_strength_Pa`, anchored so
  mixer-nominal Bond → the calibrated 2000 Pa. (The aggregates the centrifuge destroys are the
  ones the mixer built.)
- **centrifuge → reslurry → dryer:** `centrifuge_to_spray()` runs the centrifuge, dilutes the cake
  to a sprayable solids fraction, and returns the 12 stream-set spray inputs. Dryer operating
  setpoints stay independent.

---

## 4. Verified nominal run (capillary-bridge template)

| interface | key state |
|---|---|
| mixer | C_solid 0.23, gas holdup 0.095, Bond 0.101 |
| foam-wash | gas 0.095 → 0.024, D_b 46 µm, 3 atm |
| centrifuge | cake 43% solids, exit **1.23 g/cc**, floc 1271 Pa |
| reslurry | 30% solids, ρ 1154, dilution 1.76× |
| spray | Dp 19.9 µm, porosity 0.329, skin 0.798, **tapped 0.369 g/cc**, moisture 3.8% |

Both **exit density (1.23 g/cc, inside the calibrated 0.5–1.4 band)** and **tapped density
(0.369 g/cc ≈ your ~0.3 g/cc plant bulk density)** reproduce plant observations.

---

## 5. Environment / dependencies

- **`deSolve`** is required for the mixer ODE. It was installed in this session (with
  `gfortran`, `libblas-dev`, `liblapack-dev`) **but the container is ephemeral** — it will not
  persist to your next session. Reinstall, or run on a machine where `deSolve` is available.
- CRAN `sensitivity` is blocked in this sandbox → all Morris screens use a base-R OAT fallback
  (`morris_shim.R` / inline). No hard CRAN dependency beyond `deSolve`.

---

## 6. Open / next steps

1. **Foam-wash column model** — replace the `foam_wash_column()` stub in
   `full_train_mixer_to_dryer.R` (it already conforms to `f(stream, pars) -> stream`).
2. **Chained Morris screen over the full train** — not yet run. Before running, decide which of
   the ~24 mixer factors to sweep vs. fix, so the screen isn't ~60 factors at once. Locked config
   from prior discussion: fix `Tg_polymer=250K`, `T_sticky_K=493K`; outputs D_particle_um,
   sphericity, porosity, skin, tapped density; Morris r=20, group-colored panels.
3. **Note:** in the foam-wash gas-removal sweep the *final* powder porosity barely moves — the
   centrifuge already strips most gas and reslurry dilutes the rest, so the column matters more for
   centrifuge operability than final porosity. Confirm against plant experience.

## Branch lineage
- Source of mixer code: `claude/r-script-missing-factors-jmlr8b` → follow-up
  `claude/unified-model-two-modules-hvv4qd` (pulled in commit `a42beb1`).
- All chain assembly + evolved models: **`claude/polymer-centrifuge-review-teykp0`** (this branch).
