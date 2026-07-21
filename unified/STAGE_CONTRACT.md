# Stream & Stage Contract

Single reference for building the intermediate unit-operation stages (and any
new upstream/downstream module) that plug into the unified chain. If you are
starting a fresh thread to build a stage, **read this file first**, then the
source module you are wrapping. You should not need the rest of the chain's
history.

```
UP1 (mixer) ──► UP2 (foam-wash) ──► UP3 (centrifuge) ──► reslurry ──► UP4 (spray dryer)
```

- Chain wiring:            `unified_model.R` → `run_unified()`
- Stream definition/adapter: `unified/interface_stream.R`
- Upstream module:         `unified/up1_mixer_module.R`  (`up1_run_mixer`)
- Downstream module:       `unified/up2_spray_dryer_module.R` → now UP4 in the full train

---

## The stream object

A stream is a plain named `list` describing the slurry **leaving one stage and
entering the next**. Every stage is a function

```r
stage(stream, pars = list()) -> stream
```

that returns a stream with the **same field set**, overwriting only the fields
the unit operation physically changes. A stage must not add or drop fields —
downstream code indexes them by name. If a stage needs a new physical quantity,
add the field to `stream_from_up1()` (initialised to a sensible default) so the
schema stays uniform across the whole chain, and note it here.

### Field dictionary

Units are SI unless noted. "Set by" = the stage that currently establishes the
value; "Consumed by" = where it is read downstream.

#### Composition (mass fractions of total slurry unless noted)
| Field | Meaning | Set by | Consumed by |
|---|---|---|---|
| `C_solid` | latex/polymer solids fraction (diluted by template feed in-tank) | UP1 | UP4 `C_sol` (+ `C_solid_rigid`) |
| `C_solid_rigid` | rigid-template solids fraction (0 for liquid/gas templates) | UP1 | UP4 `C_sol` |
| `C_binder` | binder concentration | UP1 | UP4 rheology, permeability, softness |
| `C_monomer` | residual monomer | UP1 | UP4 softness, permeability, Fox Tg |
| `C_plasticizer` | plasticizer | UP1 | UP4 softness, permeability, Fox Tg |
| `C_surfactant` | surfactant | UP1 | UP4 surfactant stoichiometry |
| `CMC` | critical micelle concentration | UP1 | (carried; not yet read by UP4) |
| `HLB` | surfactant HLB | UP1 | UP4 (currently unused knob) |
| `MW_surfactant` | surfactant MW **[g/mol]** | UP1 | UP4 (converted to kg/mol) |
| `A_molecule` | surfactant molecule area **[nm²]** | UP1 | UP4 (converted to m²) |
| `Delta_pH` | ΔpH vs isoelectric point | UP1 | UP4 DLVO `E_rep` |
| `ionic_strength` | ionic strength [M] | UP1 | UP4 DLVO screening |

#### Physical state
| Field | Meaning | Set by | Consumed by |
|---|---|---|---|
| `rho_slurry` | slurry density [kg/m³] | UP1 | UP4 `rho_L` |
| `rho_polymer` | dry polymer density [kg/m³] (UP1 uses 1050) | UP1 | UP4 `rho_s` |
| `T_K` | stream temperature [K] | UP1 exit T | UP4 `T_feed`; **pre-heater (between UP3 and UP4) should overwrite** |
| `P_Pa` | stream pressure [Pa] | UP1 headspace | **pump (between UP1 and UP2) should overwrite** |
| `mu_exit_PaS` | mixer-exit apparent viscosity [Pa·s] (diagnostic) | UP1 | (diagnostic only) |

#### Particulate state
| Field | Meaning | Set by | Consumed by |
|---|---|---|---|
| `D_particle_um` | primary colloid particle diameter [µm] | UP1 | UP4 `a_prim`, `d_ratio` |
| `D_agg_um` | blended aggregate size at mixer exit [µm] | UP1 | UP4 `d_ratio` (`D_particle/D_agg`) |
| `sphericity` | Ω at mixer exit (diagnostic) | UP1 | (diagnostic only) |
| `WetSkin` | wet-skin fraction at mixer exit | UP1 | UP4 skin seed (× `SKIN_SURVIVAL`) |

#### Gas phase
| Field | Meaning | Set by | Consumed by |
|---|---|---|---|
| `alpha_g` | entrained gas holdup (trapped microbubbles) | UP1 `Blended_Porosity` (capped 0.75) | UP4 `alpha_g_0` (via UP3); **UP2 may remove** |
| `D_b_m` | bubble diameter [m]; **`NA` until a stage sets it** | — | UP4: uses this if finite, else falls back to the `D_b` screening factor. **UP2 should set it.** |

#### Template phase
| Field | Meaning | Set by | Consumed by |
|---|---|---|---|
| `template_type` | 1=rigid, 2=gas, 3=surface_weld, 4=capillary_bridge | UP1 equipment | UP4 gating |
| `phi_templ_free` | free (interstitial) liquid-template volume fraction | UP1 (RTF split) | UP4 `phi_emulsion` (pore templating) |
| `D_template_um` | template droplet/seed diameter [µm] | UP1 | UP4 `D_e`; **stages may ripen it** |
| `RTF` | Residual_Template_Fraction: liquid template absorbed into cores | UP1 rev38 | UP4 (via `w_core`) |
| `w_core` | core-absorbed template, mass fraction of slurry | UP1 (RTF × load) | UP4 softness, permeability, Fox Tg, burst inventory |

#### Diagnostics (carried, not consumed by the dryer closure)
`Softness_exit`, `Bond_Strength`, `Retained_Porosity`, `Mixing_Potential`.

---

## Interface constants (in `interface_stream.R`)

These are **assumptions**, not derived — revisit them when a real stage lands.

| Constant | Value | Meaning |
|---|---|---|
| `RHO_TEMPLATE_LIQ` | 750 kg/m³ | template solvent liquid density (matches UP2 `rho_solv`) |
| `SKIN_SURVIVAL` | 0.30 | fraction of mixer wet-skin surviving atomization (new surface is created at the nozzle) |
| `UP1_EPS_RCP` | 0.36 | random-close-packing void fraction; used to convert `template_dose` (fraction of **void**) into a slurry volume fraction |

The RTF split (in `stream_from_up1()`):
```
phi_templ_liq   = template_dose · phi_s · EPS_RCP/(1-EPS_RCP)   (capped 0.30)
phi_templ_free  = phi_templ_liq · (1 - RTF)     # free droplets → dryer pores
w_core          = phi_templ_liq · RTF · RHO_TEMPLATE_LIQ/rho_slurry   # core plasticizer
```

---

## Stage contracts

Each stage is `f(stream, pars = list()) -> stream`. Both are **identity
functions today** — they return the stream unchanged. Fill in the physics and
overwrite only the listed fields.

### `intermediate_stage_1` — transfer line / feed pump (suggested)
Physically: the slurry is pumped from the mixer to the nozzle feed line under
shear. Expected to set/overwrite:
- `P_Pa` — pump discharge pressure (currently carried straight from mixer headspace; the dryer's own `P_feed` factor stands in for now).
- `D_b_m` — bubble size after line shear. **This is the field the dryer most wants** (it is `NA` on the mixer exit, so the dryer currently falls back to a `D_b` screening factor). Setting it here retires that factor.
- optionally `D_template_um` — shear break-up of the emulsion.

### `intermediate_stage_2` — hold tank / pre-heater / degasser (suggested)
Physically: residence before the nozzle at possibly elevated temperature.
Expected to set/overwrite:
- `T_K` — pre-heat (feeds the dryer flash/superheat calcs and the Fox Tg state).
- `alpha_g` — gas venting/degassing over the hold.
- `D_b_m`, `D_template_um` — Ostwald ripening / coalescence over the hold time (UP2 already has an in-nozzle `t_hold` ripening model in Module 0d — decide whether hold-tank ripening belongs here or stays there, and don't double-count).

> Note on overlap: UP2 Module 0d/2 already models a pressurized hold + bubble
> coarsening *at the nozzle* (`t_hold`, `k_rip0`, `k_emu0`). If stage 2 becomes
> a real hold tank, move that residence physics upstream into the stage and
> reduce UP2's `t_hold` to the short nozzle-line residence only — otherwise the
> hold is counted twice.

---

## How to add a stage (checklist)

1. Branch from (or work on) `claude/unified-model-two-modules-hvv4qd` so
   `unified/` is present.
2. Read this file + the source physics you are wrapping.
3. Implement the stage body in `unified/interface_stream.R` (or a new
   `unified/stageN_*.R` sourced by `unified_model.R`), overwriting only the
   fields it changes; keep the field set intact.
4. If you add a stream field, add it to `stream_from_up1()` with a default and
   document it in the dictionary above.
5. Run `Rscript unified_model.R`; confirm the nominal chain still solves for all
   four templating strategies and the Morris screen completes.
6. If a stage now supplies a value the dryer took as a screening factor (e.g.
   `D_b`), drop that factor from the `factors` table in `unified_model.R`.
