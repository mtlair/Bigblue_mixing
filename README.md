# spray
spray

## Morris sensitivity analysis

`morris_sensitivity_analysis.R` implements the reduced-order two-fluid
atomization + drying model from the deep-research report (Feed Properties →
Two-Phase Conditioning → Nozzle Hydraulics → Primary Atomization → Secondary
Breakup → Drying / Particle Formation), extended with closures from the Latex
Coagulation Engine Integrated Master Specification (v30.0.0 / R-Engine v47):
Flory–Huggins free-volume swelling by residual solvent (`C_monomer`,
`C_plasticizer`, `C_binder`), a Fox-equation product glass transition with
stickiness / pore-collapse / caking above `Tg_eff`, surfactant molar
stoichiometry (`theta_surf` from `C_surfactant`, `MW_surfactant`,
`A_molecule`), DLVO electrostatics (`Delta_pH`, `I_strength`), and
Krieger–Dougherty crowding from feed solids (`phi_s`).

Spray-line extensions beyond the reactor spec: pressurized-hold bubble
coarsening (Ostwald ripening + coalescence over `t_hold` at the liquid-line
pressure `P_feed`, surfactant-retarded), an effervescent–airblast hybrid
exit plane (nothing atomizes in the confined feed line; the entrained gas
acts at the bi-fluid exit, where its throat voidage thins the annular liquid
film feeding the Rizk–Lefebvre correlation, and bubbles smaller than the
film flash-shatter the fragments on the final letdown), an explicit nozzle
pressure ladder (`P_feed` → chamber `0.8·P_system` → near-choked throat →
ambient), and a co-current dryer energy/moisture balance (`mdot_gas_dry`,
`Y_in` → outlet temperature and humidity driving the drying Péclet number
and the sticky-point state). The v47 impeller closure (`v_tip`) was removed
as belonging to an upstream unit operation. The droplet population is a
three-mode volume mixture (main atomized mode, starved-air coarse tail,
shear-strip/bubble-debris fine mode), so distribution quantiles and
bimodality are screened alongside the means, and the particle-formation
module runs per mode (gas retention, porosity, drying time), giving a
particle size distribution and a tapped density that feel the tails.

Rheology is power-law shear-thinning (`mu_app = mu·(γ̇/100 s⁻¹)^(n_flow−1)`,
Krieger–Dougherty crowded): atomization sees the apparent viscosity at the
nozzle film shear rate (~10⁵–10⁷ s⁻¹) while droplet-scale drying transport
uses the low-shear value. Cake mechanics compete the Rumpf-type yield
strength of the consolidated shell (weakened by Flory–Huggins softness and
above-Tg mobility, strengthened by coagulation bonding) against meniscus
capillary pressure to decide pore collapse vs lock-in, and the under-dried
coarse tail carries residual moisture that plasticizes the product Tg
(Fox, three components) and degrades cake strength and packing.

It screens 23 input factors with the Morris elementary-effects method for
fifteen outputs:

| Output | Nomenclature symbol |
|---|---|
| Spray droplet size | `Dv50` |
| Distribution tails | `d10`, `d90`, `d99`, span, bimodality index |
| Final particle size | `Dp50`, `Dp90` (per-mode particle distribution) |
| Particle skin formation | `theta_skin,z` |
| Particle sphericity | `Omega_struct,z` |
| Particle porosity | `phi_porosity,z` (mass-weighted across modes) |
| Powder tapped density | `rho_tapped` (bulk analogue of `rho_colloid,out` / `SG_out`) |
| Product glass transition | `Tg_eff` (Fox: polymer + solvent + moisture) |
| Residual moisture | fraction of feed water retained (wet coarse tail) |
| Cake yield strength | `sigma_y` [MPa] vs capillary collapse (`Pi_col`) |

Run with:

```sh
Rscript morris_sensitivity_analysis.R
```

It uses the `sensitivity` package when installed and otherwise falls back to a
built-in base-R Morris OAT design (no dependencies required). Results are
written to `output/morris_sensitivity_plots.png` (μ* vs σ Morris plots, one
panel per output) and `output/morris_indices.csv`.

## Unified process model (UP1 mixer → spray dryer)

`unified_model.R` connects the two standalone modules end-to-end through a
common feed-stream interface:

```
UP1 gassed/templated mixer  ──►  intermediate stage 1 (placeholder)
                                  └►  intermediate stage 2 (placeholder)
                                        └►  UP2 spray dryer / particle formation
```

- `unified/up1_mixer_module.R` — the rev38 mixer physics
  (`up1_module_rev38_dryer_risk.r`) packaged as a pure function
  `up1_run_mixer()`, with no study driver.
- `unified/up2_spray_dryer_module.R` — the spray-drying physics from
  `morris_sensitivity_analysis.R` refactored as `up2_run_dryer(feed, x)`,
  where `feed` is the stream handed over by the upstream chain.
- `unified/interface_stream.R` — the stream definition (composition,
  physical state, particulate/gas/template state, structure history), the
  UP1→stream adapter, and the two identity placeholder stages reserved for
  the unit operations to be inserted later (e.g. transfer line/pump, hold
  tank/pre-conditioner).

Factors the two standalone screens used to duplicate (solids, additives,
surface chemistry, template) are now defined once at the mixer, and the
dryer receives their transformed values through the stream: diluted solids,
slurry density, exit temperature, trapped-gas holdup → `alpha_g_0`, wet-skin
fraction → dryer skin seed, primary/aggregate size → `a_prim`/`d_ratio`,
surfactant MW/area/HLB → the dryer's surfactant stoichiometry, and the
liquid template split by the rev38 `Residual_Template_Fraction` into free
interstitial droplets (pore templating) vs core-absorbed plasticizer
(`w_core`: enters Flory–Huggins softness, shell permeability, Fox residual
solvent, and the micro-explosion volatile inventory — the rev38 collapse-risk
pathway wired into the dryer).

Run with:

```sh
Rscript unified_model.R
```

(needs `deSolve`; everything else is base R). It writes to `unified_output/`:
a nominal end-to-end run for each templating strategy
(`nominal_chain_summary.txt`, `nominal_chain_outputs.csv`) and a Morris
elementary-effects screen over the full 41-factor chain
(`unified_morris_indices.csv`, zone-classified
`unified_morris_{process,surface,polymer}_variables.png`; circles = mixer-side
factors, triangles = dryer-side factors).
