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
bimodality are screened alongside the means. It screens 22 input factors
with the Morris elementary-effects method for twelve outputs:

| Output | Nomenclature symbol |
|---|---|
| Spray droplet size | `Dv50` |
| Distribution tails | `d10`, `d90`, `d99`, span, bimodality index |
| Final particle size | `D_particle` |
| Particle skin formation | `theta_skin,z` |
| Particle sphericity | `Omega_struct,z` |
| Particle porosity | `phi_porosity,z` |
| Powder tapped density | `rho_tapped` (bulk analogue of `rho_colloid,out` / `SG_out`) |
| Product glass transition | `Tg_eff` (Fox equation on residual solvent) |

Run with:

```sh
Rscript morris_sensitivity_analysis.R
```

It uses the `sensitivity` package when installed and otherwise falls back to a
built-in base-R Morris OAT design (no dependencies required). Results are
written to `output/morris_sensitivity_plots.png` (μ* vs σ Morris plots, one
panel per output) and `output/morris_indices.csv`.
