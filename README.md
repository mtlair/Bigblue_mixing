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
coarsening (Ostwald ripening + coalescence over `t_hold`, pressure-
accelerated, surfactant-retarded), two-stage atomization (a Lund-type
effervescent stage — annular film, Rayleigh ligament breakup, flash-expansion
shattering — feeding a traditional Rizk–Lefebvre bi-fluid airblast stage),
and a co-current dryer energy/moisture balance (`mdot_gas_dry`, `Y_in` →
outlet temperature and humidity driving the drying Péclet number and the
sticky-point state). The v47 impeller closure (`v_tip`) was removed as
belonging to an upstream unit operation. It screens 21 input factors with
the Morris elementary-effects method for seven outputs:

| Output | Nomenclature symbol |
|---|---|
| Spray droplet size | `Dv50` |
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
