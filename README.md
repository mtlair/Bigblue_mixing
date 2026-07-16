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
`A_molecule`), DLVO electrostatics (`Delta_pH`, `I_strength`),
Krieger–Dougherty crowding from feed solids (`phi_s`), and impeller
dissipation (`v_tip`) conditioning bubble size and gas entrainment. It
screens 20 input factors with the Morris elementary-effects method for six
outputs:

| Output | Nomenclature symbol |
|---|---|
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
