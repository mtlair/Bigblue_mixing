# spray
spray

## Morris sensitivity analysis

`morris_sensitivity_analysis.R` implements the reduced-order two-fluid
atomization + drying model from the deep-research report (Feed Properties →
Two-Phase Conditioning → Nozzle Hydraulics → Primary Atomization → Secondary
Breakup → Drying / Particle Formation) and screens its inputs with the Morris
elementary-effects method for five outputs:

| Output | Nomenclature symbol |
|---|---|
| Final particle size | `D_particle` |
| Particle skin formation | `theta_skin,z` |
| Particle sphericity | `Omega_struct,z` |
| Particle porosity | `phi_porosity,z` |
| Powder tapped density | `rho_tapped` (bulk analogue of `rho_colloid,out` / `SG_out`) |

Run with:

```sh
Rscript morris_sensitivity_analysis.R
```

It uses the `sensitivity` package when installed and otherwise falls back to a
built-in base-R Morris OAT design (no dependencies required). Results are
written to `output/morris_sensitivity_plots.png` (μ* vs σ Morris plots, one
panel per output) and `output/morris_indices.csv`.
