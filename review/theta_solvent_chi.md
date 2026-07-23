# Template-from-chemistry predictor (Flory-Huggins χ / theta solvent)

**Module:** `theta_solvent_chi.R` · predicts a template's process behavior from its
chemical identity, no experiment required.

## Method
Flory-Huggins interaction parameter from Hildebrand solubility parameters:

```
χ = χ_S + (V_s / R T) · (δ_polymer − δ_solvent)²        χ_S ≈ 0.34
```

δ in MPa^0.5, V_s in cm³/mol (units reduce to dimensionless). **Theta point is
χ = 0.5**: below → good solvent (miscible, swells the matrix); above → poor solvent
(immiscible, phase-separates into discrete droplets). Combined with **volatility**
(boiling point vs drying temperature), this fixes the regime and the parameters the
process model already consumes.

## Regime map (χ + volatility → model)

| chemistry | regime | model template_type | feeds |
|-----------|--------|--------------------|-------|
| good + volatile | plasticize then **collapse** (avoid) | 3 (surface_weld) | w_core + collapse |
| good + non-volatile | permanent **plasticizer** | — (not a template) | φ_solvent / Softness |
| theta | marginal / mixed | 3 | — |
| poor + volatile | **clean pore template** (best) | 4 (capillary_bridge) | φ_templ_free → D_pore/porosity |
| poor + non-volatile | rigid **filler** / pore former | 1 (rigid) | rigid template load |
| gas | **PV=nRT** expansion | 2 (gas) | α_g / chi_npgas |

`RTF` (core-absorbed fraction) is predicted from miscibility —
`RTF = 1/(1+exp((χ−0.5)/0.1))` — high for good solvents (absorb into cores →
plasticize), ~0 for poor solvents (stay as free droplets → pores). This is exactly
the split `interface_stream.R` uses between `w_core` and `phi_templ_free`.

## Worked screen (polymer δ = 19.0, T_dry = 90 °C — illustrative)

| template | χ | solvency | regime | type |
|----------|----|----------|--------|------|
| water | 5.3 | very poor | clean pore template | 4 |
| toluene / MMA / acetone | 0.36 | good | plasticize→collapse (avoid) | 3 |
| styrene / DBP | 0.34 | good | permanent plasticizer | — |
| n-hexane / ethanol | 1.1–1.4 | very poor | clean pore template | 4 |
| cyclohexane | 0.51 | theta | marginal | 3 |
| paraffin oil | 1.5 | very poor + non-volatile | rigid filler | 1 |
| liquid CO₂ | — | gas | PV=nRT | 2 |

The pattern is the actionable result: **a good solvent plasticizes and collapses
(bad pore fidelity); a poor/non-solvent phase-separates and templates clean pores;
a non-volatile poor solvent stays as a rigid filler; a gas expands by the ideal-gas
law.** Theta (χ≈0.5) is the boundary between plasticizing and pore-templating.

## What this gives you now vs what still needs data
- **Now (from tabulated chemistry, no experiment):** direction and ranking of any
  candidate template — plasticize / collapse / clean-pore / rigid / gas — plus the
  `template_type`, `RTF`, `T_bp_solv` and density to drop into the model.
- **Still needs calibration:** absolute pore size / porosity / burst fraction. The
  model's template→pore closures (`D_pore`, `phi_templ`) carry fitted constants not
  yet tied to a measured *templated* product; 2–3 templated SEM samples fix the scale.
- **One input to pin:** `δ_polymer` (19.0 is a placeholder for a moderately polar
  vinyl/acrylic matrix). Set it from the real resin (group-contribution or a swelling
  test); the δ values for named solvents are literature Hildebrand values.

## Next integration step (optional)
Wire `template_from_chemistry()` into the chain so `template_type`, `RTF` and
`T_bp_solv` are *derived* from the chosen template instead of set by hand — turning
the four hard-coded template_types into a continuous, chemistry-driven input.
