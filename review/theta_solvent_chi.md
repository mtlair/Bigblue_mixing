# Template-from-chemistry predictor (Hansen HSP / theta solvent)

**Module:** `theta_solvent_chi.R` · predicts a template's process behavior from its
chemical identity — no experiment required.

## Method: Hansen Solubility Sphere

Relative Energy Distance from the polymer's Hansen sphere:

```
Ra² = 4(ΔδD)² + (ΔδP)² + (ΔδH)²        RED = Ra / r₀
```

δ in MPa^0.5. **Theta point is RED = 1**: below → good solvent (miscible,
swells/plasticizes); above → poor solvent (immiscible, phase-separates into
discrete droplets).

### Polymer Hansen parameters (pinned from characterization data)

| δD | δP | δH | r₀ | Hildebrand equiv. |
|----|----|----|-----|-------------------|
| 17.00 | 12.10 | 10.20 MPa^0.5 | 8.0 | 23.2 MPa^0.5 |

r₀ from supplier/cloud-point data; user data indicates r₀ closer to 8 than
14. Refine with a swelling or cloud-point titration if tighter RED boundaries
are needed.

### Why Hansen, not single Hildebrand

For this polymer (high δP + δH), only solvents with matching polarity **and**
H-bonding components are "good" solvents. Single Hildebrand collapses the
three components and misclassifies polar species with different δP/δH profiles:

| template | RED (Hansen) | χ (Hildebrand) | actual behaviour |
|----------|-------------|----------------|-----------------|
| DBP plasticizer | 0.91 (theta/good) | 1.16 (poor ✗) | known plasticizer ✓ |
| MMA | 0.95 (theta) | 1.30 (poor ✗) | near-theta in practice ✓ |
| toluene | 1.70 (very poor) | 1.25 (poor) | consistent — poor solvent ✓ |
| butyl acetate | 1.20 (poor) | 1.19 (poor) | consistent ✓ |

## Water-dispersibility filter

In aqueous spray-drying the template must form **stable droplets** in the
slurry, not dissolve into the continuous water phase. A candidate is viable only if:

```
water_sol < 5 g/100mL   OR   logP > 1.5
```

Acetone, ethanol, acetonitrile etc. fail this test and dissolve into the
aqueous phase regardless of their RED with the polymer.

## Regime map (RED + volatility + water_dispersible → model)

| solvency (RED) | volatile? | w-dispersible? | regime | template_type | feeds |
|----------------|-----------|---------------|--------|---------------|-------|
| good (RED<0.9) + volatile | yes | yes | plasticize → **collapse** (avoid) | 3 (surface_weld) | w_core + collapse |
| good/theta + non-volatile | no | yes | permanent **plasticizer** | — (not a template) | φ_solvent / Softness |
| theta (0.9–1.1) | — | — | marginal / mixed | 3 | — |
| poor (RED>1.1) + volatile | yes | yes | **clean pore template** (best) | 4 (capillary_bridge) | φ_templ_free → D_pore/porosity |
| poor + non-volatile | no | yes | rigid **filler** / pore-former | 1 (rigid) | rigid template load |
| any | — | NO | **not viable** — dissolves in slurry | — | — |
| gas (bp < 25 °C) | — | — | **PV=nRT** expansion | 2 (gas) | α_g / chi_npgas |

**Volatility threshold:** `bp_C < T_dry_C + 40`. For spray-drying, set
`T_dry_C` to the **inlet** temperature (≈140–160 °C) rather than outlet
(≈90 °C) when assessing whether a template escapes in-process.

`RTF` (core-absorbed fraction) is a continuous sigmoid on RED:
`RTF = 1/(1 + exp((RED − 1)/0.1))` — high for good solvents (absorb into
cores → plasticize), ≈ 0 for poor solvents (stay as free droplets → pores).

## Candidate screen (polymer: δD=17.00, δP=12.10, δH=10.20, r₀=8, T_dry=90 °C outlet)

### Viable pore-template candidates (type 4)

| template | RED | bp (°C) | logP | RTF | note |
|----------|-----|---------|------|-----|------|
| butyl acetate | 1.20 | 126 | 1.78 | 0.12 | mild poor solvent, low water sol ✓ |
| chloroform | 1.27 | 61 | 1.97 | 0.08 | fast evaporation; **Class 2 solvent — check ICH Q3C** |
| toluene | 1.70 | 111 | 2.73 | 0.01 | clean pore template, water insol ✓ |
| n-hexane | 2.05 | 69 | 3.90 | ~0 | very poor, highly volatile ✓ |
| cyclohexane | 1.96 | 81 | 3.44 | ~0 | very poor, volatile ✓ |
| styrene | 1.63 | 145 | 2.95 | ~0 | pore template but residual monomer concern |
| d-limonene* | 1.49 | 176 | 4.57 | ~0 | type **1** at T_dry=90 C; type **4** at T_dry=140 C (inlet) |

*d-limonene is borderline volatile: bp 176 °C means it evaporates readily at
spray-dryer inlet conditions (≈150 °C) but is classified as rigid filler at
the 90 °C outlet reference. Use `template_from_chemistry("d-limonene", T_dry_C=140)`
to evaluate at inlet conditions — result flips to type 4 (clean pore template),
making it an attractive green-solvent candidate.

`methylene_chloride` (DCM, RED=0.94) is theta/good → collapse path (type 3),
not a pore template; also Class 2.

### Permanent plasticizers (type NA)

| template | RED | logP | RTF | water_sol (g/100mL) |
|----------|-----|------|-----|---------------------|
| DBP | 0.91 (theta) | 4.72 | 0.50 | <0.001 |

DBP sits right at the theta boundary (RED=0.91 < r₀). RTF≈0.5 means roughly
half absorbs into particle cores and half stays as free phase. In practice DBP
is a well-known vinyl plasticizer — consistent with Hansen RED < 1.

### Rigid fillers / pore-formers (type 1, non-volatile)

| template | RED | logP | note |
|----------|-----|------|------|
| paraffin oil | 2.00 | 7.0 | stays as solid inclusions → porosity after leaching |
| isopropyl myristate | 1.51 | 6.0 | food-grade, non-volatile rigid filler |
| d-limonene (at T_dry=90 C) | 1.49 | 4.57 | reclassifies to type 4 at T_dry=140 C |

### Not viable (dissolve into aqueous phase)

acetone, ethanol, triethyl citrate (TEC), water itself.
MMA (water_sol=1.5 g/100mL) is borderline dispersible (passes the < 5 g/100mL
threshold) but its RED=0.95 puts it in the theta/collapse path at 90 °C — not
a useful pore template.

## What this gives you now vs what still needs data

- **Now (from tabulated chemistry, no experiment):** direction and ranking of
  any candidate template — plasticize / collapse / clean-pore / rigid / gas —
  plus `template_type`, `RTF`, `T_bp_solv` and density to drop into the model.
- **Recommended screen candidates:** butyl acetate (mild, safe) and toluene
  (strong poor-solvent, water-insol) for type-4 clean-pore; d-limonene at higher
  inlet temperature (green, food-contact potential); paraffin or IPM for
  rigid pore-formers.
- **Still needs calibration:** absolute pore size / porosity / burst fraction.
  The model's template→pore closures (`D_pore`, `phi_templ`) carry fitted
  constants not yet tied to a measured *templated* product; 2–3 templated SEM
  samples fix the scale.
- **r₀ to pin:** refine from a cloud-point titration or equilibrium swelling test
  in a good/poor solvent pair. Current r₀=8 gives physically consistent results
  (DBP correctly at theta, toluene at very poor) but a measured value removes
  the remaining uncertainty.

## Next integration step (optional)

Wire `template_from_chemistry()` into the chain so `template_type`, `RTF` and
`T_bp_solv` are *derived* from the chosen template instead of set by hand —
turning the four hard-coded template_types into a continuous, chemistry-driven
input. Pass `poly = POLY_HSP` and `T_dry_C` = dryer inlet temperature.
