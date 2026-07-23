# SEM morphology metrics — definitions, extraction, and model calibration

**Date:** 2026-07-23 · Samples 114641–114644 = up3_1..up3_4 (dried product).
Extractor: `sem_morphology.py` → `data/sem_morphology.csv`.

## Metric definitions

Each SEM metric is dimensionless and maps to one model morphology output
(`unified/up2_spray_dryer_module.R`).

| Metric | Source mag | Definition | Model output |
|--------|-----------|------------|--------------|
| **Ω_SEM** (sphericity) | 500× | median particle **solidity** = area / convex-hull area, over watershed-segmented particles (aspect = minor/major and circularity 4πA/P² reported alongside) | `Omega_struct_z` |
| **φ_SEM** (surface porosity) | 5k/10k× | dark interstitial **pore-pixel area / particle-body area**, on the curvature-flattened surface | `phi_porosity_z` |
| **θ_SEM** (skin/fusion) | 5k/10k× | **smooth-area fraction** = fraction of the surface whose local texture CV < τ (0.06); a fused skin is smooth, an unfused bead packing is granular. `rms_rough` (mean texture CV) and `bump_density` (primary-bead peaks) reported alongside | `theta_skin_z` |

## Model vs SEM (4 chain conditions)

| cond | Ω model / SEM | φ model / SEM | θskin model / SEM(fused) |
|------|---------------|---------------|--------------------------|
| up3_1 | 0.969 / 0.881 | 0.083 / 0.171 | 0.934 / 0.035 |
| up3_2 | 0.973 / 0.883 | 0.080 / 0.154 | 0.924 / 0.283 |
| up3_3 | 0.971 / 0.879 | 0.079 / 0.172 | 0.920 / 0.127 |
| up3_4 | 0.976 / 0.882 | 0.078 / 0.153 | 0.909 / 0.143 |

## What the images show (calibration implications)

1. **The product surface is an unfused granular packing of primary beads**
   (~100–200 nm), with open interstitial pores — not a smooth melt skin. This is
   the central finding and it splits the model's `theta_skin`:
   - The model's `theta_skin_z ≈ 0.92` is a **consolidation / Péclet-enrichment**
     measure — a packed shell IS present — and should be read that way, *not* as
     surface fusion.
   - The **fusion** the images actually show is low: θ_SEM ≈ 0.04–0.28 (mean
     ~0.15). This is the quantity the rev38 **surface-fusion route**
     (`theta_skin_fus`) predicts, and its low nominal value is consistent with
     the SEM. **Recommend calibrating `theta_skin_fus` against θ_SEM** and keeping
     the consolidation and fusion sub-terms reported separately.

2. **Sphericity:** model `Omega_struct_z ≈ 0.97` vs SEM solidity ≈ 0.88 (aspect
   ≈ 0.72). The model runs ~0.09 high; a fixed offset/gain on `Omega_struct_z`
   reconciles it. Both are flat across the four conditions, as the model predicts.

3. **Surface porosity:** SEM surface pore-area (≈0.16) is ~2× the model
   `phi_porosity_z` (≈0.08). These are different quantities — SEM sees *surface*
   interstitial porosity, the model reports *structural/bulk* porosity — so the
   map is a calibration factor (~2×), not a 1:1 comparison. A fracture-section
   SEM would be needed to measure true internal porosity.

4. **Cross-condition signal for the DoE:** up3_2 is distinctly smoother / more
   fused (θ_SEM 0.28, rms_rough 0.12) than up3_1 (θ_SEM 0.04, rms_rough 0.21).
   Surface fusion therefore *does* vary with process conditions even though d50
   and the other metrics are flat — a lead worth resolving once the SEM metrics
   are wired as calibration targets.

## Surface porosity — % open space between spheres (`sem_porosity.py`)

Per the reviewer's definition (114641_5kx annotation): **skin** = smooth, gap-free
patches where beads have merged (darker, low local-range); **open pore** = the
interstitial gap between distinct packing spheres (dark, bounded by bright caps,
high local range). Pixels inside the particle body are classified by local
intensity **range** (max−min over a bead-scale window) plus a darkness cut;
overlays (`sem/overlay_*_5kx.png`, red = pore, blue = skin) confirm the split.

The headline number is normalised to the **granular (bead-packed) region** so it
is independent of how much smooth substrate/background is in frame (the 5k views
differ — up3_1 is one large particle face, up3_2 is several granules on a stub):

| cond | **open space between spheres** | (of whole surface) | skin | granular frac |
|------|-------------------------------|--------------------|------|---------------|
| up3_1 | **43.7 %** | 26.3 % | 11.1 % | 60.3 % |
| up3_2 | **42.8 %** | 25.9 % | 14.4 % | 60.5 % |
| up3_3 | **44.5 %** | 27.1 % | 10.6 % | 61.0 % |
| up3_4 | **39.1 %** | 24.0 % | 15.0 % | 61.5 % |

**Surface open-porosity ≈ 40–45 % (mean ~42.5 %)**, fairly flat across conditions;
up3_4 is the most closed (39 %, highest skin 15 %), up3_3 the most open (44.5 %).

**Vs the model:** `phi_porosity_z ≈ 0.08` is *bulk/structural* porosity; the SEM
~42 % is the *2-D projected surface* void between the outer bead layer — a looser,
projected quantity, so >> the bulk value by definition. The qualitative reading
still matters: the surface is a very open (~42 %), lightly-fused (~12 % skin) bead
packing, so the surface-fusion route `theta_skin_fus` should calibrate LOW and
`phi_porosity` may under-represent a genuinely porous particle. Resolving bulk vs
surface porosity needs a fractured/FIB cross-section (data request).

Thresholds (`DARK_PCT`, `SMOOTH_PCT`) shift the absolute value a few %; the ~42 %
level and the cross-condition ranking are the robust, DoE-usable outputs.

## Caveats
- τ (smoothness threshold) and the flattening σ are tunable; **relative** ranking
  across conditions is robust, **absolute** θ_SEM needs a reference pair (a known
  fused and a known granular sample) to fix the scale.
- 500× gives n = 130–180 particles per condition (adequate for median sphericity);
  touching particles are watershed-split, which can slightly bias solidity low.
- φ_SEM / θ_SEM are surface metrics (top of the particle only); internal structure
  needs cross-sectioned (FIB or fractured) samples — a data request if internal
  porosity is to be calibrated.
