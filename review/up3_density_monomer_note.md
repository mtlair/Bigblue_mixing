# UP3 density re-analysis — residual monomer and cake gas

**Date:** 2026-07-23 · Follow-up to `HANDOFF.md` / `DATA_REVIEW.md` Part 9, using the
full 44-row historical set in `data/up1_2_3_solids_rho.xlsx` (not just the 4 clean
chain conditions).

## Short answer on residual monomer

The **right comparison is on the UP1 colloid feed, not the UP3 cake.** The feed
carries no sparged gas (gas is introduced *inside* UP1), so it is a clean
three-component mix — colloid solid + water + monomer — with none of the
solid-vs-gas degeneracy that plagues the cake. Referencing the measured feed
density against the solids+water theoretical isolates monomer by volume balance:

```
ρ_theo = 1 / (s/ρ_s + (1−s)/ρ_w)                       # solid + water only
w_m    = (1/ρ_meas − 1/ρ_theo) / (1/ρ_m − 1/ρ_w)       # monomer mass frac of feed
```

I ran this over the 4 clean chain feeds and the historical set. The method is
mathematically valid but **its answer is set almost entirely by two densities the
repo does not specify** — the bare-colloid solid density ρ_s and the monomer
density ρ_m — and by the measurement noise floor. Two hard limits:

**1. The ρ_s lever is enormous.** Monomer holdup for the up3_1 feed (1117.9 kg/m³
@ 25 % solid, ρ_m = 950):

| ρ_s (kg/m³) | 1700 | 1900 | 2100 | 2400 | 2650 |
|-------------|------|------|------|------|------|
| w_m (wt %)  | −10  | +22  | +47  | +77  | +97  |

A ±5–7 kg/m³ feed-density deficit, swung over this lever, gives anything from 0 %
to ~100 % monomer. You cannot read monomer off feed density without ρ_s known to
≈ ±1 %.

**2. At the model's own ρ_s = 1700, there is no deficit — the sign is wrong for
monomer.** Predicted feed density (solid@1700 + water@997) = 1112; measured =
1117.9 → feed is **6–7 kg/m³ *denser*** than the two-component prediction. Monomer
is lighter than water, so it would make the feed *lighter*, not heavier. Referenced
to the dried-product skeletal density (1700), monomer holdup reads as ~zero (in
fact slightly negative). A positive, physical monomer number requires assuming the
bare colloid is denser than the dried product (ρ_s ≳ 1900) — plausible if drying
adds closed porosity to the skeleton, but unquantified here.

**3. Noise floor ≈ 7–15 wt % monomer.** 1 wt % monomer shifts feed density by only
0.6–1.4 kg/m³ (for ρ_m = 950–900), while the historical feed-density scatter at
25 % solid is σ ≈ 9 kg/m³. Even with ρ_s known exactly, only a *large* holdup would
clear the noise.

**Conclusion:** the feed-density comparison is the correct experiment and I did it,
but the density signal is too small and too dependent on an unknown ρ_s to promote
the "guess" to a number. Close it with an independent measurement:
- **He pycnometry** on the *bare dried colloid* → ρ_s directly (kills the lever).
- **GC / solvent-extraction or TGA** → residual monomer wt% directly.

With bare-colloid ρ_s in hand, the `w_m` formula above turns the existing feed
densities straight into monomer holdup per condition.

## What is now well-established: ρ_solid ≈ 1700 kg/m³

This is triangulated three independent ways and should be treated as settled:
1. UP1 feed density (1117.9 kg/m³ @ 25 % solid) → ρ_solid ≈ 1730.
2. Dry-powder bulk density (~0.30 g/cc) + SEM porosity work → 1.70 g/cc skeletal
   (already baked into the model: `unified/up1_mixer_module.R:234`, `rho_polymer <- 1700`).
3. Consistent with UP3 cake densities **once trapped gas is included** (below).

## What the new UP3 data *does* change

### 1. Cake gas holdup is ~14 %, not 7.6 %
The handoff's 7.6 % comes from the single up3_1 chain point, which is a **low
outlier**. Across the 23 historical UP3 cake measurements, the implied trapped gas
(at ρ_solid = 1700) is **6.0–18.6 %, median ≈ 14 %**. The 40–52 % solid cakes
measure 1035–1170 kg/m³, well below the gas-free prediction (~1255 kg/m³ @ 48 %).
Typical trapped gas is ~2× the handoff headline.

### 2. Trapped gas barely responds to centrifugal force
If the deficit were free gas, higher g-force should expel it and drive measured
density toward the gas-free value. It doesn't: `corr(Fg, ρ_solid-ignoring-gas) =
0.36`, and even the ~1900 g cakes stay at 1035–1150 kg/m³. Two readings, both
worth a follow-up:
- the gas is finely trapped inside particles/pores (survives centrifugation), consistent with UP1 sparging; and/or
- **Fg and weir are confounded** in this dataset — every high-Fg row is `weir = low`
  and every low-Fg row is `weir = hi`, so the current data cannot separate the
  g-force effect from the weir (pool-depth) effect on cake wetness. A DoE must
  cross these two.

### 3. The feed→ρ_solid inference is noisier than "1729–1742" suggests
Over the full historical set, `up1feed_rho` spans 977–1122 kg/m³ and back-solves to
ρ_solid **916–1771** (median 1716). The tight "1729–1742" is only the 4 clean
chain conditions. Many historical rows are startup/transient (near-water feed
density, zero or negative UP2 flows) and must be filtered to steady state before
any fit.

### 4. Reject loss is not universally < 0.5 %
`up3_reject_solid%` reaches 5–12 % in several historical upset rows. The handoff's
"~0.28 % solid, negligible" holds for the 4 clean chain conditions only — do not
generalize it when using the wider set.

## Bottom line
The new UP3 data firms up ρ_solid ≈ 1700 and **raises** the modelled cake gas to
~14 %, but it does not — and cannot — isolate residual monomer. Close that item
with pycnometry/TGA, not more density arithmetic. If UP3 cake moisture/gas is going
into the model, use ~13–14 % trapped gas as the central value and decouple Fg from
weir in the next DoE.
