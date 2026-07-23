# UP3 density re-analysis — residual monomer and cake gas

**Date:** 2026-07-23 · Follow-up to `HANDOFF.md` / `DATA_REVIEW.md` Part 9, using the
full 44-row historical set in `data/up1_2_3_solids_rho.xlsx` (not just the 4 clean
chain conditions).

## Resolved: monomer is a gas — quantify it with PV=nRT, not a density deficit

Two facts (confirmed by the process owner) settle the residual-monomer question:

1. **Bare colloid density = 1.70 g/cc, and the colloid is polymer** (no mineral
   loading). This equals the dried-product skeletal density already in the model
   (`up1_mixer_module.R:234`). So the DATA_REVIEW narrative — "1730 is *below* a
   pure mineral skeleton because residual monomer dilutes it" — is **wrong**: there
   is no mineral skeleton and no dilution to infer. ρ_s = 1700 is just the colloid.

2. **Monomer is gaseous at room temperature (MW ≈ 70 g/mol).** It therefore never
   appears as a *liquid* density deficit. This is exactly why the clean feed check
   shows no deficit: predicted feed density (colloid@1700 + water@997) = 1112 vs
   measured 1117.9 — the feed is if anything 6 kg/m³ *denser*, i.e. **zero
   dissolved/liquid monomer**, consistent with a gas-phase monomer. My earlier
   liquid-holdup formula was the wrong model; the correct tool is the ideal-gas law
   applied to the measured **gas holdup**.

### PV=nRT monomer holdup (MW = 70)

For a gas void fraction α at stream pressure P and temperature T,
`n = P·V_gas/(R·T)` and monomer mass `= n·0.070 kg/mol`, expressed per unit solid:

| Basis | Measurement | Gas frac | Monomer (wt % of solid) |
|-------|-------------|----------|--------------------------|
| **Sparge input** (up1_scfh as monomer) | flow meter | — | **0.59–0.60** |
| UP1 exit, all-entrained (upper bound) | up1_scfh expanded to exit P,T | ~30 % | 0.59–0.60 |
| UP1 exit, retained (back-calc from UP3) | Boyle from cake holdup | ~4.5 % | 0.06–0.07 |
| **UP3 cake, retained** | measured cake density | 7.6 % (up3_1) / ~14 % (hist) | **0.05–0.08** |

**The gas is foaming air, not monomer** (confirmed: `up1_scfh` is an air sparge).
So the numbers in the table above quantify *air*, not monomer — the ~0.60 wt %
"input" is the foaming air, and the ~0.05–0.08 wt % "retained" is trapped air in
the exit foam / UP3 cake. They stand as a useful **air** balance (≈85–90 % of the
foaming air escapes before the product exits; up3_1: ~0.37 lb/hr in, ~0.03 lb/hr
retained) but say nothing directly about monomer.

### Bottom line on residual monomer
It **cannot be quantified from the current data.** Monomer is gaseous (MW ≈ 70), so:
- it leaves **no liquid-density signature** — the feed check confirms zero dissolved
  monomer (colloid+water fully explains 1118 kg/m³); and
- it is **not** the measured void gas — that void is foaming air.

The most one can say is a loose upper bound: *if* any residual monomer shared the
trapped-air voids, it would be ≤ the ~0.05–0.08 wt % air holdup. A real number
needs data not in these columns:
- a **monomer feed / off-gas flow** (then PV=nRT with MW 70 gives it directly), or
- a **TGA / GC assay** on the product.

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
