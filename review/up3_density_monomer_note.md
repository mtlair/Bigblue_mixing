# UP3 density re-analysis — residual monomer and cake gas

**Date:** 2026-07-23 · Follow-up to `HANDOFF.md` / `DATA_REVIEW.md` Part 9, using the
full 44-row historical set in `data/up1_2_3_solids_rho.xlsx` (not just the 4 clean
chain conditions).

## Monomer: free-gas holdup in the *upstream feed*, via PV=nRT

**Process topology (confirmed).** `up1_rho` is measured on the stream *entering*
UP1 — colloid + gaseous monomer + liquid (mostly water) — **upstream of the foaming
air** (`up1_scfh` is an air sparge added inside UP1, and is *not* in `up1_rho`).
Two confirmed facts frame the estimate:

1. **Bare colloid = 1.70 g/cc, polymer** (no mineral loading); equals the dried-
   product skeletal density in the model (`up1_mixer_module.R:234`). So the old
   DATA_REVIEW story — "1730 sits below a mineral skeleton diluted by monomer" — is
   wrong; ρ_s = 1700 is just the colloid.
2. **Monomer is gaseous at room T (MW ≈ 70).** As free gas in the feed it occupies
   gas-phase volume, so it lowers `up1_rho` below the colloid+water baseline. That
   deficit → moles → mass by the ideal-gas law is the correct "guess":

```
ρ_theo   = 1/(s/ρ_s + (1−s)/ρ_w)            # colloid(1700) + water, no gas
ρ_gas    = P·MW/(R·T)                        # ≈ 2.86 kg/m³ at 1 atm, 25 °C
x_mono   = (1/ρ_meas − 1/ρ_theo)/(1/ρ_gas − 1/ρ_w)   # mass frac of feed
```

### Result: free gaseous monomer in the feed is **< ~30 ppm (effectively none)**

| cond | up1_rho | colloid+water | free monomer |
|------|---------|---------------|--------------|
| up3_1 | 1117.9 | 1112.0 | **−14 ppm** |
| up3_2–4 | 1119.1 | 1112.0 | **−16 ppm** |
| historical (n=28) | 1085–1122 | 1112 | median −14, range −23…+64 ppm |

The deficit is ~zero (in fact `up1_rho` sits a few kg/m³ *above* colloid+water), so
the implied free monomer is nil to slightly negative.

**Why the bound is so tight — the hypersensitivity.** Monomer gas is ~350× lighter
than water, so **1 kg/m³ of density deficit = only 2.3 ppm of free monomer** (1 atm,
25 °C; 4.6 ppm at 2 atm). Equivalently, **100 ppm of free gaseous monomer would drop
`up1_rho` by ~43 kg/m³** (1118 → 1075). Nothing remotely like that is seen — the feed
density hugs the gas-free colloid+water line — so free gaseous monomer is bounded to
a few tens of ppm at most.

**The one caveat is the baseline, not the method.** Because the real deficit is only
single-digit kg/m³, the sign flips with a ±5 kg/m³ move in the assumed colloid+water
baseline (water density at true feed T, dissolved surfactant/solids, exact solid
fraction). To turn the bound into a signed number, pin the baseline with a *measured*
monomer-free feed density on the same rig (or the feed T + liquid-phase density).

**What this does *not* catch:** monomer that is molecularly **dissolved** in the
water/polymer occupies liquid-like volume, not gas volume, so it leaves no PV=nRT
signature and is invisible to density at any realistic level (needs TGA / GC). The
PV=nRT feed method bounds only the *free-gas* fraction.

### Aside — the foaming-air balance (not monomer)
For completeness, applying the same ideal-gas conversion to the **foaming air**
(`up1_scfh`, added inside UP1) and the downstream void fractions: ~0.60 wt % of
solids goes in as air, ~0.05–0.08 wt % is retained as trapped air in the exit foam /
UP3 cake (≈85–90 % escapes). This is an *air* balance and is unrelated to monomer.

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
