# Review of `HANDOFF.md` (main)

**Reviewed:** 2026-07-23 · **Scope:** the session handoff merged into `main` via PR #4
(`1f8e663 Add session handoff document`).

Every quantitative claim in the handoff was checked against the committed code,
data, and Morris outputs. The handoff is accurate and well-supported. One
substantive data-integrity issue and a handful of cosmetic nits are noted below.

---

## Verified accurate

| Claim | Check | Result |
|-------|-------|--------|
| theta_skin two-route closure (Péclet + surface-fusion) | `unified/up2_atomizer_dryer_module.R:262–290` | matches the snippet |
| `SKIN_SURVIVAL = 0.30` seeding | `unified/interface_stream.R:59` | ✓ |
| C_monomer / C_plasticizer → surface_fusion at ranks **14 / 16**, positive μ* | ranked `unified_morris_indices.csv` for `up2_theta_skin_z` | ✓ exact |
| Morris: 41 factors, r=30 → **1260** runs, 1260/1260 valid | `(41+1)×30 = 1260`; `morris_run.log` | ✓ |
| Top-5 driver tables (all 5 outputs) | `unified_morris_indices.csv` vs `nominal_chain_summary.txt` | ✓ match (spot-checked D_particle, X_moisture) |
| Solid particle density 1730–1742 kg/m³ | `cond_up1234.csv` → 1729, 1742, 1742, 1742 | ✓ |
| UP3 total mass balance 249.84 ≈ 152.8 + 97.2 | cake+reject = 250.0 | ✓ (and solid balance closes: 62.2+0.28 ≈ 62.5) |
| UP2 wash water (49 lb/hr) bypasses separator | UP3 balance closes on UP1 flow alone | ✓ consistent inference |
| UP3 cake gas holdup 7.6% | (1209−1117.3)/(1209−1) = 7.59% | ✓ |
| UP1 exit density (both estimates) | `cond_up1234.csv` cols match text | ✓ |
| UP4 feed table (mb / heat-duty / efficiency) | efficiency = mb/heatduty, matches CSV | ✓ |
| All files named in the directory tree | present on disk | ✓ |

The `theta_skin` sign fix is real and does what the handoff says: after the fix
`C_monomer`/`C_plasticizer` carry positive μ* for `up2_theta_skin_z` (more
solvent → more surface_fusion), consistent with the stated physical rationale.

---

## Substantive issue — SEM sample numbers do not match the images on disk

**Next-step #5** and the directory listing cite SEM sample numbers **114641–114644**
mapping to conditions up3_1–up3_4. That transcription is faithful to the
`sem_sample_no_filename` sheet in `data/up1_2_3_4_visc_sem.xlsx`, which does list
114641→up3_1 … 114644→up3_4.

**But no image in `sem/` carries those sample numbers.** The 95 files on disk span
12 distinct samples — `X109384`, `X110041`–`X110049`, `X110494`, `X110495` — and
`SEM_X114641*` … `SEM_X114644*` do not exist. Per the sheet's own filename
convention (`SEM_X<sample no>_…`), the mapped conditions have no matching files.

Consequence: the SEM analysis next-step is **not executable as written**. Either
the images for the four chain conditions were never committed (the folder holds a
different sample set), or the mapping sheet's sample numbers are stale. This
should be resolved before anyone attempts next-step #5 — the handoff presents the
mapping as usable, and it isn't.

---

## Cosmetic nits (no action required)

1. **Code snippet vs. source.** The handoff's `theta_skin_fus` snippet writes
   `Tg_polymer` where the source clamps it: `inv_Tg_plas <- (1 - phi_solvent) /
   max(Tg_pol, 200) + …` (`up2_atomizer_dryer_module.R:278`). Math-equivalent except
   for the `Tg_pol < 200 K` guard; the snippet is a faithful simplification.
2. **Column label.** `up1_exit_rho_5pctgas_kgm3` holds the Boyle back-calc value
   (1063–1070 kg/m³) that the handoff text describes as ~4.4–4.9% gas, not 5%.
   Label rounds the gas fraction; value is correct.
3. **Header branch.** The handoff lists `claude/psd-viscosity-data-review-khqg39`
   (its authoring branch, since merged). Accurate as a historical record.
4. **Gas-free cake density.** 1209 kg/m³ vs. ~1207 recomputed with ρ_solid=1730,
   ρ_water=1000 (water at 55 °C is ~986). ~0.15% rounding; does not move the 7.6%
   holdup materially.

---

## Bottom line

The handoff is a reliable record of the session: the model change, the Morris
re-run, and every derived process number reproduce from the committed artifacts.
The only thing a follow-on session should not take at face value is the SEM
sample-number mapping — reconcile it against the actual `sem/` filenames first.
