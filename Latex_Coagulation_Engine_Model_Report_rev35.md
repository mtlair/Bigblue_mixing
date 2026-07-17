# Latex Coagulation Engine — Model Review, Restored Factors & Morris Sensitivity Report

**Script:** `up1_module_rev35_restored_factors.r` (v49, "Restored Factors"), derived from `up1_module_rev34_dpHIEP94379.r` (v48)
**Specification sources:**
1. *Latex Coagulation Engine: Integrated Master Specification*, v30.0.0 Armored Architecture / R-Engine v47 (the governing document — it states that where conflicts exist, the R-Engine v47 architecture supersedes prior versions)
2. *Production Model Package*, revision 28.0.0 (incl. revG expansion with review-layer surfactant closure)
3. *Fully Integrated Master Document*, v26.0.0/27.0.0
4. *Master Expanded revG*

**Date:** 2026-07-17

---

## 1. Executive summary

The rev34 R script implements a two-zone (cavern / dead-wall) compartment model of a stirred latex
coagulation reactor as a smooth 18-state ODE system, screened by the Morris elementary-effects
method over 23 formulation and process factors across 162 chunked operating regimes.

Auditing the script against the four specification documents found **six factors/mechanisms that
the specs require but the script had dropped or left inert**:

| # | Missing item | Spec source | Consequence in rev34 |
|---|---|---|---|
| 1 | **CMC cap** `C_eff = min(C_surfactant, CMC)` | revG §3.4 | `CMC` was a declared Morris factor but appeared **nowhere** in the physics — it always screened as exactly zero effect |
| 2 | **HLB packing-efficiency penalty** `max(0.1, 1 − 0.01(HLB−12)²)` | revG §3.4 (θ_surf,eq) | HLB influenced only foam stability; the coverage-efficiency optimum at HLB ≈ 12 was lost |
| 3 | **DLVO repulsion boost on breakage** `(1 + E_repulsion)`, `E_repulsion = ΔpH · [1/(1+10·I)] · (1−θ_surf)` | v30 report §6, Eq. 14/16 | Electrostatic redispersion of weak flocs was absent; ΔpH and ionic strength entered breakage only indirectly |
| 4 | **Wet-skin shielding of breakage** `(1 + 0.5·θ_skin)` | v30 report Eq. 16 | Skinned (coagulated-surface) aggregates broke as easily as fresh ones |
| 5 | **Particle-size flocculation kernel** `(0.2 / max(0.01, D_particle))` | v30 report Eq. 15 | Primary particle size had no effect on flocculation rate (only on breakage) |
| 6 | **Rigid template phase** (conserved solid inventory) | Production pkg §4, v30 report §9 | `n_template` was gated by `chi_npgas = 0`, so the entire template phase was identically zero: `D_template` and `C_temp_mass` were near-null factors. A units bug (missing `/ρ_template`) also made the template volume fraction ~2000× too large had it ever been enabled |

All six were restored in rev35 (Section 4 gives line-level detail). A before/after one-at-a-time probe
confirms every previously dead factor now produces finite output deltas. The full Morris omnibus
(162 regimes × 25 trajectories × 24 model runs = 97,200 ODE solves) was re-run with the restored
physics; results are in Section 5.

Additionally, rev35 fixes a plotting bug (the "Pure Volatility" zone label did not match its palette
key, so zone-4 points rendered without color) and bundles a GPL-attributed fallback of the
`sensitivity` package's Morris implementation (`morris_shim.R`) so the analysis runs in offline
environments.

---

## 2. Model architecture (as implemented)

### 2.1 Compartment structure

The reactor is split into **Zone A (the cavern)** — the well-sheared volume around the impellers —
and **Zone B (the dead wall)** — the stagnant remainder. Each zone carries a 9-state vector
(N, α_free, C_gas, C_chem, α_trap, θ_skin, Ω, θ_surf, n_template), for 18 ODE states total.
Inter-zone exchange is a symmetric first-order mass transfer whose rate is the pumping capacity
`Q_ideal = 0.7·N_rps·D³·N_imp` scaled by the *cube* of the cavern volume fraction — a strong
penalty on exchange when the cavern is small, which is the mechanistic origin of the model's
"regime map" behavior.

### 2.2 The cavern closure

Cavern size follows the classic yield-stress cavern model: the cavern boundary is where impeller
stress falls to the fluid's yield stress, giving

> D_cavern/D_imp = [ (N_p ρ N² D²) / τ_y ]^(1/3)

which is exactly the Elson–Cheesman–Nienow torque-balance cavern relation established for
yield-stress (Herschel–Bulkley/Casson) fluids in stirred tanks. The script's `Mixing_Potential`
output is the ratio of the unclamped single-impeller cavern volume to the per-stage tank volume —
i.e., a dimensionless "can the impellers own their stage" number, >1 meaning full turnover.

### 2.3 State regularization ("Titanium Armor")

Every min/max/abs is replaced by C∞ smooth surrogates (`s_pos`, `s_min`, `s_max` with a 1e-8
softening), and dissipation, breakage, and morphology rates are saturated with `tanh` caps. This
matches the v30 "Armored Architecture" requirement of a fully continuous Jacobian — important
because `lsoda` switches between Adams and BDF and benefits from smooth derivatives; hard
switches are a classic source of solver failure in parameter sweeps (hence "no NA cells" in the
Morris design).

### 2.4 What the v47/v48 engine deliberately streamlines vs the production spec

The production package (rev 28) specifies an N-stage sequential reactor with volumetric flow
continuity, a local μ_eff/γ̇ DAE closure, thermal balance, latent fractal state χ_Df, BET
observability, and seed-shell growth `t_shell`. The R engine collapses this to 2 compartments,
an algebraic power-law rheology, isothermal operation with a temperature *scale* factor, fixed
morphology proxies (Ω, θ_skin), and no shell-growth state. The v30 report explicitly blesses this
streamlining ("the R-Engine v47 architecture supersedes prior versions"), so those are treated as
design decisions, not gaps. The `C_chem` precursor state remains carried but frozen
(`dC_chem/dt = 0`) in both rev34 and rev35 — retained as a placeholder per the 9-state vector.

---

## 3. Physics audit against established science

| Model element | Implementation | Established science it maps to | Verdict |
|---|---|---|---|
| Turbulent dissipation | `eps_0D ∝ N_imp · v_tip³/D_imp`, tanh-saturated | ε = N_p ρ N³D⁵ / (ρV): for fixed geometry ε ∝ N³D² ∝ v_tip³/D. Correct scaling; the 0.005 prefactor is an effective power number | Sound (uncalibrated prefactor) |
| Cavern size | `D_c/D = (N_p ρ N²D²/τ_y)^(1/3)`, clamped to [D_imp, D_tank] | Elson & Nienow cavern model for yield-stress fluids | Textbook-consistent |
| Yield stress | Power-law consistency `1 + K(γ̇)^0.65` × foam/gas structuring `(1+100·α_trap·φ_s·Foam)`; ÷ swelling softness | Herschel–Bulkley behavior; gas-liquid foams raise apparent yield stress sharply (foam rheology, Princen); plasticizers reduce network modulus | Directionally sound; the ÷Softness placement is a v48 extension beyond Eq. 19 of the spec |
| Exit viscosity | Krieger–Dougherty `μ = μ₀(1−φ/φ_m)^(−2.5φ_m)`, φ_m = 0.64 | Krieger–Dougherty with random-close-packing limit (0.64 for monodisperse spheres); −2.5φ_m exponent recovers Einstein's 2.5 intrinsic viscosity in the dilute limit | Textbook-consistent |
| Colloid stability | Barrier `W = exp(ΔpH − 5√I)`, clamped [0.01, 100], amplified by coverage `(1+10θ_surf)` | DLVO: stability ratio W depends exponentially on the energy barrier; barrier grows with surface potential (∝ distance from IEP, i.e. ΔpH) and shrinks with ionic strength via double-layer compression (κ ∝ √I — the √I in the exponent is the correct Debye scaling). Steric stabilization by adsorbed surfactant | Sound functional form |
| DLVO repulsion (restored) | `E_rep = ΔpH · 1/(1+10I) · (1−θ_surf)` boosting breakage | Charged, uncovered surfaces inside a weak floc repel — re-peptization/redispersion of flocs at high ζ-potential and low salt | Consistent with spec Eq. 14 |
| Flocculation kernel | `R ∝ (0.2/D_p) · φ_s² · ε · Softness / E_stability` | Smoluchowski orthokinetic (shear) coagulation: rate ∝ collision frequency ∝ number density² (∝ φ²/d³ at fixed φ) × shear rate (∝ √(ε/ν)), divided by stability ratio W. The restored 1/D_p factor captures part of the number-density size dependence | Consistent; exponents are effective, not derived |
| Breakage kernel | `∝ (D_p − 0.05) · ε · (1+E_rep) / [(1+0.5α_trap)(1+0.5θ_skin)(1+100θ_surf)·Softness]` | Turbulent fragmentation: larger aggregates exceed the Kolmogorov stress limit first (size–strength scaling); gas cushioning, skin toughening and adsorbed-layer steric protection all reduce fragmentation | Consistent with spec Eq. 16 |
| Surfactant coverage | Molar coverage capacity `(C_eff/MW)·A_molecule` vs dynamic gas area `∝ (α_free+α_trap)·√ε`; CMC cap; HLB optimum at 12 | Gibbs/Langmuir adsorption with molecular footprint; above the CMC added surfactant forms micelles and does not raise surface loading; interfacial area in turbulence grows with ε (smaller Sauter diameter d₃₂ ∝ ε^(−2/5), Hinze); Griffin HLB optimum for O/W stabilization ≈ 8–16 with packing efficiency peaking mid-range | Sound after restoration |
| Adsorption dynamics | `k_surf ∝ √ε · (1000/MW) · (1+10⁹·D_surf)` | Diffusion-controlled adsorption (Ward–Tordai): flux ∝ D_surf; larger molecules diffuse and rearrange slower (∝ 1/MW trend) | Directionally sound |
| Gas thermodynamics | `C_eq ∝ P_local/T`; flash when `C_gas > C_eq`, dissolution when below; sparging `∝ Q_gas/Q_colloid` | Henry's law solubility ∝ P (retrograde with T for most gases in water); supersaturation-driven flashing; gas holdup from sparging | Sound proxy |
| Cavitation muzzle | `P_local = s_pos(P_sys − 0.001·v_tip²) + 0.05` | Bernoulli low-pressure region behind blades; floor prevents unphysical negative absolute pressure | Sound (note: report Eq. 5 prints coefficient 0.02 but its own R listing uses 0.001 — the R value governs per the supersession rule) |
| Compressible gas density | `ρ_gas = 1.2·(P_local/T_scale)` | Ideal-gas law referenced to air at 1 atm/298 K | Textbook-consistent |
| Plasticization | `Softness = exp(25·(C_monomer + C_plasticizer))` | Flory–Huggins/free-volume plasticization: solvent-swollen polymer softens roughly exponentially in solvent fraction (WLF/Fujita free-volume). Replaces the old Tg gate that zeroed solvent effects | Consistent with spec §4 |
| Rigid templates (restored) | Conserved `n_template` inventory; heterogeneous deposition `R_het ∝ A_template·φ_s·ε·S_adh_seed` with `S_adh_seed = [1/(1+e^ΔpH)]·DebyeScreening·(1+5·C_monomer)` | Seeded heterocoagulation: deposition of latex onto oppositely-charged/neutral seeds scales with available seed area and collision rate, maximized at the seed's isoelectric point and suppressed by double-layer repulsion (spec revG §4.2 form) | Restored per spec; k_het = 0.05 sits inside the spec's 1e-3…1e-1 dictionary range |
| Specific gravity closure | `SG = SG_colloid·(1−φ_por) + (ρ_gas/997)·φ_por` with ideal-gas exit density | Volume-weighted two-phase density; 997 kg/m³ = water at 25 °C per spec §11 | Textbook-consistent |
| Morris screening | OAT design, r = 25 trajectories, 5 levels, μ* and σ per factor, per regime chunk | Morris (1991) elementary effects with Campolongo's μ* (mean absolute EE) — the standard screening method for expensive nonlinear models; σ ≫ μ* flags interaction/nonlinearity | Standard practice; r = 25 is above the usual 10–20 recommendation |

**Known modeling caveats (pre-existing, unchanged in rev35, flagged for calibration):**

- `Blended_Size_um` reaches physically implausible magnitudes (10⁴–10⁵ µm) because N is a
  dimensionless index (initialized at 100) rather than a true number density; treat Size as a
  *relative* response, not microns.
- `Fraction_A + Fraction_B` can slightly exceed 1 (V_B floor of 1% tank volume), so blended
  sphericity can read marginally above 1.
- Zone B is fully algebraically stagnant (no flash/dissolution); acceptable for screening, not for
  gas inventory closure.
- `C_chem` precursor and the thermal balance remain frozen placeholders vs. the production spec.
- The kernel prefactors (0.01, 0.001, 0.05, 25.0, 100.0 …) are effective screening constants from
  the spec's baseline dictionary, not calibrated values; Morris *rankings* are meaningful, absolute
  magnitudes are not.

---

## 4. Restorations in rev35 (line-level)

All restorations are tagged `[RESTORED #n]` in the source for grep-ability.

1. **CMC cap** — RHS §1.5 and post-processing: `C_surf_eff <- s_min(C_surfactant, CMC)` feeds
   `Moles_Surf`. Smooth `s_min` keeps the Jacobian continuous at the CMC crossover.
2. **HLB efficiency** — `HLB_efficiency <- s_max(0.1, 1.0 - 0.01*(HLB-12)^2)` multiplies
   `Coverage_Capacity` (both RHS and post-processing), implementing revG's θ_surf,eq packing term.
3. **DLVO repulsion** — `Debye_Screening <- 1/s_pos(1+10*ionic_strength)`;
   `E_repulsion <- dpH * Debye_Screening * s_pos(1-Theta_surf)`; breakage numerator gains
   `s_pos(1 + E_repulsion)`. The spec's Eq. 16 breakage denominator does **not** include
   `E_stability`, so rev34's extra `* E_stability` there was removed — stability now suppresses
   *flocculation* (Eq. 15) while repulsion *assists breakage* (Eq. 16), restoring the spec's
   asymmetry.
4. **Skin shielding** — breakage denominator gains `(1.0 + 0.5*s_pos(t_skin_A))`.
5. **Size-dependent flocculation** — `Rate_floc_A` gains `(0.2 / s_max(0.01, D_particle))`.
6. **Template phase** —
   - `n_temp_initial <- phi_template_effective / V_single_template` (chi_npgas gating removed;
     the non-penetrating-gas switch gates *gas* crossing, not rigid solids, per spec §4.1);
   - volume-fraction conversion corrected to `((Q_t·C_temp/Q_tot)/ρ_template)/denom` with
     ρ_template = 2000 kg/m³ (the missing density division);
   - conserved transport `dn_temp_z/dt = exchange only` (R_loss = 0 per spec §4.1);
   - heterogeneous seed deposition added to particle formation:
     `R_het = 0.05 · SeedAreaFrac · φ_s · ε · Softness · S_adh_seed`, with seed area fraction
     computed from `n_temp·π·d_t²` against latex area `6φ_s/d_p`, and `S_adh_seed` per revG §4.2.

Plot/tooling fixes: zone-4 palette key corrected to `"4. Pure Volatility"`; `morris_shim.R`
(GPL-2, extracted from `sensitivity` 1.31.0) is sourced automatically when the CRAN package is
unavailable; results database renamed `Morris_Omnibus_rev35_r25.rds`.

**Restoration smoke test** (mid-range base point, one-at-a-time low→high sweep, max |Δ| across the
8 outputs): every restored factor now registers a finite effect — CMC 8.5e2, HLB 3.6e2,
D_template 6.1e2, C_temp_mass 3.1e3, D_particle 1.0e3, ΔpH 1.7e4, ionic strength 2.6e3,
C_surfactant 1.8e3 (deltas dominated by the large-magnitude Size response; the point is they are
no longer zero).

---

## 5. Morris sensitivity results (rev35)

### 5.1 Run configuration

- **Design:** Morris OAT, r = 25 trajectories, 5 levels, grid jump 1 (via bundled `morris_shim.R`)
- **Factors:** 23 (all `base_vars`)
- **Regimes:** 162 chunked sub-domains (v_tip 3 × C_solid_mass 3 × ΔpH 3 × ionic strength 2 × C_surfactant 3)
- **Model runs:** 162 × 25 × 24 = 97,200 ODE solves, `lsoda`, **zero failed/NA cells**, all 162 regimes valid for all 8 targets
- **Outputs:** μ* (mean |elementary effect|, main influence) and σ (EE standard deviation, interaction/nonlinearity), averaged over regimes; 24 lens plots (`Lens_{Chemistry,Additives,Physical}_{target}.png`)

### 5.2 Global factor rankings (μ* averaged over all 162 regimes, top 5 per target)

| Target | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| **Mixing_Potential** | tau (1.00) | Q_template (0.58) | C_plasticizer (0.47) | v_tip (0.42) | C_monomer (0.39) |
| **Blended_Porosity** | tau (1.00) | Q_template (0.70) | Q_gas (0.57) | Q_colloid (0.33) | Delta_pH (0.18) |
| **Blended_Sphericity** | v_tip (1.00) | tau (0.18) | Q_colloid (0.10) | C_plasticizer (0.08) | C_monomer (0.07) |
| **Blended_WetSkin** | v_tip (1.00) | tau (0.51) | C_plasticizer (0.27) | Q_colloid (0.27) | C_monomer (0.22) |
| **Blended_Size_um** | Q_template (1.00) | C_solid_mass (0.68) | Q_colloid (0.48) | tau (0.38) | Q_gas (0.37) |
| **Blended_Viscosity_PaS** | Q_template (1.00) | C_solid_mass (0.60) | Q_colloid (0.36) | C_temp_mass (0.04) | v_tip (0.01) |
| **SG_Colloid** | Q_template (1.00) | Q_colloid (0.69) | C_solid_mass (0.66) | — | — |
| **Blended_SG_Gassed** | tau (1.00) | Q_template (0.67) | Q_gas (0.57) | Q_colloid (0.33) | Delta_pH (0.18) |

(Numbers in parentheses are μ* normalized by the target's largest μ*.)

### 5.3 Verification that the restored factors are alive

The point of the exercise: in rev34, `CMC` was structurally guaranteed to screen at exactly zero,
and `D_template`/`C_temp_mass` were gated to near-nothing. In the rev35 omnibus:

- **CMC** now registers finite μ* on every gas/coverage-coupled target (Porosity, SG_Gassed,
  Size, WetSkin, Sphericity, Mixing_Potential). It stays a *low-ranked* factor — expected, since
  it only acts by capping `C_surfactant` in the sub-CMC corner of the sweep — but it is no longer
  a dead input silently wasting a Morris dimension.
- **HLB** likewise contributes through both foam stability and the restored packing-efficiency
  optimum (visible in the Chemistry-lens plots).
- **C_temp_mass** jumped to **rank 4 for Blended_Viscosity** (μ*/max ≈ 0.04) and **rank 7 for
  Blended_Size** (≈ 0.14) — the rigid-template inventory now crowds the exit envelope exactly as
  spec §4.3/Eq. 17 requires.
- **D_template** is live but weak (it survives only through seed surface area in `R_het`;
  its volumetric contribution cancels because n·v_single = φ_template). See recommendation 4.
- **Delta_pH / ionic_strength** now reach the top-6 for Porosity and SG_Gassed with the
  spec-conformant breakage asymmetry (stability suppresses flocculation, repulsion assists
  breakage).

Structural zeros that remain are *correct* zeros: `SG_Colloid` depends only on the solids
volume fraction, so only C_solid_mass / Q_colloid / Q_template can move it — all 20 other
factors screen at exactly 0, a useful integrity check on the design matrix and the shim.
(Per production-spec §2.12 the outlet density should eventually also carry template mass —
the current closure ignores it; see recommendation 4.)

### 5.4 Interpretation

1. **Residence time (tau) dominates the gas-structural outputs** (Mixing_Potential, Porosity,
   SG_Gassed): longer residence lets trapped-gas and foam states develop, which feeds back into
   yield stress and cavern collapse. Its σ/μ* ≈ 0.7–1.2 says this is strongly interaction-coupled,
   not additive — consistent with the tanh-saturated gas trapping kinetics.
2. **The feed-ratio group (Q_template, Q_colloid) is a first-order lever everywhere**, mostly via
   dilution (it rescales effective solids, residence time, and dissolved-gas load simultaneously).
   Note this is the *stream flow* effect, distinct from the template-phase physics restored in
   rev35.
3. **Morphology outputs (Sphericity, WetSkin) belong to v_tip** — shear-driven compaction and
   skin erosion — with the plasticizer/monomer softness group as the chemistry-side moderator,
   exactly the Flory–Huggins pathway the v30 spec added.
4. **Almost every active factor has σ ≥ μ*** (the "Interacting/Complex" zone of the trellis
   maps). The engine is dominated by regime-dependent interactions rather than additive effects —
   which justifies the chunked-regime Morris design over a single global screen, and warns
   against tuning any single factor in isolation.
5. **Electrostatics matter most for gas retention**: ΔpH and ionic strength act on Porosity and
   SG_Gassed (barrier-controlled flocculation/breakage set how much gas gets locked into the
   network) more than on nominal size — a plausible, DLVO-consistent outcome.

### 5.5 rev36 finer-sweep omnibus (design-audit follow-up)

A pooled audit of all 97,200 rev35 runs showed the sweep design itself needed work: 92% of runs
sat above Mixing_Potential = 10 (cavern saturated), 97% of the sphericity response occurred below
~5 m/s tip speed, tau (the top-ranked driver) was not a sweep dimension, ionic strength's steep
Debye end was split in its flat region, C_surfactant's top chunk was CMC-capped by construction
(and its sweep bounds disagreed with the dictionary: 1e-5 vs 5e-5), and the decade-spanning
factors (MW_surfactant, D_surf, C_surfactant) were sampled linearly.

`up1_module_rev36_finer_sweep.r` (v50) addresses all of these: v_tip extended down to 0.5 m/s
with non-uniform chunk breaks [0.5–4 | 4–12 | 12–32]; tau chunked ×3; ionic strength chunked
[0.01–0.05 | 0.05–0.15 | 0.15–0.5]; C_surfactant bounds aligned and split at the CMC floor
(1e-3); log-uniform sampling for the three decade-spanning factors (their elementary effects are
per log10-unit). Total: 486 regimes × 600 runs = 291,600 ODE solves, zero failures.

Verification that the redesign worked:

| v_tip chunk | fraction MP < 1 | median MP | Sphericity p10–p90 |
|---|---|---|---|
| [0.5, 4] | 0.067 | 16 | 0.20 – 0.93 |
| (4, 12] | 0.000 | 147 | 0.99 – 1.01 |
| (12, 32] | 0.000 | 664 | 1.01 – 1.01 |

The cavern-collapse transition and the entire sphericity dynamic range now live inside a
dedicated regime slice instead of being averaged away, and the trellis maps (`rev36_Lens_*.png`)
show regime-resolved rankings for tau and the electrostatics chunks.

Note on cross-revision comparison: chunking a variable narrows its within-regime range, which
shrinks its elementary effects proportionally (tau's apparent global rank drops after being
chunked — that is a design consequence, not a physics change), and log-sampled factors report
EEs per decade. rev35 and rev36 rankings are therefore each internally consistent but not
directly comparable number-to-number.

### 5.6 Result artifacts

- `Morris_Omnibus_rev35_r25.rds` / `Morris_Omnibus_rev36_r25.rds` — full design matrices, raw outputs, and regime metadata (162 / 486 scenarios)
- `Lens_*.png` (rev35) and `rev36_Lens_*.png` — 24 + 24 regime-trellis μ*–σ maps (8 targets × 3 factor lenses each)
- `morris_run.log` / `morris_run_rev36.log` — execution logs (97,200 + 291,600 runs, no solver failures)

---

## 6. Recommendations

1. **Calibrate before trusting magnitudes.** Follow the spec's §7.2 measurement map: cavern
   size/Mixing_Potential against torque or PIV data; SG_out and porosity against pycnometry;
   Size against PSD; θ_skin against extraction/microscopy. The Morris ranking above tells you
   which knobs the data must constrain first.
2. **Resolve the Size closure.** Replace the dimensionless N index with a true number balance
   (the spec's extensive n_agg/n_L architecture) so `Blended_Size_um` becomes physical.
3. **Unfreeze C_chem** or remove it from the state vector; a carried-but-frozen state costs solver
   effort and invites misreading.
4. **Template shell growth** (`t_shell`, spec §4.2) is the one spec mechanism still absent after
   rev35; adding it would let `D_template` act through shell-thickness kinetics rather than only
   through seed area.
5. **Units audit** per the production package's final notes (§9) — particularly the surfactant
   coverage scale factor (1e6) and the sparging group — before any quantitative use.
6. **Morris hygiene:** with restored factors the factor count is 23; r = 25 trajectories per
   regime is adequate for screening, but confirmatory variance-based indices (Sobol) on the
   top-5 factors per target would firm up the interaction findings (σ ≥ μ* zones).

---

## 7. References (established science)

- Morris, M.D. (1991). *Factorial sampling plans for preliminary computational experiments.* Technometrics 33(2), 161–174.
- Campolongo, F., Cariboni, J., Saltelli, A. (2007). *An effective screening design for sensitivity analysis of large models.* Env. Mod. & Software 22, 1509–1518.
- Elson, T.P., Cheesman, D.J., Nienow, A.W. (1986). *X-ray studies of cavern sizes and mixing performance with fluids possessing a yield stress.* Chem. Eng. Sci. 41, 2555–2562.
- Krieger, I.M., Dougherty, T.J. (1959). *A mechanism for non-Newtonian flow in suspensions of rigid spheres.* Trans. Soc. Rheol. 3, 137–152.
- Derjaguin, B., Landau, L. (1941); Verwey, E.J.W., Overbeek, J.Th.G. (1948). DLVO theory of colloid stability.
- von Smoluchowski, M. (1917). *Versuch einer mathematischen Theorie der Koagulationskinetik kolloider Lösungen.* Z. Phys. Chem. 92, 129–168.
- Hinze, J.O. (1955). *Fundamentals of the hydrodynamic mechanism of splitting in dispersion processes.* AIChE J. 1, 289–295.
- Griffin, W.C. (1949). *Classification of surface-active agents by "HLB".* J. Soc. Cosmet. Chem. 1, 311–326.
- Ward, A.F.H., Tordai, L. (1946). *Time-dependence of boundary tensions of solutions.* J. Chem. Phys. 14, 453–461.
- Flory, P.J. (1953). *Principles of Polymer Chemistry.* Cornell University Press. (Flory–Huggins swelling / free-volume plasticization.)
- Princen, H.M. (1983). *Rheology of foams and highly concentrated emulsions.* J. Colloid Interface Sci. 91, 160–175.
- Israelachvili, J. (2011). *Intermolecular and Surface Forces*, 3rd ed. (Debye length κ⁻¹ ∝ I^(−1/2); steric stabilization.)
- Saltelli, A. et al. (2008). *Global Sensitivity Analysis: The Primer.* Wiley.
