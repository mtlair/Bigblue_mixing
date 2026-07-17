# Liquid Templating of Colloid Particles — Physics Summary

**Purpose:** Standalone reasoning note on using a liquid template to wrap/bond colloid solid
particles, prepared for hand-off to a separate process-engine thread.
**Origin:** Latex Coagulation Engine analysis (branch `claude/r-script-missing-factors-jmlr8b`).
**File location:** `/home/user/Bigblue_mixing/Liquid_Template_Summary.md`
(repo `mtlair/bigblue_mixing`, GitHub path `Liquid_Template_Summary.md`).
**Companion files:** `up1_module_rev37_templating.r` (implementation), `Pendular_Dose_Window.csv`
and `rev37_Pendular_Dose_Window.png` (dose sweep), `rev37_templating_summary.md` (Morris rankings).

---

## 1. Target structure

Goal: **discrete colloid solids that stick to each other lightly at contact points while
preserving the interstitial void space between packed spheres.** In powder/agglomeration
terminology this is a **pendular-state agglomerate** — solid spheres joined by discrete liquid
(or, after solidification, solid) bridges at their contact points, with the pore network still
open. It is deliberately *short of* the funicular and capillary states, where liquid fills the
pores and the structure densifies.

Two physics families govern it:
- **Capillary bridging / spherical agglomeration** (Capes; Rumpf agglomerate strength) — a
  wetting bridging liquid pulls particles together with force `F ≈ 2πrγ·cosθ` per contact.
- **Viscous sintering of necks** (Frenkel) — where the bridging liquid dissolves/softens the
  polymer, colliding particles weld: neck growth `x²/r ≈ σ·t/η_local`.

---

## 2. The 2×2 template matrix (miscibility × volatility)

The template liquid is **immiscible with the bulk water phase** in all cases (prepared as a
dispersion). The two discriminating axes are (a) whether it dissolves/softens the colloid
**polymer**, and (b) its **boiling point relative to the bulk liquid** (above or below).

### Case 1 — Immiscible with water, SOLUBLE in the polymer ("surface-weld" liquid)
The liquid partitions preferentially onto particle **surfaces** and plasticizes a thin shell while
the core stays hard (transfer is surface-limited on mixing timescales). Colliding particles then
**weld** — the same mechanism as latex film formation, but deliberately starved so it stops at
necks rather than proceeding to full coalescence. Bonds are **permanent polymer necks** at
contact points.

- **BP below water (volatile) → BEST FIT for the target.** Necks form, then residual solvent
  strips out during devolatilization/drying, freezing the structure before coalescence completes.
  Result: discrete spheres, welded necks (~0.2–0.3 neck-to-radius ratio for green strength), most
  of the ~36% random-close-packing porosity retained.
- **BP above water (non-volatile) → voids lost.** Solvent stays in the polymer permanently;
  coalescence continues during storage (cold flow), and micron-scale voids collapse under
  capillary pressure (2σ/r ≈ several atm at 1 µm) during drying. Dense, soft agglomerates. Only
  viable with a downstream extraction step.

### Case 2 — Immiscible with water AND with the polymer ("capillary-bridge" liquid)
No dissolution; bonding is purely **capillary**. Critical caveat: *immiscible ≠ non-wetting.*
Wetting decides everything.

- **If the liquid wets the polymer** (contact angle through the bridge < 90°; true for e.g.
  alkanes on most latex polymers): forms **pendular capillary bridges**. Rumpf agglomerate
  strength at d = 1 µm, γ ≈ 30 mN/m: `σ ≈ (9/8)·((1−ε)/ε)·(F/d²) ≈ 10^5 Pa` — orders of magnitude
  above turbulent cavern stresses, so bridged agglomerates easily survive agitation. But the bond
  is **reversible**: remove the liquid and the glue is gone.
  - **BP below water:** bridges hold in the mixer, evaporate on drying → green body crumbles back
    to discrete spheres unless the particles have intrinsic tack (T near Tg) or a brief thermal
    consolidation step precedes stripping.
  - **BP above water:** agglomerates stay intact, but the "voids" are **liquid-filled** — porosity
    only appears after an extraction step. This is classical spherical agglomeration (particle
    rounding / granulation). Fine if the liquid is meant to stay (oil-extended product).
- **If the liquid does NOT wet the polymer:** no bridging. It rides along as inert occluded
  droplets, templating isolated inclusions rather than interparticle bonds.

### Summary table

| Template | Water | Polymer | BP | Bond type | Voids after drying | Fit to goal |
|---|---|---|---|---|---|---|
| Surface-weld, volatile | immiscible | soluble | < water | welded neck (permanent) | **preserved** | **best** |
| Surface-weld, non-volatile | immiscible | soluble | > water | neck → full coalescence | collapsed | poor (needs extraction) |
| Capillary-bridge, volatile | immiscible | immiscible (wetting) | < water | capillary (reversible) | preserved *if* consolidated | good with consolidation step |
| Capillary-bridge, non-volatile | immiscible | immiscible (wetting) | > water | capillary (reversible) | liquid-filled until extracted | good only if liquid stays or is extracted |
| Non-wetting droplet | immiscible | immiscible (non-wetting) | any | none (occlusion) | isolated inclusions | off-target |

---

## 3. The controlling variable is DOSE, not chemistry

Bond state is set by how much of the interstitial **void volume** the template liquid occupies:

| Saturation regime | Liquid / void volume | Structure | Relevance |
|---|---|---|---|
| Pendular | ~0 – 25 % | discrete bridges at contacts, **open pores** | **TARGET** |
| Funicular | ~25 – 80 % | bridges + partially filled pores | transition (voids closing) |
| Capillary | ~80 – 100 % | pores liquid-filled, surface menisci | densified granule (voids gone) |
| Droplet/slurry | > 100 % | particles suspended in template | no agglomerate |

For a random close packing (ε ≈ 0.36), interstitial void ≈ 0.56 × solid volume, so the pendular
window corresponds to **liquid ≈ 4 – 9 wt% of solids** (with polymer/template density ~1). Above
~30% void fill the structure densifies — the opposite of the goal. The engine's existing
`C_temp_mass` range (1–10 wt%) conveniently brackets this window.

---

## 4. Two process couplings to watch

1. **In-situ flashing.** The engine's operating window reaches T = 360 K and P = 0.34 atm, plus a
   Bernoulli low-pressure zone behind the impeller blades. A below-BP template can **flash-boil in
   place** (e.g. pentane boils ~10 °C at 0.34 atm) and become a *gas* template, with the
   associated bubble instabilities and a possible pressure/foaming excursion. Either design it as
   a deliberate in-situ blowing agent or keep the template comfortably above the local bubble point
   everywhere in the vessel — do not leave it to chance.
2. **Surface vs bulk softening.** For the surface-weld case, the template must add tack at the
   *interface* (raising adhesion/neck strength) **without** feeding the bulk plasticization pool.
   This is the key distinction from the generic monomer/plasticizer softening already in the model:
   generic plasticizer softens the whole matrix and lowers yield stress everywhere; a surface-weld
   template should only bond contacts. Conflating the two over-predicts coalescence and void
   collapse.

---

## 5. Recommendation

For "stick a bit, keep the voids," use a **volatile (BP < bulk water), water-immiscible,
polymer-soluble liquid at pendular dose (≈ 4–9 wt% of solids)**: welded necks that survive solvent
removal, open porosity retained. The both-immiscible variant works only with a separate
consolidation step before the bridge liquid is stripped; any above-BP template leaves either
plasticized collapse (polymer-soluble) or liquid-filled voids (polymer-immiscible) unless followed
by extraction.

---

## 6. How this maps into the coagulation engine (rev37)

Implemented as a `template_type` switch with four modes — `rigid` (baseline silica seed), `gas`
(non-penetrating gas, `chi_npgas`=1), `surface_weld`, `capillary_bridge` — plus a continuous
`template_dose` factor (template-liquid volume as a fraction of solids volume) and a pendular-fill
calculation. New model outputs:

- **`Bond_Strength`** — neck (Frenkel, surface-weld) or capillary (Rumpf, bridge) bond strength,
  0–1, rising with fill and saturating in the pendular window; ≈0 for gas, low for rigid.
- **`Retained_Porosity`** — post-drying/post-extraction void fraction; preserved up to the
  pendular/funicular boundary, then declining as fill densifies the granule (surface-weld also
  loses voids to coalescence when over-dosed and non-volatile).
- **`Template_Target_Score`** = `Bond_Strength × Retained_Porosity` — peaks inside the pendular
  window; the single-number objective for the "stick-but-porous" goal.

The rev37 deterministic dose sweep (`Pendular_Dose_Window.csv`,
`rev37_Pendular_Dose_Window.png`) shows these three curves vs dose for each template type, locating
the pendular optimum; the rev37 Morris screen (`rev37_templating_summary.md`) ranks which
formulation/process factors most move the target score within each templating strategy.

*Model caveat carried from the parent report: the engine's kernel prefactors are effective
screening constants, not calibrated values, so rankings and qualitative dose trends are meaningful
but absolute magnitudes require calibration against pendular-bridge/BET/porosimetry data.*
