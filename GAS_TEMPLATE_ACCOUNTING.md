# Gas & Template Accounting — Three-Stream Split

**Status:** Open design item (not yet implemented in the model)
**Applies to:** UP1 mixer → UP2 foam-wash → UP3 separator → UP4 dryer chain
**Related fields:** `alpha_g`, `D_b_m`, `template_type`, `D_template_um`, `RTF`, `w_core`
(see `unified/STAGE_CONTRACT.md` field dictionary)

---

## The problem

The template gas and the foam gas are **not the same population**, and today the
model tends to fold them into a single entrained-gas holdup (`alpha_g`). They
have different fates, different length scales, and different downstream effects,
so they must be tracked as **separate streams**.

`D_template` in particular is the size of gas (or liquid) that gets **trapped
inside the very small solid particles** as coagulum grows around it. It is an
**intra-particle** quantity — it is not the interstitial foam holdup and it is
not the gas that escapes. It must be accounted for on its own, sized by
`D_template`, not lumped into the bulk foam gas.

---

## The three gas streams (gas-template case, `template_type = 2`)

When the template is gas, the gas in the system partitions into three
physically distinct populations:

### 1. Trapped template gas — *inside* the solid particles
- Gas that seeded a coagulation nucleus; the small solid particle **grows
  coagulum around it**, locking the bubble inside as an intra-particle cavity.
- Sized by `D_template` (the trapped-cavity diameter), analogous to the
  liquid-template `RTF`/`w_core` path but for gas.
- **Fate:** stays with the solids all the way to the dryer; becomes
  intra-particle porosity / a templated void in the final powder.
- **Accounting:** a dedicated field (a gas analogue of `w_core`/`RTF`), *not*
  part of `alpha_g`. This is the stream `D_template` belongs to.

### 2. Interstitial holdup gas — *between* the particles (the foam)
- Gas held up in the interstitial space between particles: this is the foam
  holdup proper, "the gas in foam."
- This is what `alpha_g` (`Blended_Porosity`, capped 0.75) currently represents,
  and what the foam-wash column (UP2) is designed to wash out.
- **Fate:** removed/compressed in UP2; whatever survives biases downstream gas
  holdup into UP3/UP4.

### 3. Slug gas — escaping the foam system
- Large gas slugs that coalesce and break out of the foam and **escape the
  system** overhead — either as excess template gas or as foam make-up/sparge
  gas that has done its job.
- **Fate:** leaves the mixer/column headspace; not carried downstream in the
  slurry (it is vented / washed overhead in UP2).
- **Accounting:** an escape/vent term in the mixer (and UP2 wash) balance, not
  a slurry-stream field.

---

## Why the split matters

| Stream | Length scale | Location | Downstream effect | Removed by |
|---|---|---|---|---|
| Trapped template gas | `D_template` | intra-particle | templated void / intra-particle porosity in powder | (survives to dryer) |
| Interstitial foam holdup | `D_b` (bubble) | inter-particle | UP4 `alpha_g_0` seed | UP2 foam-wash |
| Slug escape | slug (mm+) | headspace / overhead | none (leaves system) | vents / UP2 overhead |

Collapsing these into one `alpha_g` misattributes gas that should stay inside
the particle (porosity template) to the interstitial foam that UP2 washes out —
so the model would (a) let the wash column strip gas it physically cannot reach,
and (b) lose the `D_template`-sized intra-particle porosity signal in the powder.

---

## Mass-balance requirement

Total gas fed to the mixer must close across the three sinks:

```
G_in  =  G_trapped(template, ~D_template)     # -> stays in solids (porosity)
       + G_holdup(interstitial foam, alpha_g) # -> to UP2, partly washed
       + G_slug(escape)                        # -> vents / overhead
```

- **Gas template:** all three terms are gas; `G_trapped` is sized by
  `D_template`.
- **Liquid template:** the "trapped" term is the existing liquid `RTF`/`w_core`
  split (template absorbed into cores); only `G_holdup` and `G_slug` are gas.
- **Rigid template:** no trapped-gas term from the template; `G_holdup` and
  `G_slug` still apply to any entrained/sparge gas.

---

## Proposed implementation (open)

1. Add a dedicated **trapped-template** stream field (gas analogue of `w_core`),
   set at UP1 from `template_type`, `D_template_um`, and the trapping fraction —
   kept separate from `alpha_g`.
2. Keep `alpha_g` as **interstitial foam holdup only**; ensure UP2 washes only
   this stream, never the trapped fraction.
3. Add a **slug-escape / vent** term to the UP1 (and UP2) gas balance so total
   gas closes across trapped + holdup + escape.
4. In UP4, route the trapped-template stream to intra-particle porosity sized by
   `D_template`, distinct from the interstitial-gas porosity path.
