# =============================================================================
# INTER-MODULE STREAM INTERFACE
# =============================================================================
# The unified model chains unit operations through a single "stream" object —
# a named list describing the slurry leaving one stage and entering the next:
#
#   UP1 mixer -> [intermediate stage 1] -> [intermediate stage 2] -> spray dryer
#
# The two intermediate stages are IDENTITY PLACEHOLDERS for the unit
# operations to be added later (e.g. transfer line / pump, hold tank,
# pre-heater, degasser, in-line conditioner). Each takes the stream and
# returns a stream with the same fields, so a future stage only has to
# overwrite the fields it physically changes.
#
# Stream fields
# -------------
# Composition (mass fractions of total slurry unless noted):
#   C_solid        solids (latex polymer) mass fraction, diluted by template feed
#   C_solid_rigid  rigid-template solids mass fraction (0 for liquid/gas types)
#   C_binder, C_monomer, C_plasticizer, C_surfactant   formulation additives
#   CMC, HLB, MW_surfactant [g/mol], A_molecule [nm2], Delta_pH, ionic_strength
# Physical state:
#   rho_slurry [kg/m3], rho_polymer [kg/m3], T_K [K], P_Pa [Pa],
#   mu_exit_PaS (mixer exit apparent viscosity, diagnostic)
# Particulate state:
#   D_particle_um  primary colloid particle diameter (mixer INPUT / reference)
#   D_primary_exit_um  primary size AFTER the UP1 three-regime v_tip closure:
#                  == D_particle_um below the critical tip speed, milled smaller
#                  past it. This (not D_particle_um) is the primary size
#                  downstream stages should use.
#   v_tip_crit     critical tip speed [m/s] (aggregation onset; diagnostic)
#   v_tip_ratio    v_tip / v_tip_crit (>1 aggregation, >>1 milling; diagnostic)
#   D_agg_um       aggregate (blended) size at mixer exit (enlarged in the
#                  aggregation regime, broken back down under milling)
#   sphericity     Omega at mixer exit (diagnostic)
#   WetSkin        wet-skin fraction at mixer exit (partially survives
#                  atomization; seeds the dryer skin state)
# Gas phase:
#   alpha_g        entrained gas holdup (trapped microbubbles)
#   D_b_m          bubble diameter [m]; NA until set — the transfer
#                  line / pump shear (future intermediate stage) determines
#                  it, so the dryer treats it as a screening factor for now
# Template phase:
#   template_type  1=rigid, 2=gas, 3=surface_weld, 4=capillary_bridge
#   phi_templ_free free (interstitial) liquid-template volume fraction of
#                  the slurry -> becomes the dryer's pore-templating emulsion
#   D_template_um  template droplet / seed diameter
#   RTF            Residual_Template_Fraction: fraction of the liquid
#                  template absorbed into particle CORES by mixer exit
#   w_core         core-absorbed template, mass fraction of total slurry
#                  -> acts as core plasticizer + confined volatile in dryer
# Diagnostics carried along (not consumed by the dryer closure):
#   Softness_exit, Bond_Strength, Retained_Porosity, Mixing_Potential
# =============================================================================

RHO_TEMPLATE_LIQ <- 750    # template solvent liquid density [kg/m3] (UP2 spec)
SKIN_SURVIVAL    <- 0.30   # fraction of mixer wet-skin surviving atomization
                           # (new surface is created at the nozzle) — interface
                           # assumption, revisit when intermediate stages land

stream_from_up1 <- function(up1_res, pars, equipment) {
  p  <- as.list(pars)
  ex <- up1_res$extras
  out <- up1_res$outputs

  t_type <- ex$template_type

  # Liquid template volume fraction of the slurry: template_dose is defined
  # relative to the interstitial void volume of the packed solids
  # (V_void = phi_s * eps_rcp / (1 - eps_rcp)), rev37 definition.
  phi_templ_liq <- if (t_type %in% c(3, 4)) {
    min(0.30, ex$template_fill * ex$phi_s_exit * UP1_EPS_RCP / (1 - UP1_EPS_RCP))
  } else 0.0

  # RTF splits the liquid template between free interstitial droplets
  # (pore templating in the dryer) and core-absorbed plasticizer.
  RTF <- unname(out[["Residual_Template_Fraction"]])
  phi_templ_free <- phi_templ_liq * (1 - RTF)
  w_core <- phi_templ_liq * RTF * RHO_TEMPLATE_LIQ / ex$rho_slurry

  list(
    # composition
    C_solid        = ex$C_solid_exit,
    C_solid_rigid  = if (t_type == 1) ex$C_temp_rigid else 0.0,
    C_binder       = p$C_binder,
    C_monomer      = p$C_monomer,
    C_plasticizer  = p$C_plasticizer,
    C_surfactant   = p$C_surfactant,
    CMC            = p$CMC,
    HLB            = p$HLB,
    MW_surfactant  = p$MW_surfactant,
    A_molecule     = p$A_molecule,
    Delta_pH       = p$Delta_pH,
    ionic_strength = p$ionic_strength,
    # physical state
    rho_slurry     = ex$rho_slurry,
    rho_polymer    = ex$rho_polymer,
    T_K            = ex$T_exit_K,
    P_Pa           = ex$P_exit_Pa,
    mu_exit_PaS    = unname(out[["Blended_Viscosity_PaS"]]),
    # particulate state
    D_particle_um      = p$D_particle,           # input primary size (reference)
    D_primary_exit_um  = ex$D_primary_exit_um,   # UP1 three-regime milled primary size
    v_tip_crit         = ex$v_tip_crit,           # critical tip speed [m/s] (diagnostic)
    v_tip_ratio        = ex$v_tip_ratio,          # regime indicator >1=milling (diagnostic)
    D_agg_um       = unname(out[["Blended_Size_um"]]),
    sphericity     = unname(out[["Blended_Sphericity"]]),
    WetSkin        = unname(out[["Blended_WetSkin"]]),
    # gas phase
    alpha_g        = min(max(unname(out[["Blended_Porosity"]]), 0), 0.75),
    D_b_m          = NA_real_,
    # template phase
    template_type  = t_type,
    phi_templ_free = phi_templ_free,
    D_template_um  = p$D_template,
    RTF            = RTF,
    w_core         = w_core,
    # diagnostics
    Softness_exit  = unname(out[["Swelling_Softness_exit"]]),
    Bond_Strength  = unname(out[["Bond_Strength"]]),
    Retained_Porosity = unname(out[["Retained_Porosity"]]),
    Mixing_Potential  = unname(out[["Mixing_Potential"]])
  )
}

# -----------------------------------------------------------------------------
# PLACEHOLDER STAGES — future unit operations slot in here.
# Contract: f(stream, pars = list()) -> stream (same fields, possibly updated).
# -----------------------------------------------------------------------------

# Intermediate stage 1 (e.g. transfer line / feed pump):
# would set stream$P_Pa (pump discharge), stream$D_b_m (line-shear bubble
# size), and could coarsen/shear the emulsion (stream$D_template_um).
intermediate_stage_1 <- function(stream, pars = list()) {
  stream
}

# Intermediate stage 2 (e.g. hold tank / pre-heater / degasser):
# would update stream$T_K, stream$alpha_g (venting), stream$D_b_m and
# stream$D_template_um (Ostwald ripening over the hold time).
intermediate_stage_2 <- function(stream, pars = list()) {
  stream
}

print_stream <- function(stream, label = "stream") {
  cat(sprintf("--- %s ---\n", label))
  num <- vapply(stream, function(v) is.numeric(v) && length(v) == 1, logical(1))
  for (n in names(stream)[num])
    cat(sprintf("  %-16s %s\n", n, format(signif(stream[[n]], 4))))
  invisible(stream)
}
