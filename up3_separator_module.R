# =============================================================================
# UP3 — DECANTING SEPARATOR MODULE  (stream-contract closure)
# =============================================================================
# Concentrates the washed foam from UP1/UP2 (~25% solid) to the measured cake
# solids (~38-43%), consistent with the confirmed mass balance in DATA_REVIEW:
#   * UP2 wash water exits UP2 separately -> UP3 receives the UP1 flow unchanged
#     in total; it only CONCENTRATES (removes centrate liquid).
#   * A small solid fraction is lost to the centrate reject (<0.5%, up3_1 = 0.28%).
#   * Centrifugation removes most free gas but leaves a trapped-gas floor in the
#     cake (measured 7.6% for up3_1; ~14% median across the historical set).
#
# Conforms to the stream contract f(stream, pars) -> stream, so it slots into
# unified/interface_stream.R as intermediate_stage_2.
#
# VISCOSITY CAVEAT (read before trusting the eta output)
# ------------------------------------------------------
# The post-UP3 slurry measures eta ~ 12-27 Pa.s at ~40% solid (up3_viscometry.csv)
# -- roughly 60x the UP1-exit value and ~40x what Krieger-Dougherty predicts for
# the concentration jump. The concentrated cake is a structured/flocculated paste
# (likely yield-stress + thixotropic) that the current KD/power-law framework does
# NOT capture from first principles. So this module sets eta from a CALIBRATED
# concentration anchor (mu_ref_cake at Cs_ref), not a mechanistic rheology model.
# The 4 available points do not correlate cleanly with solids (up3_4 at 38.8% is
# the most viscous), so the anchor is a representative value, not a fitted law.
# A mechanistic concentrated-paste closure needs eta measured across a solids
# sweep -- flagged as a data request in review/up3_density_monomer_note.md.
# =============================================================================

up3_separator <- function(stream, pars = list()) {
  p <- modifyList(list(
    Cs_target        = 0.405,   # target cake solids mass frac (measured 38.8-42.7%)
    Fg               = 430,     # centrifugal number [g]; reference up3_1 = 430
    Fg_ref           = 430,     # reference g for the gas/solids closures
    reject_solid_frac= 0.003,   # solid mass lost to centrate (<0.5%; up3_1 = 0.28%)
    alpha_gas_cake   = 0.076,   # trapped-gas floor surviving centrifugation (up3_1)
    rho_s            = 1700,    # bare colloid / dried-product skeletal density
    rho_liq          = 1000,
    # --- calibrated viscosity anchor (NOT first-principles; see header) ---
    mu_ref_cake      = 17.0,    # representative post-UP3 eta at Cs_ref [Pa.s]
    Cs_ref           = 0.405,   # solids the anchor was taken at
    mu_conc_exp      = 8.0      # eta sensitivity to solids about the anchor [-]
  ), pars)

  Cs_in  <- stream$C_solid
  # concentrate (never dilute); mild g-force dependence about the reference
  Cs_out <- max(Cs_in, p$Cs_target * (p$Fg / p$Fg_ref)^0.05)
  Cs_out <- min(Cs_out, 0.60)

  # gas-free condensed (solid+liquid) density at the new concentration
  rho_nogas <- 1 / (Cs_out / p$rho_s + (1 - Cs_out) / p$rho_liq)

  # trapped gas: remove free gas down to the Fg-set cake floor (higher g -> less)
  alpha_floor <- p$alpha_gas_cake * (p$Fg_ref / p$Fg)
  alpha_out   <- min(stream$alpha_g, alpha_floor)

  # calibrated concentration-viscosity anchor (structured paste; see header)
  mu_out <- p$mu_ref_cake * exp(p$mu_conc_exp * (Cs_out - p$Cs_ref))

  stream$C_solid     <- Cs_out * (1 - p$reject_solid_frac)   # reject solid loss
  stream$rho_slurry  <- rho_nogas
  stream$alpha_g     <- alpha_out
  stream$mu_exit_PaS <- mu_out
  stream$P_Pa        <- 1.013e5    # UP3 discharges to atmosphere (gas can escape)
  stream
}
