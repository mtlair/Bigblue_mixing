# =============================================================================
# theta_solvent_chi.R — predict TEMPLATE behavior from chemistry (no experiment)
# =============================================================================
# Maps a template's chemical identity to its process behavior via the
# Flory-Huggins interaction parameter chi, computed from solubility parameters:
#
#     chi = chi_S + (V_s / R T) * (delta_polymer - delta_solvent)^2
#
# with delta = Hildebrand solubility parameter [MPa^0.5], V_s = solvent molar
# volume [cm3/mol], chi_S ~ 0.34 (Flory entropic term). Units reduce cleanly:
# V_s[cm3]*(dDelta[MPa^0.5])^2 = J/mol, /(R T) -> dimensionless.
#
# The theta point is chi = 0.5:
#   chi < 0.5  good solvent  -> swells/dissolves the matrix (miscible)
#   chi = 0.5  theta          -> marginal
#   chi > 0.5  poor solvent  -> phase-separates as discrete droplets (immiscible)
#
# Combined with VOLATILITY (boiling point vs the drying temperature) this predicts
# the template regime and the parameters the process model already consumes
# (template_type, RTF split, T_bp_solv):
#
#   good  + volatile      -> plasticize then COLLAPSE  (dissolves in, escapes, matrix closes) [avoid]
#   good  + non-volatile  -> permanent PLASTICIZER     (stays, softens; not a pore former)
#   theta                 -> marginal / mixed
#   poor  + volatile      -> clean PORE TEMPLATE        (discrete droplets evaporate -> voids)  [best]
#   poor  + non-volatile  -> rigid FILLER / pore former (discrete phase stays)
#   gas                   -> PV=nRT expansion (handled by template_type 2)
#
# NOTE: delta values below are literature Hildebrand values [MPa^0.5] for the
# named species. delta_polymer is the ONE input to set for your system (19.0 is
# a placeholder for a moderately polar vinyl/acrylic matrix). The METHOD is
# chemistry-correct regardless; only the polymer number needs pinning.
# =============================================================================

R_GAS <- 8.314   # J/mol/K

# species | delta [MPa^0.5] | Vm [cm3/mol] | bp [C] | volatile | rho [kg/m3]
solubility_db <- function() data.frame(
  species  = c("polymer_matrix","water","styrene","MMA","toluene","n-hexane",
               "cyclohexane","acetone","ethanol","DBP_plasticizer","paraffin_oil","CO2_liq"),
  delta    = c(19.0, 47.8, 19.0, 18.2, 18.2, 14.9, 16.8, 19.9, 26.5, 19.0, 15.8, 14.6),
  Vm       = c(NA,   18.0, 115,  106,  107,  131,  108,  74,   58.5, 266,  330,  55),
  bp_C     = c(NA,   100,  145,  100,  111,  69,   81,   56,   78,   340,  350,  -57),
  volatile = c(NA,   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE),
  rho      = c(1700, 1000, 906,  940,  867,  655,  779,  784,  789,  1043, 850,  770),
  stringsAsFactors = FALSE)

# --- Flory-Huggins chi from solubility parameters ----------------------------
chi_from_delta <- function(delta_poly, delta_solv, Vm_solv, T_K = 363.15, chi_S = 0.34)
  chi_S + Vm_solv * (delta_poly - delta_solv)^2 / (R_GAS * T_K)

# --- predict template regime + model parameters ------------------------------
# T_dry_C: representative particle/drying temperature that decides "escapes".
template_from_chemistry <- function(template, delta_poly = 19.0,
                                     T_dry_C = 90, db = solubility_db()) {
  p <- db[db$species == "polymer_matrix", ]
  s <- db[db$species == template, ]
  if (nrow(s) == 0) stop("unknown template: ", template)

  gas <- is.finite(s$bp_C) && s$bp_C < 25            # gas at room T
  chi <- if (gas) NA_real_ else
           chi_from_delta(delta_poly, s$delta, s$Vm, T_K = T_dry_C + 273.15)

  solvency <- if (gas) "gas"
              else if (chi < 0.45) "good"
              else if (chi <= 0.55) "theta"
              else if (chi < 1.0)  "poor"
              else                 "very_poor"
  escapes <- isTRUE(s$volatile) && s$bp_C < T_dry_C + 40   # volatile enough to leave in-process

  # RTF (fraction absorbed into particle cores = plasticizing) from miscibility:
  # good solvent -> absorbed (high RTF); poor -> stays as free droplets (low RTF).
  RTF <- if (gas) 0 else 1 / (1 + exp((chi - 0.5) / 0.1))

  regime <- if (gas) "gas_expansion (PV=nRT)"
    else if (solvency %in% c("good","theta") &&  escapes) "plasticize_then_collapse (avoid)"
    else if (solvency %in% c("good","theta") && !escapes) "permanent_plasticizer"
    else if (solvency %in% c("poor","very_poor") &&  escapes) "clean_pore_template (best)"
    else "rigid_filler / pore_former"

  ttype <- if (gas) 2L                                   # gas
    else if (regime == "clean_pore_template (best)") 4L  # capillary_bridge / immiscible emulsion
    else if (regime == "plasticize_then_collapse (avoid)") 3L  # surface_weld / collapse path
    else if (regime == "rigid_filler / pore_former") 1L  # rigid
    else NA_integer_                                      # permanent plasticizer = not a template

  # where the template feeds the process model
  mechanism <- if (gas) "gas holdup -> PV=nRT expansion (up1 chi_npgas, alpha_g)"
    else if (is.na(ttype)) "matrix modifier -> phi_solvent / Softness (NOT a discrete template)"
    else if (ttype == 4L) "free droplet -> phi_templ_free -> D_pore / porosity (up2)"
    else if (ttype == 3L) "core-absorbed then escapes -> w_core + collapse (up2)"
    else "discrete solid -> rigid template load (up1)"

  list(template = template, chi = chi, solvency = solvency, escapes = escapes,
       regime = regime, template_type = ttype, RTF = RTF, mechanism = mechanism,
       T_bp_solv_K = s$bp_C + 273.15, rho_template = s$rho, delta = s$delta)
}

# --- demo / self-test --------------------------------------------------------
if (sys.nframe() == 0) {
  cat("Flory-Huggins chi template screen  (polymer delta = 19.0 MPa^0.5, T_dry = 90 C)\n")
  cat(sprintf("theta point chi = 0.5;  chi<0.5 good/miscible, chi>0.5 poor/immiscible\n\n"))
  cat(sprintf("%-16s %7s %-10s %-8s %-32s %5s %5s\n",
              "template","chi","solvency","escapes","predicted regime","type","RTF"))
  for (t in c("water","styrene","MMA","toluene","n-hexane","cyclohexane",
              "acetone","ethanol","DBP_plasticizer","paraffin_oil","CO2_liq")) {
    r <- template_from_chemistry(t)
    cat(sprintf("%-16s %7s %-10s %-8s %-32s %5s %5.2f\n",
        t, ifelse(is.na(r$chi),"gas",formatC(r$chi,format="f",digits=2)),
        r$solvency, ifelse(r$escapes,"yes","no"), r$regime,
        ifelse(is.na(r$template_type),"-",as.character(r$template_type)), r$RTF))
  }
  cat("\ntemplate_type: 1=rigid 2=gas 3=surface_weld 4=capillary_bridge ('-' = plasticizer, not a template).\n")
  cat("RTF feeds interface_stream.R (core-absorbed plasticizer vs free pore-templating droplet).\n")
  cat("\nMechanism (where each feeds the model):\n")
  for (t in c("MMA","n-hexane","DBP_plasticizer","paraffin_oil","CO2_liq")) {
    r <- template_from_chemistry(t)
    cat(sprintf("  %-16s %s\n", t, r$mechanism))
  }
}
