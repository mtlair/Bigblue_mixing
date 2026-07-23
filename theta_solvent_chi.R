# =============================================================================
# theta_solvent_chi.R — predict TEMPLATE behavior from chemistry (Hansen HSP)
# =============================================================================
# Maps a template's chemical identity to its process behavior via the Hansen
# solubility sphere (RED) and, for informational continuity, the Flory-Huggins
# interaction parameter chi.
#
# HANSEN RELATIVE ENERGY DISTANCE (RED):
#     Ra^2 = 4*(DdD)^2 + (DdP)^2 + (DdH)^2     [HSP components in MPa^0.5]
#     RED  = Ra / r0
#
# Polymer Hansen sphere (pinned from characterization data):
#     dD = 17.00, dP = 12.10, dH = 10.20 MPa^0.5  |  r0 = 8.0
#     Equivalent total Hildebrand: delta = sqrt(dD^2+dP^2+dH^2) = 23.2 MPa^0.5
#
# The theta point is RED = 1:
#   RED < 1   good solvent  -> swells / dissolves the matrix (miscible)
#   RED = 1   theta          -> marginal
#   RED > 1   poor solvent  -> phase-separates as discrete droplets (immiscible)
#
# WHY Hansen, not single Hildebrand:
#   DBP (known plasticizer): Hildebrand chi = 1.16 (wrong - flags as poor)
#                            Hansen RED    = 0.91 (correct - theta/good)
#   MMA:                     Hildebrand chi = 1.30 (poor)
#                            Hansen RED    = 0.95 (theta - more accurate)
#   For a polar polymer (high dP + dH), only solvents with matching polarity
#   and H-bonding are "good". Hildebrand collapses the three components and
#   misclassifies polar species with different dP/dH profiles.
#
# WATER DISPERSIBILITY:
#   In aqueous spray-drying the template must form stable droplets in the slurry
#   (not dissolve into the continuous water phase). Criterion:
#     water_sol < 5 g/100mL  OR  logP > 1.5
#   Templates that fail this dissolve into the aqueous phase and cannot template
#   pores regardless of their RED with the polymer.
#
# REGIME MAP (RED + volatility + water_dispersible -> template_type):
#   good  + volatile      -> plasticize then COLLAPSE  (avoid)     [type 3]
#   good  + non-volatile  -> permanent PLASTICIZER      (softener)  [NA]
#   theta                 -> marginal / mixed                       [type 3]
#   poor  + volatile      -> clean PORE TEMPLATE        (best)      [type 4]
#   poor  + non-volatile  -> rigid FILLER / pore-former             [type 1]
#   gas                   -> PV=nRT expansion                       [type 2]
#
# Volatility threshold: bp_C < T_dry_C + 40 (40 C headroom for evaporation
# kinetics). For spray-drying, use T_dry_C = inlet temperature (~140-160 C)
# rather than outlet (~90 C) to assess whether a template escapes in-process.
# =============================================================================

R_GAS <- 8.314   # J/mol/K

# Polymer Hansen sphere (update r0 if cloud-point or swelling data available;
# user data indicates r0 closer to 8 than 14)
POLY_HSP <- list(dD = 17.00, dP = 12.10, dH = 10.20, r0 = 8.0)

# Database: species | Hansen triplet [MPa^0.5] | Vm [cm3/mol] | bp [C] |
#           volatile | rho [kg/m3] | water_sol [g/100mL, 20 C] | logP (Kow)
# NA in dD/dP/dH for CO2 (gas/supercritical -- RED not defined in same sense).
# NA in water_sol/logP for polymer_matrix (reference entry only).
solubility_db <- function() data.frame(
  species          = c("polymer_matrix",
                       "water", "styrene", "MMA", "toluene",
                       "n-hexane", "cyclohexane", "acetone", "ethanol",
                       "DBP_plasticizer", "paraffin_oil", "CO2_liq",
                       "d-limonene", "butyl_acetate", "chloroform",
                       "methylene_chloride", "triethyl_citrate", "isopropyl_myristate"),
  dD               = c(17.00,
                       15.5,  18.6, 15.8, 18.0,
                       14.9,  16.8, 15.5, 15.8,
                       17.8,  15.8,   NA,
                       17.2,  15.8, 17.8,
                       18.2,  16.5, 14.3),
  dP               = c(12.10,
                       16.0,   1.0,  6.5,  1.4,
                        0.0,   0.0, 10.4,  8.8,
                        8.6,   0.0,   NA,
                        1.8,   3.7,  3.1,
                        6.3,   6.3,  3.7),
  dH               = c(10.20,
                       42.3,   4.1,  5.7,  2.0,
                        0.0,   0.2,  7.0, 19.4,
                        4.0,   0.0,   NA,
                        4.3,   6.3,  5.7,
                        6.1,  12.0,  3.4),
  Vm               = c(  NA,
                        18.0, 115,  106,  107,
                        131,  108,   74,   58.5,
                        266,  330,   55,
                        137,  132,   81,
                         64,  228,  294),
  bp_C             = c(  NA,
                        100,  145,  100,  111,
                         69,   81,   56,   78,
                        340,  350,  -57,
                        176,  126,   61,
                         40,  294,  192),
  volatile         = c(  NA,
                        TRUE,  TRUE, TRUE, TRUE,
                        TRUE,  TRUE, TRUE, TRUE,
                        FALSE, FALSE, TRUE,
                        TRUE,  TRUE, TRUE,
                        TRUE,  FALSE, FALSE),
  rho              = c(1700,
                       1000,  906,  940,  867,
                        655,  779,  784,  789,
                       1043,  850,  770,
                        841,  882, 1479,
                       1325, 1135,  853),
  water_sol_g100mL = c(  NA,
                       100.0, 0.03,  1.5,  0.05,
                       0.001, 0.006, 100.0, 100.0,
                       0.001, 0.001,  NA,
                       0.001,  0.7,   0.8,
                         2.0,  7.0,  0.001),
  logP             = c(  NA,
                       -1.38,  2.95,  1.38,  2.73,
                        3.90,  3.44, -0.24, -0.31,
                        4.72,  7.00,   NA,
                        4.57,  1.78,  1.97,
                        1.25,  0.90,  6.00),
  stringsAsFactors = FALSE)

# --- Hansen RED relative to the polymer sphere --------------------------------
RED_from_hansen <- function(dD_s, dP_s, dH_s, poly = POLY_HSP) {
  Ra2 <- 4*(poly$dD - dD_s)^2 + (poly$dP - dP_s)^2 + (poly$dH - dH_s)^2
  sqrt(Ra2) / poly$r0
}

# --- chi from total Hildebrand (informational; primary classifier is RED) ----
# delta_total = sqrt(dD^2 + dP^2 + dH^2); both polymer and solvent.
chi_from_hansen <- function(dD_s, dP_s, dH_s, Vm_s,
                             poly = POLY_HSP, T_K = 363.15, chi_S = 0.34) {
  d_poly <- sqrt(poly$dD^2 + poly$dP^2 + poly$dH^2)
  d_solv <- sqrt(dD_s^2  + dP_s^2  + dH_s^2)
  chi_S + Vm_s * (d_poly - d_solv)^2 / (R_GAS * T_K)
}

# --- can the template survive as a discrete phase in the aqueous slurry? ------
is_water_dispersible <- function(water_sol, logP)
  !is.na(water_sol) && (water_sol < 5 | (!is.na(logP) && logP > 1.5))

# --- predict template regime + model parameters ------------------------------
# T_dry_C: use INLET temperature (~140-160 C for spray drying) to assess
#          whether a template's bp is low enough to escape in-process.
template_from_chemistry <- function(template, poly = POLY_HSP,
                                     T_dry_C = 90, db = solubility_db()) {
  s <- db[db$species == template, ]
  if (nrow(s) == 0) stop("unknown template: ", template)

  gas <- is.finite(s$bp_C) && s$bp_C < 25   # gas at room temperature

  RED  <- if (gas || anyNA(c(s$dD, s$dP, s$dH))) NA_real_ else
            RED_from_hansen(s$dD, s$dP, s$dH, poly)
  chi  <- if (gas || is.na(s$Vm) || anyNA(c(s$dD, s$dP, s$dH))) NA_real_ else
            chi_from_hansen(s$dD, s$dP, s$dH, s$Vm, poly, T_K = T_dry_C + 273.15)

  solvency <- if (gas) "gas"
    else if (RED <  0.9) "good"
    else if (RED <= 1.1) "theta"
    else if (RED <  1.5) "poor"
    else                 "very_poor"

  escapes <- isTRUE(s$volatile) && is.finite(s$bp_C) && s$bp_C < T_dry_C + 40
  w_disp  <- if (gas) TRUE else is_water_dispersible(s$water_sol_g100mL, s$logP)

  # RTF = core-absorbed fraction; continuous sigmoid on RED (theta RED=1 -> RTF=0.5).
  # good solvent (RED << 1): absorbed into cores -> plasticizes.
  # poor solvent (RED >> 1): stays as free droplets -> pore template.
  RTF <- if (gas) 0 else 1 / (1 + exp((RED - 1) / 0.1))

  regime <- if (gas) "gas_expansion (PV=nRT)"
    else if (solvency %in% c("good","theta") &&  escapes) "plasticize_then_collapse (avoid)"
    else if (solvency %in% c("good","theta") && !escapes) "permanent_plasticizer"
    else if (solvency %in% c("poor","very_poor") &&  escapes) "clean_pore_template (best)"
    else "rigid_filler / pore_former"

  ttype <- if (gas) 2L
    else if (regime == "clean_pore_template (best)") 4L
    else if (regime == "plasticize_then_collapse (avoid)") 3L
    else if (regime == "rigid_filler / pore_former") 1L
    else NA_integer_   # permanent plasticizer -- not a discrete template

  mechanism <- if (gas) "gas holdup -> PV=nRT expansion (chi_npgas, alpha_g)"
    else if (is.na(ttype)) "matrix modifier -> phi_solvent / Softness (NOT a discrete template)"
    else if (ttype == 4L)  "free droplet -> phi_templ_free -> D_pore / porosity (up2)"
    else if (ttype == 3L)  "core-absorbed then escapes -> w_core + collapse (up2)"
    else                   "discrete solid -> rigid template load (up1)"

  list(template          = template,
       RED               = RED,
       chi               = chi,
       solvency          = solvency,
       water_dispersible = w_disp,
       escapes           = escapes,
       regime            = regime,
       template_type     = ttype,
       RTF               = RTF,
       mechanism         = mechanism,
       T_bp_solv_K       = s$bp_C + 273.15,
       rho_template      = s$rho,
       dD = s$dD, dP = s$dP, dH = s$dH,
       logP              = s$logP,
       water_sol_g100mL  = s$water_sol_g100mL)
}

# --- demo / self-test --------------------------------------------------------
if (sys.nframe() == 0) {
  poly   <- POLY_HSP
  d_hild <- sqrt(poly$dD^2 + poly$dP^2 + poly$dH^2)
  cat(sprintf("Polymer HSP: dD=%.2f  dP=%.2f  dH=%.2f  r0=%.1f MPa^0.5\n",
              poly$dD, poly$dP, poly$dH, poly$r0))
  cat(sprintf("Equivalent Hildebrand: %.2f MPa^0.5\n", d_hild))
  cat(sprintf("Theta point: RED=1.0 | T_dry=90 C (outlet; set to ~150 C inlet for spray-dry)\n\n"))

  templates <- c("water","styrene","MMA","toluene","n-hexane","cyclohexane",
                 "acetone","ethanol","DBP_plasticizer","paraffin_oil","CO2_liq",
                 "d-limonene","butyl_acetate","chloroform","methylene_chloride",
                 "triethyl_citrate","isopropyl_myristate")

  cat(sprintf("%-22s %5s %5s %-10s %-6s %-6s %-34s %4s %5s\n",
              "template","RED","chi","solvency","w_dis","escap","regime","type","RTF"))
  cat(strrep("-", 108), "\n")
  for (t in templates) {
    r <- template_from_chemistry(t)
    note <- if (!r$water_dispersible && !isTRUE(r$template_type == 2L))
              " [*]" else ""
    cat(sprintf("%-22s %5s %5s %-10s %-6s %-6s %-34s %4s %5.2f%s\n",
        t,
        ifelse(is.na(r$RED), " gas", formatC(r$RED, format="f", digits=2)),
        ifelse(is.na(r$chi), " gas", formatC(r$chi, format="f", digits=2)),
        r$solvency,
        ifelse(r$water_dispersible, "yes", "NO"),
        ifelse(r$escapes, "yes", "no"),
        r$regime,
        ifelse(is.na(r$template_type), "  -", as.character(r$template_type)),
        r$RTF, note))
  }
  cat("[*] water_sol > 5 g/100mL AND logP < 1.5: dissolves into aqueous phase -- not viable\n")

  cat("\n=== VIABLE CANDIDATES for aqueous spray-drying ===\n")

  cat("\n--- Type 4: clean pore templates (poor solvent + volatile + water-dispersible) ---\n")
  cat("    (Use T_dry_C = spray-dryer inlet ~150 C for border-line bp solvents)\n")
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (isTRUE(r$template_type == 4L) && r$water_dispersible)
      cat(sprintf("  %-22s RED=%4.2f  bp=%3.0f C  logP=%4.1f  RTF=%.3f  water_sol=%.3f g/100mL\n",
                  t, r$RED, r$T_bp_solv_K - 273, r$logP, r$RTF, r$water_sol_g100mL))
  }
  # show d-limonene at inlet T too (borderline at outlet T)
  r_lim_hot <- template_from_chemistry("d-limonene", T_dry_C = 140)
  cat(sprintf("  %-22s RED=%4.2f  bp=%3.0f C  logP=%4.1f  RTF=%.3f  [type %s at T_dry=140 C]\n",
              "d-limonene*", r_lim_hot$RED, r_lim_hot$T_bp_solv_K - 273,
              r_lim_hot$logP, r_lim_hot$RTF,
              ifelse(is.na(r_lim_hot$template_type), "-", as.character(r_lim_hot$template_type))))

  cat("\n--- Type 1: rigid fillers / pore-formers (poor + non-volatile + water-dispersible) ---\n")
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (isTRUE(r$template_type == 1L) && r$water_dispersible)
      cat(sprintf("  %-22s RED=%4.2f  bp=%3.0f C  logP=%4.1f  RTF=%.3f\n",
                  t, r$RED, r$T_bp_solv_K - 273, r$logP, r$RTF))
  }

  cat("\n--- Permanent plasticizers (type NA: good/theta + non-volatile + water-dispersible) ---\n")
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (is.na(r$template_type) && !r$escapes && r$water_dispersible)
      cat(sprintf("  %-22s RED=%4.2f  logP=%4.1f  RTF=%.3f  water_sol=%.4f g/100mL\n",
                  t, r$RED, r$logP, r$RTF, r$water_sol_g100mL))
  }

  cat("\n--- Type 2: gas template ---\n  CO2_liq: gas expansion (PV=nRT)\n")

  cat("\n--- Hansen vs Hildebrand: key discrepancies for this polymer ---\n")
  for (t in c("DBP_plasticizer","MMA","toluene","butyl_acetate")) {
    r <- template_from_chemistry(t)
    cat(sprintf("  %-22s RED=%4.2f (%s)  chi(Hild)=%4.2f  -> %s\n",
                t, r$RED, r$solvency,
                ifelse(is.na(r$chi), NA, r$chi),
                r$regime))
  }

  cat("\ntemplate_type: 1=rigid  2=gas  3=surface_weld  4=capillary_bridge  -=plasticizer\n")
  cat("RED: Hansen relative energy distance (RED=1 = theta, r0=8.0)\n")
  cat("RTF: core-absorbed fraction [feeds interface_stream.R w_core vs phi_templ_free]\n")
  cat("chi: Flory-Huggins from Hildebrand equivalent (informational; RED used for regime)\n")
}
