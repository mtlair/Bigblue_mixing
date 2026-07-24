# =============================================================================
# theta_solvent_chi.R — predict TEMPLATE behavior from chemistry (Hansen HSP)
# =============================================================================
# Maps a template's chemical identity to its process behavior via the Hansen
# solubility sphere (RED) and Flory-Huggins interaction parameter chi.
#
# HANSEN RELATIVE ENERGY DISTANCE (RED):
#     Ra^2 = 4*(DdD)^2 + (DdP)^2 + (DdH)^2     [HSP components in MPa^0.5]
#     RED  = Ra / r0
#
# Polymer Hansen sphere (pinned from characterization data):
#     dD = 17.00, dP = 12.10, dH = 10.20 MPa^0.5  |  r0 = 8.0
#     Equivalent total Hildebrand: delta = sqrt(dD^2+dP^2+dH^2) = 23.2 MPa^0.5
#
# NOTE: The uploaded HSP_Calculations_Ruben-Manuel.xlsx contains a pre-fitted
# 1-sphere result (dD=17.83, dP=9.76, dH=8.07, r0=7.07) and a 2-sphere result
# (Sphere2: dD=17.76, dP=11.42, dH=8.17, r0=5.26) from solvent dissolution
# tests. Sphere 2 is close to the user-supplied values above. Both are stored
# as alternatives in POLY_HSP_ALT below; the primary POLY_HSP uses the
# user-supplied values.
#
# WATER DISPERSIBILITY:
#   Template must form stable droplets in the aqueous slurry (not dissolve into
#   the continuous phase). Criterion: water_sol < 5 g/100mL OR logP > 1.5.
#
# REGIME MAP:
#   good  + volatile      -> plasticize then COLLAPSE  (avoid)     [type 3]
#   good  + non-volatile  -> permanent PLASTICIZER      (softener)  [NA]
#   theta                 -> marginal / mixed                       [type 3]
#   poor  + volatile      -> clean PORE TEMPLATE        (best)      [type 4]
#   poor  + non-volatile  -> rigid FILLER / pore-former             [type 1]
#   gas                   -> PV=nRT expansion                       [type 2]
#
# VOLATILITY THRESHOLD: bp_C < T_dry_C + 40.
# Use T_dry_C = dryer inlet (~140-160 C) to classify border-line solvents.
# =============================================================================

R_GAS <- 8.314   # J/mol/K

# Primary polymer Hansen sphere (user-supplied, confirmed close to xlsx fit)
POLY_HSP <- list(dD = 17.00, dP = 12.10, dH = 10.20, r0 = 8.0)

# Alternative fitted spheres from HSP_Calculations_Ruben-Manuel.xlsx:
#   1-sphere: dD=17.83, dP=9.76, dH=8.07, r0=7.07
#   2-sphere Sphere2: dD=17.76, dP=11.42, dH=8.17, r0=5.26 (closer to primary)
POLY_HSP_1SPH <- list(dD = 17.83, dP =  9.76, dH = 8.07, r0 = 7.07)
POLY_HSP_2SPH <- list(dD = 17.76, dP = 11.42, dH = 8.17, r0 = 5.26)

# Database: species | Hansen triplet [MPa^0.5] | Vm [cm3/mol] | bp [C] |
#           volatile | rho [kg/m3] | water_sol [g/100mL, 20 C] | logP
# Vm sourced from literature; NA if unavailable.
solubility_db <- function() data.frame(
  species          = c("polymer_matrix",
                       # original set
                       "water", "styrene", "MMA", "toluene",
                       "n-hexane", "cyclohexane", "acetone", "ethanol",
                       "DBP_plasticizer", "paraffin_oil", "CO2_liq",
                       "d-limonene", "butyl_acetate", "chloroform",
                       "methylene_chloride", "triethyl_citrate", "isopropyl_myristate",
                       # expanded set (from HSP database screen)
                       "ethyl_butyrate", "isoamyl_acetate", "amyl_acetate",
                       "hexyl_acetate", "dibutyl_sebacate",
                       "chlorobenzene", "xylene",
                       "ethyl_levulinate", "butyl_levulinate",
                       # additional candidates (user request)
                       "butyl_butyrate", "octyl_acetate"),
  dD               = c(17.00,
                       15.5, 18.6, 15.8, 18.0, 14.9, 16.8, 15.5, 15.8,
                       17.8, 15.8,   NA, 17.2, 15.8, 17.8, 18.2, 16.5, 14.3,
                       15.5, 15.3, 15.8, 15.8, 16.7, 19.0, 17.8, 14.6, 15.7,
                       15.6, 15.8),
  dP               = c(12.10,
                       16.0,  1.0,  6.5,  1.4,  0.0,  0.0, 10.4,  8.8,
                        8.6,  0.0,   NA,  1.8,  3.7,  3.1,  6.3,  6.3,  3.7,
                        5.6,  3.1,  3.3,  2.9,  4.5,  4.3,  1.0, 10.5,  9.7,
                        2.9,  2.9),
  dH               = c(10.20,
                       42.3,  4.1,  5.7,  2.0,  0.0,  0.2,  7.0, 19.4,
                        4.0,  0.0,   NA,  4.3,  6.3,  5.7,  6.1, 12.0,  3.4,
                        5.0,  7.0,  6.1,  5.9,  4.1,  2.0,  3.1,  7.0,  5.8,
                        5.6,  5.1),
  Vm               = c(  NA,
                        18.0, 115, 106, 107, 131, 108,  74,  58.5,
                        266,  330,  55, 137, 132,  81,  64, 228, 294,
                        132,  149, 148, 168, 340,  102, 124, 130, 158,
                        165, 199),
  bp_C             = c(  NA,
                        100, 145, 100, 111,  69,  81,  56,  78,
                        340, 350, -57, 176, 126,  61,  40, 294, 192,
                        121, 142, 149, 171, 345,  132, 139, 206, 238,
                        166, 210),
  volatile         = c(  NA,
                       TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
                       FALSE,FALSE,TRUE, TRUE, TRUE, TRUE, TRUE,FALSE,FALSE,
                       TRUE, TRUE, TRUE, TRUE,FALSE, TRUE, TRUE,FALSE,FALSE,
                       FALSE,FALSE),
  rho              = c(1700,
                       1000, 906, 940, 867, 655, 779, 784, 789,
                       1043, 850, 770, 841, 882,1479,1325,1135, 853,
                        879, 872, 879, 870, 936,1106, 864, 994, 970,
                        875, 865),
  water_sol_g100mL = c(  NA,
                       100.0, 0.03, 1.5, 0.05, 0.001, 0.006, 100.0, 100.0,
                       0.001, 0.001,  NA, 0.001,  0.7,   0.8,   2.0,   7.0, 0.001,
                        0.55,  0.20,  0.17, 0.02, 0.001,  0.05, 0.014,   8.0,  1.5,
                        0.04,  0.01),
  logP             = c(  NA,
                       -1.38, 2.95, 1.38, 2.73,  3.90,  3.44, -0.24, -0.31,
                        4.72,  7.00,  NA,  4.57,  1.78,  1.97,  1.25,  0.90,  6.00,
                        1.68,  2.27,  2.30,  2.83,  5.27,  2.84,  3.15,  0.90,  1.80,
                        2.83,  3.72),
  stringsAsFactors = FALSE)

# --- Hansen RED relative to the polymer sphere --------------------------------
RED_from_hansen <- function(dD_s, dP_s, dH_s, poly = POLY_HSP) {
  Ra2 <- 4*(poly$dD - dD_s)^2 + (poly$dP - dP_s)^2 + (poly$dH - dH_s)^2
  sqrt(Ra2) / poly$r0
}

# --- chi from total Hildebrand (informational; primary classifier is RED) ----
chi_from_hansen <- function(dD_s, dP_s, dH_s, Vm_s,
                             poly = POLY_HSP, T_K = 363.15, chi_S = 0.34) {
  d_poly <- sqrt(poly$dD^2 + poly$dP^2 + poly$dH^2)
  d_solv <- sqrt(dD_s^2  + dP_s^2  + dH_s^2)
  chi_S + Vm_s * (d_poly - d_solv)^2 / (R_GAS * T_K)
}

# --- water dispersibility check ----------------------------------------------
is_water_dispersible <- function(water_sol, logP)
  !is.na(water_sol) && (water_sol < 5 | (!is.na(logP) && logP > 1.5))

# --- predict template regime + model parameters ------------------------------
# T_dry_C: use INLET temperature (~140-160 C for atomizer drying) for volatility.
template_from_chemistry <- function(template, poly = POLY_HSP,
                                     T_dry_C = 90, db = solubility_db()) {
  s <- db[db$species == template, ]
  if (nrow(s) == 0) stop("unknown template: ", template)

  gas <- is.finite(s$bp_C) && s$bp_C < 25

  RED  <- if (gas || anyNA(c(s$dD, s$dP, s$dH))) NA_real_ else
            RED_from_hansen(s$dD, s$dP, s$dH, poly)
  chi  <- if (gas || is.na(s$Vm) || anyNA(c(s$dD, s$dP, s$dH))) NA_real_ else
            chi_from_hansen(s$dD, s$dP, s$dH, s$Vm, poly, T_K = T_dry_C + 273.15)

  solvency <- if (gas) "gas"
    else if (RED <  0.9) "good"
    else if (RED <= 1.1) "theta"
    else if (RED <  1.5) "poor"
    else                 "very_poor"

  # escapes: reconciled with the up2 dryer's constant-rate + falling-rate physics.
  # volatile=FALSE only used for species where evaporation is structurally
  # impossible (DBP, paraffin, IPM).
  truly_nonvolatile <- identical(s$volatile, FALSE) && (is.na(s$bp_C) || s$bp_C > 280)
  # The constant-rate droplet surface sits near the WET-BULB (~adiabatic
  # saturation, ~45 C for aqueous drying) -- NOT the hot inlet gas. The old
  # `bp < T_dry+40` boolean used the inlet T and over-promised ~10x for
  # high-boiling solvents (e.g. d-limonene bp 176: old flag TRUE, but up2 gives
  # f_esc = 0.06). Here we score the reduced vapour pressure at the wet-bulb
  # (constant-rate escape) OR a falling-rate boil check (particle reaches
  # ~T_dry-40), mirroring up2's f_cr + S_bs. up2's f_esc remains authoritative;
  # this is a screening estimate.
  T_WB_C <- 45                                         # constant-rate wet-bulb [C]
  p_r_wb <- if (is.finite(s$bp_C)) {
              T_bp_K <- s$bp_C + 273.15
              min(1, exp(88 * T_bp_K / R_GAS * (1/T_bp_K - 1/(T_WB_C + 273.15))))
            } else 0
  boils  <- is.finite(s$bp_C) && s$bp_C < (T_dry_C - 40)     # falling-rate boiling
  escape_index <- max(p_r_wb, if (boils) 1 else 0)
  escape_class <- if (escape_index >= 0.10) "clean"
                  else if (escape_index >= 0.05) "partial" else "retained"
  escapes <- !truly_nonvolatile && (escape_index >= 0.10)
  w_disp  <- if (gas) TRUE else is_water_dispersible(s$water_sol_g100mL, s$logP)

  # RTF = core-absorbed fraction; continuous sigmoid on RED (RED=1 -> RTF=0.5)
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
    else NA_integer_

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
       escape_index      = escape_index,
       escape_class      = escape_class,
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

# --- screen the Abbott HSP database ------------------------------------------
# Reads HSP_Calculations_Ruben-Manuel.xlsx and ranks all 1269 chemicals by RED.
# Requires: readxl (install.packages("readxl")) or the xlsx path readable by
# python3 (called via system2 if readxl absent). Falls back gracefully.
#
# water_immiscible_names: character vector of known water-immiscible species
# (case-insensitive substring match) to flag as w_disp = TRUE in the output.
screen_hsp_database <- function(
    xlsx_path = "references/HSP_Calculations_Ruben-Manuel.xlsx",
    poly      = POLY_HSP,
    RED_min   = 1.0,
    RED_max   = 2.5,
    top_n     = 50,
    water_immiscible_names = c(
      "acetate","butyrate","hexanoate","caproate","caprylate","laurate",
      "palmitate","stearate","linoleate","oleate","sebacate","phthalate",
      "toluene","benzene","xylene","styrene","limonene","cyclohexane",
      "hexane","heptane","octane","decane","dodecane","paraffin","isooctane",
      "chlorobenzene","chloroform","dichloromethane","methylene chloride",
      "coconut","mineral","isopar","dodecanol","lauryl","cetyl","stearyl"
    )) {

  # try readxl first
  db <- tryCatch({
    if (!requireNamespace("readxl", quietly=TRUE)) stop("no readxl")
    df <- readxl::read_excel(xlsx_path, sheet="Chemicals")
    df <- df[, 1:6]; names(df) <- c("CAS","Chemical","dD","dP","dH","Vm")
    df
  }, error = function(e) {
    # fallback: python3 to parse
    tmp <- tempfile(fileext=".csv")
    cmd <- sprintf(
      "python3 -c \"
import openpyxl, warnings, csv
warnings.filterwarnings('ignore')
wb=openpyxl.load_workbook('%s',data_only=True)
rows=list(wb['Chemicals'].iter_rows(values_only=True))
import sys; w=csv.writer(sys.stdout)
w.writerow(['CAS','Chemical','dD','dP','dH','Vm'])
for r in rows[1:]:
  if all(isinstance(r[i],(int,float)) for i in [2,3,4]):
    w.writerow([r[0] or '',r[1] or '',r[2],r[3],r[4],r[5] or ''])
\" > %s", xlsx_path, tmp)
    system(cmd)
    df <- read.csv(tmp, stringsAsFactors=FALSE)
    df$dD <- as.numeric(df$dD); df$dP <- as.numeric(df$dP)
    df$dH <- as.numeric(df$dH); df$Vm  <- suppressWarnings(as.numeric(df$Vm))
    df
  })

  db <- db[complete.cases(db[,c("dD","dP","dH")]), ]
  db$RED <- mapply(function(dD,dP,dH) RED_from_hansen(dD,dP,dH,poly),
                   db$dD, db$dP, db$dH)
  db <- db[db$RED >= RED_min & db$RED <= RED_max, ]
  db <- db[order(db$RED), ]

  # flag known water-immiscible species
  lname <- tolower(db$Chemical)
  db$w_disp_flag <- sapply(lname, function(n)
    any(sapply(water_immiscible_names, function(p) grepl(p, n, fixed=TRUE))))

  head(db[, c("Chemical","dD","dP","dH","Vm","RED","w_disp_flag")], top_n)
}

# --- NxtLevel green solvent summary ------------------------------------------
nxtlevel_summary <- function(poly = POLY_HSP) {
  nxt <- data.frame(
    name = c("Ethyl levulinate","Butyl levulinate","EtLPK","EtLGK","EHL"),
    dD   = c(14.60, 15.69, 15.17, 15.27, 15.27),
    dP   = c(10.46,  9.66,  4.68, 12.72,  6.43),
    dH   = c( 6.98,  5.82,  6.80, 10.93,  3.85),
    # literature / estimated bp; EtLPK/EHL/EtLGK: confirm with supplier
    bp_C      = c(206,  238,   NA,    NA,    NA),
    water_sol = c(8.0,  1.5,   NA,    NA,    NA),
    logP      = c(0.90, 1.80,  NA,    NA,    NA),
    note      = c("partially water-soluble","limited water-sol",
                  "bp unknown-confirm supplier","bp unknown-confirm supplier",
                  "bp unknown-confirm supplier"),
    stringsAsFactors=FALSE)
  nxt$RED <- mapply(function(dD,dP,dH) RED_from_hansen(dD,dP,dH,poly),
                    nxt$dD, nxt$dP, nxt$dH)
  nxt$w_disp <- mapply(is_water_dispersible, nxt$water_sol, nxt$logP)
  nxt$regime <- ifelse(nxt$RED < 0.9, "good->plasticizer",
                ifelse(nxt$RED <= 1.1, "theta/marginal",
                ifelse(nxt$RED < 1.5, "poor->pore-template", "very_poor")))
  nxt[, c("name","dD","dP","dH","RED","regime","bp_C","w_disp","note")]
}

# --- demo / self-test --------------------------------------------------------
if (sys.nframe() == 0) {
  poly   <- POLY_HSP
  d_hild <- sqrt(poly$dD^2 + poly$dP^2 + poly$dH^2)
  cat(sprintf("Polymer HSP: dD=%.2f  dP=%.2f  dH=%.2f  r0=%.1f MPa^0.5\n",
              poly$dD, poly$dP, poly$dH, poly$r0))
  cat(sprintf("Equiv Hildebrand: %.2f MPa^0.5 | theta: RED=1.0\n", d_hild))
  cat(sprintf("Xlsx-fitted sphere (1-sph): dD=%.2f dP=%.2f dH=%.2f r0=%.2f\n",
              POLY_HSP_1SPH$dD, POLY_HSP_1SPH$dP, POLY_HSP_1SPH$dH, POLY_HSP_1SPH$r0))
  cat(sprintf("Xlsx-fitted sphere (2-sph2): dD=%.2f dP=%.2f dH=%.2f r0=%.2f\n\n",
              POLY_HSP_2SPH$dD, POLY_HSP_2SPH$dP, POLY_HSP_2SPH$dH, POLY_HSP_2SPH$r0))

  db <- solubility_db()
  templates <- db$species[db$species != "polymer_matrix"]

  cat(sprintf("%-24s %5s %5s %-10s %-6s %-6s %-34s %4s %5s\n",
              "template","RED","chi","solvency","w_dis","escap","regime","type","RTF"))
  cat(strrep("-", 112), "\n")
  for (t in templates) {
    r    <- template_from_chemistry(t)
    note <- if (!r$water_dispersible && !isTRUE(r$template_type == 2L)) " [*]" else ""
    cat(sprintf("%-24s %5s %5s %-10s %-6s %-6s %-34s %4s %5.2f%s\n",
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
  cat("[*] water_sol > 5 g/100mL AND logP < 1.5 -> dissolves in aqueous phase\n")

  cat("\n=== VIABLE CANDIDATES (water-dispersible, aqueous atomizer-drying) ===\n")
  cat("\n--- Type 4: clean pore templates (poor + volatile + w-dispersible) ---\n")
  cat("    Tip: use T_dry_C = inlet temp (~150 C) for high-bp solvents\n")
  cat(sprintf("  %-24s %5s %5s %5s %5s %5s\n","template","RED","bp C","logP","RTF","Vmol"))
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (isTRUE(r$template_type == 4L) && r$water_dispersible) {
      s <- db[db$species == t, ]
      cat(sprintf("  %-24s %5.2f %5.0f %5.2f %5.3f %5s\n",
                  t, r$RED, r$T_bp_solv_K-273, r$logP,
                  r$RTF, ifelse(is.na(s$Vm),"-",as.character(s$Vm))))
    }
  }
  # d-limonene at inlet T
  r_lim <- template_from_chemistry("d-limonene", T_dry_C=140)
  if (isTRUE(r_lim$template_type == 4L))
    cat(sprintf("  %-24s %5.2f %5.0f %5.2f %5.3f   137  [T_dry=140 C]\n",
                "d-limonene*", r_lim$RED, 176, r_lim$logP, r_lim$RTF))
  # amyl/hexyl/butyl_butyrate at inlet T
  for (t in c("amyl_acetate","hexyl_acetate","butyl_butyrate")) {
    r2 <- template_from_chemistry(t, T_dry_C=140)
    if (isTRUE(r2$template_type==4L)) {
      s <- db[db$species==t,]
      cat(sprintf("  %-24s %5.2f %5.0f %5.2f %5.3f %5s  [T_dry=140 C]\n",
                  paste0(t,"*"), r2$RED, r2$T_bp_solv_K-273, r2$logP, r2$RTF,
                  ifelse(is.na(s$Vm),"-",as.character(s$Vm))))
    }
  }
  # octyl_acetate needs very high inlet T (bp=210); show for reference
  r_oct <- template_from_chemistry("octyl_acetate", T_dry_C=170)
  cat(sprintf("  %-24s %5.2f %5.0f %5.2f %5.3f %5s  [T_dry=170 C only]\n",
              "octyl_acetate*", r_oct$RED, 210, r_oct$logP, r_oct$RTF, "199"))

  cat("\n--- Permanent plasticizers (good/theta + non-volatile + w-dispersible) ---\n")
  cat(sprintf("  %-24s %5s %5s %5s\n","template","RED","logP","RTF"))
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (is.na(r$template_type) && !r$escapes && r$water_dispersible)
      cat(sprintf("  %-24s %5.2f %5.2f %5.3f\n", t, r$RED, r$logP, r$RTF))
  }

  cat("\n--- Type 1: rigid fillers (poor + non-volatile + w-dispersible) ---\n")
  for (t in templates) {
    r <- template_from_chemistry(t)
    if (isTRUE(r$template_type==1L) && r$water_dispersible)
      cat(sprintf("  %-24s RED=%4.2f  bp=%3.0f C  logP=%4.1f\n",
                  t, r$RED, r$T_bp_solv_K-273, r$logP))
  }

  cat("\n=== NXTLEVEL GREEN SOLVENTS (from uploaded xlsx) ===\n")
  nxt <- nxtlevel_summary()
  print(nxt, row.names=FALSE)
  cat("EtLPK / EHL / EtLGK: bp unknown -- confirm boiling point with NxtLevel\n")
  cat("EtLPK (RED=1.12) and EHL (RED=1.15) are in the pore-template range\n")
  cat("if bp < T_dry+40 they classify as type 4; if non-volatile -> type 1\n")

  cat("\n=== HSP DATABASE SCREEN (top 30, RED 1.1-2.0, water-immiscible solvents) ===\n")
  cat("Screening all 1269 chemicals in HSP_Calculations_Ruben-Manuel.xlsx ...\n")
  scr <- tryCatch(
    screen_hsp_database(RED_min=1.1, RED_max=2.0, top_n=30),
    error=function(e) { cat("  [screen failed:", conditionMessage(e), "]\n"); NULL })
  if (!is.null(scr)) {
    cat(sprintf("  %-45s %5s %5s %5s %5s  %5s  %s\n",
                "Chemical","dD","dP","dH","Vm","RED","w_disp?"))
    for (i in seq_len(nrow(scr))) {
      r <- scr[i,]
      cat(sprintf("  %-45s %5.1f %5.1f %5.1f %5s  %5.3f  %s\n",
                  r$Chemical, r$dD, r$dP, r$dH,
                  ifelse(is.na(r$Vm),"-",as.character(round(r$Vm))),
                  r$RED, ifelse(r$w_disp_flag,"likely","check")))
    }
  }

  cat("\ntemplate_type: 1=rigid  2=gas  3=surface_weld  4=capillary_bridge  -=plasticizer\n")
  cat("RED: Hansen relative energy distance (RED=1 = theta, r0=8.0)\n")
  cat("RTF: core-absorbed fraction [feeds interface_stream.R]\n")
}
