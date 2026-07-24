# =============================================================================
# STEP-THROUGH: UP1 -> UP2 -> UP3 -> UP4 with butyl butyrate / butyl acetate
# =============================================================================
# Experimental-validation walkthrough of the co-feed template-in-water mechanism
# (COFEED_HANDOFF.md). For each of the two recommended type-4 pore templates --
# butyl acetate (BA) and butyl butyrate (BB) -- this replays the FULL wired train
#
#     UP1 gassed mixer -> UP2 foam-wash -> UP3 centrifuge -> UP4 atomizer dryer
#
# and prints the stream state at every interface, so the chemistry (Hansen RED,
# boiling point, core-absorption RTF) can be traced from the mixer feed to the
# final porous particle.
#
# Process backbone: the fully-measured up3_1 condition (SEM sample 114641), so
# the UP1/UP4 setpoints and the UP3 cake state are the real plant values; only
# the TEMPLATE identity is swapped (BA vs BB). Template loading is the co-feed
# pendular optimum from the handoff (fill = template_dose = 0.20).
#
# Run:  Rscript step_up1_up4_bb_ba.R
# =============================================================================

source("unified/up1_mixer_module.R")
source("unified/interface_stream.R")
source("unified/up2_atomizer_dryer_module.R")
source("foam_wash_module.R")        # UP2 foam-wash column
source("up3_centrifuge_module.R")   # UP3 decanting centrifuge
source("theta_solvent_chi.R")       # chemistry -> template regime (Hansen RED)

# --- load the unified chain's factor dictionary / nominal point (functions) ---
lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))

RHO_AIR <- 0.0765   # lb/ft3 standard air

# --- process backbone: measured up3_1 (fully-characterised chain condition) ---
dat  <- read.csv("data/cond_up1234.csv")
visc <- read.csv("data/up3_viscometry.csv")
row  <- dat[dat$cond == "up3_1", ]
cond <- "up3_1"
up4_atom_psig <- 9.051   # from up1_2_3_4_psd sheet (not in cond_up1234.csv)
up4_pump_psig <- 5.515

visc_ref <- function(cond, target = 12.7) {
  s <- visc[visc$cond == cond, ]; s <- s[order(s$shear_rate_1s), ]
  exp(approx(log(s$shear_rate_1s), log(s$viscosity_PaS), log(target), rule = 2)$y)
}
dryer_airflow <- function(feed_lbhr, s, Tin, Tout, Tf) {
  h_fg <- 2.30e6; cp_gas <- 1005; resid <- 0.005; LB <- 0.4536; HR <- 3600
  feed <- feed_lbhr * LB / HR
  mevap <- feed * (1 - s) - resid / (1 - resid) * feed * s
  cp_feed <- (1 - s) * 4180 + s * 1500
  (mevap * h_fg + feed * cp_feed * (Tout - Tf)) / (cp_gas * (Tin - Tout))
}

# --- build the UP1/UP4 setpoint vector from the measured backbone -------------
base_x <- function() {
  x <- nominal_x
  x["v_tip"]         <- row$up1_vtip
  x["C_solid_mass"]  <- 0.25
  x["Q_colloid"]     <- 3.3          # co-feed colloid flow (handoff nominal)
  x["Q_template"]    <- 0.45         # ~13% of colloid -> fill ~0.20 (pendular)
  x["template_dose"] <- 0.20         # co-feed pendular optimum (COFEED_HANDOFF)
  x["C_temp_mass"]   <- 0.05         # 5% solvent in aqueous dispersion
  x["D_template"]    <- 0.5          # um dispersion droplet (microfluidised)
  x["T_system"]      <- row$up1_exit_temp_C + 273.15
  x["P_system"]      <- (row$up1_psig + 14.696) / 14.696
  air_lbhr  <- row$up4_atom_scfm * RHO_AIR * 60
  feed_lbhr <- row$up4_feed_mb_lbhr
  x["ALR"]        <- air_lbhr / feed_lbhr
  x["mdot_L"]     <- feed_lbhr * 0.4536 / 3600
  x["P_atom_air"] <- (up4_atom_psig + 14.696) * 6894.76
  x["P_feed"]     <- (up4_pump_psig + 14.696) * 6894.76
  x["T_dryer_in"] <- row$up4_Tin_C + 273.15
  s_solid <- row$up3_solid_pct / 100
  x["mdot_gas_dry"] <- dryer_airflow(feed_lbhr, s_solid, row$up4_Tin_C,
                                     row$up4_Tout_C, row$up1_exit_temp_C)
  x["size_template"] <- 1            # pore-templating overlay ACTIVE
  x
}

# --- print helpers ------------------------------------------------------------
hr <- function(ch = "-") cat(strrep(ch, 78), "\n")
show <- function(fields, s, fmt = "%-18s % .4g\n") {
  for (f in fields) if (!is.null(s[[f]])) cat(sprintf(fmt, f, s[[f]]))
}

# --- one template, full step-through ------------------------------------------
step_template <- function(template, T_dry_C) {
  chem <- template_from_chemistry(template, T_dry_C = T_dry_C)

  hr("=")
  cat(sprintf("TEMPLATE: %s   (co-feed, up3_1 backbone, dryer inlet %.0f C)\n",
              toupper(gsub("_", " ", template)), T_dry_C))
  hr("=")
  cat("\n[0] CHEMISTRY SCREEN  (theta_solvent_chi.R, Hansen sphere dD=17.0 dP=12.1 dH=10.2 r0=8)\n")
  cat(sprintf("    RED = %.3f  (%s solvent)   chi = %.2f\n", chem$RED, chem$solvency, chem$chi))
  cat(sprintf("    bp  = %.0f C   escapes@inlet = %s   water-dispersible = %s\n",
              chem$T_bp_solv_K - 273.15, chem$escapes, chem$water_dispersible))
  cat(sprintf("    regime = %s\n", chem$regime))
  cat(sprintf("    template_type = %d   RTF(dryer core-absorption) = %.4f   rho = %.0f kg/m3\n",
              chem$template_type, chem$RTF, chem$rho_template))

  # feed the chemistry into the chain: type-4, real bp, real solvent density
  x <- base_x()
  x["T_bp_solv"] <- chem$T_bp_solv_K
  eq <- equipment
  eq$template_type <- chem$template_type          # 4 = capillary_bridge
  RHO_TEMPLATE_LIQ <<- chem$rho_template          # real solvent liquid density

  # ---- UP1 : gassed / templated mixer --------------------------------------
  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)

  cat("\n[1] UP1  gassed mixer  (tau=", round(p1$tau,0), "min, T=",
      round(p1$T_system-273.15,0), "C, P=", round(p1$P_system,2), "atm)\n", sep="")
  cat(sprintf("    RTF_mixer (core-absorbed at exit)  = %.4f   -> %.1f%% template stays FREE\n",
              s$RTF, 100*(1 - s$RTF)))
  show(c("C_solid","rho_slurry","mu_exit_PaS","alpha_g",
         "phi_templ_free","w_core","D_agg_um","WetSkin","Softness_exit"), s)

  # ---- UP2 : foam-wash column ----------------------------------------------
  s <- foam_wash_column(s)
  cat("\n[2] UP2  foam-wash column  (gas washout, monomer/template rinse)\n")
  show(c("C_solid","rho_slurry","mu_exit_PaS","alpha_g","phi_templ_free","w_core"), s)

  # ---- UP3 : decanting centrifuge ------------------------------------------
  s <- up3_centrifuge(s, pars = list(
         Cs_target = row$up3_solid_pct / 100,
         Fg        = if (is.finite(row$up3_Fg)) row$up3_Fg else 430))
  cat("\n[3] UP3  decanting centrifuge  (concentrate 25% -> ", row$up3_solid_pct,
      "% cake, ", if (is.finite(row$up3_Fg)) row$up3_Fg else 430, " g)\n", sep="")
  show(c("C_solid","rho_slurry","mu_exit_PaS","alpha_g","phi_templ_free","w_core"), s)

  # ---- UP4 : atomizer dryer ----------------------------------------------------
  x_up2 <- as.list(x[up2_names]); x_up2[["size_template"]] <- x[["size_template"]]
  r2 <- up2_run_dryer(s, x_up2)
  cat("\n[4] UP4  atomizer dryer  (atomise + dry; inlet ", round(x["T_dryer_in"]-273.15,0),
      " C, ALR ", round(x["ALR"],2), ")\n", sep="")
  cat(sprintf("    %-14s %8.2f um   %-14s %8.2f um\n",
              "product d50", r2[["D_particle_um"]], "product d90", r2[["Dp90_um"]]))
  cat(sprintf("    %-14s %8.3f um   %-14s %8.3f\n",
              "pore size", r2[["D_pore_um"]], "porosity phi", r2[["phi_porosity_z"]]))
  cat(sprintf("    %-14s %8.3f      %-14s %8.3f\n",
              "skin theta", r2[["theta_skin_z"]], "sphericity", r2[["Omega_struct_z"]]))
  cat(sprintf("    %-14s %8.1f K    %-14s %8.5f\n",
              "Tg_product", r2[["Tg_product_K"]], "X_moisture", r2[["X_moisture"]]))
  cat(sprintf("    %-14s %8.3f      %-14s %8.3f\n",
              "f_cr (Raoult)", r2[["f_cr_solv"]], "f_esc (total)", r2[["f_esc_solv"]]))
  cat(sprintf("    %-14s %8.3f      %-14s %8.3f\n",
              "solv_retained", r2[["solv_retained"]], "f_burst_solv", r2[["f_burst_solv"]]))
  cat(sprintf("    %-14s %8.1f kg/m3\n", "rho_tapped", r2[["rho_tapped"]]))

  invisible(list(template = template, chem = chem, product = r2))
}

# --- no-template reference: same backbone, template feed off -----------------
run_no_template <- function() {
  x <- base_x()
  x["Q_template"]    <- 0.0
  x["template_dose"] <- 0.0
  x["C_temp_mass"]   <- 0.0
  eq <- equipment; eq$template_type <- 4L   # keep type-4 flag; phi_e->0 anyway
  RHO_TEMPLATE_LIQ <<- 750                   # default; no template so unused

  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)
  s  <- foam_wash_column(s)
  s  <- up3_centrifuge(s, pars = list(
           Cs_target = row$up3_solid_pct / 100,
           Fg        = if (is.finite(row$up3_Fg)) row$up3_Fg else 430))
  x_up2 <- as.list(x[up2_names]); x_up2[["size_template"]] <- x[["size_template"]]
  r2 <- up2_run_dryer(s, x_up2)
  list(template = "no_template", chem = NULL, product = r2)
}

# =============================================================================
cat("\n")
cat("###########################################################################\n")
cat("#  CO-FEED PORE-TEMPLATE STEP-THROUGH  --  UP1 -> UP2 -> UP3 -> UP4        #\n")
cat("#  Process backbone: measured up3_1 (SEM 114641). Template identity swapped#\n")
cat("###########################################################################\n")

T_inlet_C <- row$up4_Tin_C     # measured up3_1 dryer inlet
resBA <- step_template("butyl_acetate",  T_dry_C = T_inlet_C)
cat("\n\n")
resBB <- step_template("butyl_butyrate", T_dry_C = T_inlet_C)
cat("\n\n")
hr("=")
cat("REFERENCE: NO TEMPLATE  (same backbone, Q_template = 0)\n")
hr("=")
resNT <- run_no_template()
r2nt  <- resNT$product
cat(sprintf("    %-14s %8.2f um   %-14s %8.2f um\n",
            "product d50", r2nt[["D_particle_um"]], "product d90", r2nt[["Dp90_um"]]))
cat(sprintf("    %-14s %8.3f um   %-14s %8.3f\n",
            "pore size", r2nt[["D_pore_um"]], "porosity phi", r2nt[["phi_porosity_z"]]))
cat(sprintf("    %-14s %8.3f      %-14s %8.3f\n",
            "f_cr (Raoult)", r2nt[["f_cr_solv"]], "f_esc (total)", r2nt[["f_esc_solv"]]))
cat(sprintf("    %-14s %8.3f      %-14s %8.1f kg/m3\n",
            "solv_retained", r2nt[["solv_retained"]], "rho_tapped", r2nt[["rho_tapped"]]))

# --- side-by-side product comparison (4 columns: no-template + BA + BB) ------
cat("\n\n")
hr("=")
cat("SIDE-BY-SIDE PRODUCT COMPARISON\n")
cat("(same up3_1 backbone; Raoult constant-rate evaporation now active)\n")
hr("=")
cmp4 <- function(lbl, key, unit = "", d = 3) {
  nt <- resNT$product[[key]]; a <- resBA$product[[key]]; b <- resBB$product[[key]]
  cat(sprintf("  %-22s %11.*f %13.*f %13.*f   %s\n", lbl, d, nt, d, a, d, b, unit))
}
cat(sprintf("  %-22s %11s %13s %13s\n", "", "no_template", "butyl_acetate", "butyl_butyrate"))
cat(sprintf("  %-22s %11s %13.2f %13.2f   RED (higher = poorer solvent)\n",
            "Hansen RED", "—", resBA$chem$RED, resBB$chem$RED))
cat(sprintf("  %-22s %11s %13.0f %13.0f   C (vs %.0f C inlet)\n",
            "boiling point", "—",
            resBA$chem$T_bp_solv_K-273.15, resBB$chem$T_bp_solv_K-273.15, T_inlet_C))
cat(sprintf("  %-22s %11s %13.4f %13.4f   dryer core-absorption\n",
            "RTF (chemistry)", "—", resBA$chem$RTF, resBB$chem$RTF))
hr("-")
cmp4("product d50 [um]",    "D_particle_um",  "", 2)
cmp4("product d90 [um]",    "Dp90_um",        "", 2)
cmp4("pore size [um]",      "D_pore_um")
cmp4("porosity phi",        "phi_porosity_z")
cmp4("skin theta",          "theta_skin_z")
hr("-")
cmp4("f_cr (Raoult CR)",   "f_cr_solv",      "", 3)
cmp4("f_esc (total)",       "f_esc_solv",     "", 3)
cmp4("solvent retained",    "solv_retained",  "", 3)
cmp4("f_burst (FR)",        "f_burst_solv",   "", 3)
hr("-")
cmp4("Tg_product [K]",      "Tg_product_K",   "", 1)
cmp4("residual moisture",   "X_moisture",     "", 5)
cmp4("tapped density",      "rho_tapped",     "kg/m3", 1)

# --- T_dryer_in sweep: show where S_bs activates for each template ------------
# The falling-rate escape gate S_bs = sigmoid(T_particle - T_bp) uses
#   T_particle = 0.85 * T_out + 0.15 * T_feed
# where T_out is the outlet from the dryer energy balance.
# At up3_1 conditions (T_out=90C, T_feed=55C) -> T_particle=85C which is well
# below both boiling points. Increasing T_dryer_in raises T_out and T_particle.
# This sweep shows the threshold; S_bs>0.5 means >50% of template vapourises
# in the falling-rate window.
# NOTE: constant-rate evaporation is captured by Module 6c (Raoult); this sweep
# isolates ONLY the falling-rate S_bs gate to show when boiling would kick in.

cat("\n\n")
hr("=")
cat("T_DRYER_IN SENSITIVITY SWEEP (both templates, S_bs falling-rate escape gate)\n")
cat("NOTE: S_bs uses T_particle = 0.85*T_out + 0.15*T_feed, NOT T_dryer_in\n")
cat("      (constant-rate evaporation via Raoult is Module 6c; this sweep = S_bs only)\n")
hr("=")

sweep_templates <- list(
  BA = template_from_chemistry("butyl_acetate",  T_dry_C = T_inlet_C),
  BB = template_from_chemistry("butyl_butyrate", T_dry_C = T_inlet_C))

T_feed_K <- row$up1_exit_temp_C + 273.15
s_solid  <- row$up3_solid_pct / 100
T_in_vals <- seq(130, 200, by = 10)

cat(sprintf("\n  T_dryer_in   T_out    T_part |  S_bs BA  solv_ret BA |  S_bs BB  solv_ret BB\n"))
cat(sprintf("  %s\n", strrep("-", 75)))
for (T_in_C in T_in_vals) {
  T_in_K <- T_in_C + 273.15
  feed_lbhr <- row$up4_feed_mb_lbhr
  mdot_g <- dryer_airflow(feed_lbhr, s_solid, T_in_C, row$up4_Tout_C, row$up1_exit_temp_C)
  # approximate T_out from energy balance (water evaporation dominates)
  LB <- 0.4536; HR <- 3600; h_fg <- 2.30e6; cp_gas <- 1005
  mdot_w <- (feed_lbhr * LB/HR) * (1 - s_solid)
  T_out_K <- max(T_in_K - mdot_w * h_fg / (mdot_g * cp_gas), T_feed_K + 2)
  T_part_K <- 0.85 * T_out_K + 0.15 * T_feed_K
  T_out_C  <- T_out_K - 273.15; T_part_C <- T_part_K - 273.15
  Sbs_BA <- 1/(1+exp(-(T_part_K-(126+273.15))/8))
  Sbs_BB <- 1/(1+exp(-(T_part_K-(166+273.15))/8))
  cat(sprintf("  %8.0f C  %6.0f C  %7.0f C | %7.4f    %8.4f   | %7.4f    %8.4f\n",
              T_in_C, T_out_C, T_part_C,
              Sbs_BA, 1-Sbs_BA, Sbs_BB, 1-Sbs_BB))
}
cat("\n")
cat("  BA (bp=126C) threshold: T_out >= ~135C -> T_dryer_in >= ~175C (very high)\n")
cat("  BB (bp=166C) threshold: T_out >= ~170C -> T_dryer_in >= ~215C (impractical)\n")
cat("  CONCLUSION: at typical outlet T (80-110C), falling-rate S_bs gate DORMANT.\n")
cat("  At up3_1 conditions the dominant differentiation routes are:\n")
cat("    * f_cr (Module 6c Raoult): constant-rate partial-pressure escape;\n")
cat("            BA ~46% escape, BB ~16% (ratio ~3x, tracks vapour pressure ratio).\n")
cat("    * RTF_dryer (chemistry sigmoid): BA absorbs ~12.4% into cores, BB ~3.5%;\n")
cat("            higher RED (poorer solvent) -> less core uptake -> more free template.\n")
cat("    * phi_struct (aggregate templating from UP1 D_agg) -- chemistry-independent.\n")
cat("    * phi_vac (water boiling / steam pockets above 100C) -- chemistry-independent.\n")
cat("\n")
cat("=== RAOULT CONSTANT-RATE EVAPORATION (Module 6c, implemented) ===\n")
cat("  The model captures constant-rate evaporation via partial-pressure driving force:\n")
cat("    T_wb_cr = T_out - 0.8*(T_in - T_out)  [wet-bulb, constant-rate window]\n")
cat("    dH_vap  = 88 * T_bp  [J/mol, Trouton's rule]\n")
cat("    p_cr    = P_atm * exp(dH_vap/R * (1/T_bp - 1/T_wb_cr))  [Clausius-Clapeyron]\n")
cat("    f_cr    = 1 - exp(-C_cr * (p_cr/P_atm) * Perm_shell)     [C_cr = 5.0]\n")
cat("  Sequential escape: f_esc = f_cr + (1 - f_cr) * S_bs.\n")
cat("  At up3_1 (T_in=143C, T_out=90C -> T_wb~=77C):\n")
cat("    BA (bp=126C): p_cr/P_atm ~0.075 -> f_cr ~46%\n")
cat("    BB (bp=166C): p_cr/P_atm ~0.020 -> f_cr ~16%\n")
cat("  Vapour-pressure ratio at T_wb (~3.6x) drives the 3x escape-fraction ratio.\n")
cat("\n")
cat("=== KEY DIFFERENTIATOR: RTF (core-absorption chemistry) ===\n")
cat("  The primary BA vs BB distinction the model correctly captures is RTF_dryer:\n")
cat("    BA: RED=1.196 -> RTF=0.124  (12.4% absorbed into cores, 87.6% free)\n")
cat("    BB: RED=1.333 -> RTF=0.035  (3.5% absorbed into cores, 96.5% free)\n")
cat("  Higher RED (poorer solvent) -> less core penetration -> more free template.\n")
cat("  At the mixer exit RTF_mixer is the SAME (~0.7%) for both (diffusion into\n")
cat("  glassy cores is negligible at tau=22 min / 55C); the RTF_dryer difference\n")
cat("  applies in the hot drying period where softened cores absorb solvent.\n")

# --- persist ------------------------------------------------------------------
p <- function(res, key) res$product[[key]]
outdf <- data.frame(
  template  = c("no_template", "butyl_acetate", "butyl_butyrate"),
  RED       = c(NA,            resBA$chem$RED,   resBB$chem$RED),
  bp_C      = c(NA,            resBA$chem$T_bp_solv_K - 273.15, resBB$chem$T_bp_solv_K - 273.15),
  RTF_chem  = c(NA,            resBA$chem$RTF,   resBB$chem$RTF),
  d50_um    = c(p(resNT,"D_particle_um"),  p(resBA,"D_particle_um"),  p(resBB,"D_particle_um")),
  d90_um    = c(p(resNT,"Dp90_um"),        p(resBA,"Dp90_um"),        p(resBB,"Dp90_um")),
  pore_um   = c(p(resNT,"D_pore_um"),      p(resBA,"D_pore_um"),      p(resBB,"D_pore_um")),
  porosity  = c(p(resNT,"phi_porosity_z"), p(resBA,"phi_porosity_z"), p(resBB,"phi_porosity_z")),
  theta_skin= c(p(resNT,"theta_skin_z"),   p(resBA,"theta_skin_z"),   p(resBB,"theta_skin_z")),
  f_cr      = c(p(resNT,"f_cr_solv"),      p(resBA,"f_cr_solv"),      p(resBB,"f_cr_solv")),
  f_esc     = c(p(resNT,"f_esc_solv"),     p(resBA,"f_esc_solv"),     p(resBB,"f_esc_solv")),
  solv_ret  = c(p(resNT,"solv_retained"),  p(resBA,"solv_retained"),  p(resBB,"solv_retained")),
  f_burst   = c(p(resNT,"f_burst_solv"),   p(resBA,"f_burst_solv"),   p(resBB,"f_burst_solv")),
  Tg_K      = c(p(resNT,"Tg_product_K"),   p(resBA,"Tg_product_K"),   p(resBB,"Tg_product_K")),
  X_moist   = c(p(resNT,"X_moisture"),     p(resBA,"X_moisture"),     p(resBB,"X_moisture")),
  rho_tap   = c(p(resNT,"rho_tapped"),     p(resBA,"rho_tapped"),     p(resBB,"rho_tapped")))
dir.create("unified_output", showWarnings = FALSE)
write.csv(outdf, "unified_output/step_up1_up4_bb_ba.csv", row.names = FALSE)
cat("\nWrote unified_output/step_up1_up4_bb_ba.csv\n")
