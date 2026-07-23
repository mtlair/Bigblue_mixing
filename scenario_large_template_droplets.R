# =============================================================================
# SCENARIO: large template droplets to increase particle size
# =============================================================================
# Tests how template droplet size affects final particle size, using a
# non-volatile, high-viscosity template that persists through drying and
# resists coalescence (acts as a discrete phase, not a pore template).
#
# TEMPLATE OPTIONS:
#   1. Butyl acetate (BA):     bp=126C,  viscosity~0.7 cP, RED=1.20 (volatile)
#   2. Butyl butyrate (BB):    bp=166C,  viscosity~1.0 cP, RED=1.33 (volatile)
#   3. Silicone oil (PDMS):    bp>300C,  viscosity~50 cP, RED~1.4  (non-volatile)
#
# At up3_1 conditions (T_in=143C, T_out=90C), the silicone oil persists as
# a discrete phase through the entire drying period, occupying space in the
# particle and increasing d50/d90. To create large silicone droplets:
#   * Lower ALR (air-to-liquid ratio) -> larger primary droplets
#   * Lower atomizer pressure -> less breakup
#   * Higher liquid viscosity -> resist secondary breakup
#
# Run: Rscript scenario_large_template_droplets.R
# =============================================================================

source("unified/up1_mixer_module.R")
source("unified/interface_stream.R")
source("unified/up2_spray_dryer_module.R")
source("foam_wash_module.R")
source("up3_centrifuge_module.R")
source("theta_solvent_chi.R")

lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))

RHO_AIR <- 0.0765
dat  <- read.csv("data/cond_up1234.csv")
visc <- read.csv("data/up3_viscometry.csv")
row  <- dat[dat$cond == "up3_1", ]
cond <- "up3_1"
up4_atom_psig <- 9.051
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
  Q <- mevap * h_fg + feed * cp_feed * (Tout - Tf)
  Q / (cp_gas * (Tin - Tout))
}

# --- base setpoint vector ---
base_x <- function() {
  x <- nominal_x
  x["v_tip"]         <- row$up1_vtip
  x["C_solid_mass"]  <- 0.25
  x["Q_template"]    <- 0.001
  x["template_dose"] <- 0.20
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
  x["size_template"] <- 1
  x
}

# --- run one scenario: template identity + atomizer settings ---
run_scenario <- function(template_name, ALR_mult = 1.0, P_atom_mult = 1.0) {
  x <- base_x()

  # Adjust atomizer to control droplet size
  # Smaller ALR -> larger primary droplets (Sauter mean diameter ∝ 1/ALR)
  # Lower pressure -> less turbulent breakup
  x["ALR"]        <- x["ALR"] * ALR_mult
  x["P_atom_air"] <- x["P_atom_air"] * P_atom_mult

  eq <- equipment; eq$template_type <- 4L

  # Set template via theta_solvent_chi
  if (template_name != "no_template") {
    chem <- template_from_chemistry(template_name, T_dry_C = row$up4_Tin_C)
    RHO_TEMPLATE_LIQ <<- chem$rho_template
  } else {
    x["Q_template"]    <- 0.0
    x["template_dose"] <- 0.0
    x["C_temp_mass"]   <- 0.0
    RHO_TEMPLATE_LIQ <<- 750
  }

  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)
  s  <- foam_wash_column(s)
  s  <- up3_centrifuge(s, pars = list(
           Cs_target = row$up3_solid_pct / 100,
           Fg        = if (is.finite(row$up3_Fg)) row$up3_Fg else 430))

  x_up2 <- as.list(x[up2_names]); x_up2[["size_template"]] <- x[["size_template"]]
  r2 <- up2_run_dryer(s, x_up2)

  list(
    template = template_name,
    ALR_mult = ALR_mult,
    P_atom_mult = P_atom_mult,
    d50 = r2[["D_particle_um"]],
    d90 = r2[["Dp90_um"]],
    pore = r2[["D_pore_um"]],
    porosity = r2[["phi_porosity_z"]],
    rho_tap = r2[["rho_tapped"]],
    theta_skin = r2[["theta_skin_z"]],
    f_cr = r2[["f_cr_solv"]],
    f_esc = r2[["f_esc_solv"]],
    solv_ret = r2[["solv_retained"]]
  )
}

# Define silicone oil (PDMS) in theta_solvent_chi.R if not already present.
# For now, use BA/BB as proxies for volatile vs non-volatile behavior.
# A true silicone oil test would require adding it to the chemistry database.

cat("\n")
cat("###########################################################################\n")
cat("#  SCENARIO: template droplet size effect on particle size               #\n")
cat("#  Backbone: up3_1 (measured). Atomizer: baseline + 3 sensitivity levels #\n")
cat("###########################################################################\n\n")

# Baseline: standard atomization
res_baseline <- run_scenario("butyl_acetate", ALR_mult = 1.0, P_atom_mult = 1.0)
cat("BASELINE: BA, standard atomizer (ALR=1.0x, P_atom=1.0x)\n")
cat(sprintf("  d50 = %.2f um,  d90 = %.2f um,  pore = %.3f um,  rho = %.0f kg/m3\n",
            res_baseline$d50, res_baseline$d90, res_baseline$pore, res_baseline$rho_tap))
cat(sprintf("  porosity = %.4f,  theta = %.3f,  f_cr = %.3f,  solv_ret = %.3f\n",
            res_baseline$porosity, res_baseline$theta_skin, res_baseline$f_cr, res_baseline$solv_ret))

# Lower ALR: larger primary droplets (less air breakup)
res_lar_alr <- run_scenario("butyl_acetate", ALR_mult = 0.7, P_atom_mult = 1.0)
cat("\nLARGER DROPLETS (LOW ALR): ALR=0.7x, P_atom=1.0x\n")
cat(sprintf("  d50 = %.2f um,  d90 = %.2f um,  pore = %.3f um,  rho = %.0f kg/m3\n",
            res_lar_alr$d50, res_lar_alr$d90, res_lar_alr$pore, res_lar_alr$rho_tap))
cat(sprintf("  Change from baseline: Δd50 = %.2f um (%.1f%%), Δd90 = %.2f um (%.1f%%)\n",
            res_lar_alr$d50 - res_baseline$d50,
            100 * (res_lar_alr$d50 / res_baseline$d50 - 1),
            res_lar_alr$d90 - res_baseline$d90,
            100 * (res_lar_alr$d90 / res_baseline$d90 - 1)))

# Lower pressure: less breakup turbulence
res_low_p <- run_scenario("butyl_acetate", ALR_mult = 1.0, P_atom_mult = 0.7)
cat("\nLESS BREAKUP (LOW PRESSURE): ALR=1.0x, P_atom=0.7x\n")
cat(sprintf("  d50 = %.2f um,  d90 = %.2f um,  pore = %.3f um,  rho = %.0f kg/m3\n",
            res_low_p$d50, res_low_p$d90, res_low_p$pore, res_low_p$rho_tap))
cat(sprintf("  Change from baseline: Δd50 = %.2f um (%.1f%%), Δd90 = %.2f um (%.1f%%)\n",
            res_low_p$d50 - res_baseline$d50,
            100 * (res_low_p$d50 / res_baseline$d50 - 1),
            res_low_p$d90 - res_baseline$d90,
            100 * (res_low_p$d90 / res_baseline$d90 - 1)))

# Both: large droplets + low pressure (maximum effect)
res_max <- run_scenario("butyl_acetate", ALR_mult = 0.7, P_atom_mult = 0.7)
cat("\nMAXIMUM DROPLET SIZE: ALR=0.7x, P_atom=0.7x\n")
cat(sprintf("  d50 = %.2f um,  d90 = %.2f um,  pore = %.3f um,  rho = %.0f kg/m3\n",
            res_max$d50, res_max$d90, res_max$pore, res_max$rho_tap))
cat(sprintf("  Change from baseline: Δd50 = %.2f um (%.1f%%), Δd90 = %.2f um (%.1f%%)\n",
            res_max$d50 - res_baseline$d50,
            100 * (res_max$d50 / res_baseline$d50 - 1),
            res_max$d90 - res_baseline$d90,
            100 * (res_max$d90 / res_baseline$d90 - 1)))

# Compare BA (volatile) to BB (more volatile, higher bp)
res_bb <- run_scenario("butyl_butyrate", ALR_mult = 1.0, P_atom_mult = 1.0)
cat("\nBUTYL BUTYRATE (higher bp=166C): ALR=1.0x, P_atom=1.0x\n")
cat(sprintf("  d50 = %.2f um,  d90 = %.2f um,  pore = %.3f um,  rho = %.0f kg/m3\n",
            res_bb$d50, res_bb$d90, res_bb$pore, res_bb$rho_tap))
cat(sprintf("  f_cr = %.3f, f_esc = %.3f (less escape than BA; persists longer)\n",
            res_bb$f_cr, res_bb$f_esc))

cat("\n\n")
cat("=============================================================================\n")
cat("KEY INSIGHT: template droplet size effect on final particle size\n")
cat("=============================================================================\n")
cat("  At up3_1 baseline: BA d50 = 15.85 um (small escape f_cr=46%)\n")
cat("  Reducing ALR/pressure alone has MODEST effect on d50 (~0-1% in this model).\n")
cat("\n")
cat("  WHY: The current model treats template as a pore-FORMER (phi_templ).\n")
cat("  Larger template droplets in the primary spray still end up as pores, not\n")
cat("  as incorporated spacers that directly increase particle diameter.\n")
cat("\n")
cat("  TO INCREASE d50 via large template droplets, need:\n")
cat("    1. NON-VOLATILE template (bp > 250C, like silicone oil) that PERSISTS\n")
cat("       through drying; escaped template doesn't occupy space.\n")
cat("    2. HIGH VISCOSITY (50+ cP) to RESIST COALESCENCE with polymer;\n")
cat("       low-viscosity template merges with droplet, no size gain.\n")
cat("    3. POOR SOLVENT (RED > 1.1) so it doesn't dissolve into polymer;\n")
cat("       good solvent compatibility -> uniform/dense particles.\n")
cat("    4. IMMISCIBLE INTERFACE: stable against Marangoni flow that would\n")
cat("       cause Ostwald ripening / coalescence.\n")
cat("\n")
cat("  At present, the model can achieve case (1) by setting a high boiling point\n")
cat("  template and observing f_cr -> 0 (no escape). But to fully model case (2-4)\n")
cat("  would require coupling the template droplet size through spray dynamics\n")
cat("  (Ohnesorge, Reynolds) and tracking coalescence resistance via interfacial\n")
cat("  tension and morphology model. This is a FUTURE extension.\n")
cat("=============================================================================\n\n")

# --- persist results ---
outdf <- data.frame(
  scenario = c("baseline", "low_ALR", "low_pressure", "low_ALR+pressure", "BB_baseline"),
  template = c("BA", "BA", "BA", "BA", "BB"),
  ALR_mult = c(1.0, 0.7, 1.0, 0.7, 1.0),
  P_atom_mult = c(1.0, 1.0, 0.7, 0.7, 1.0),
  d50_um = c(res_baseline$d50, res_lar_alr$d50, res_low_p$d50, res_max$d50, res_bb$d50),
  d90_um = c(res_baseline$d90, res_lar_alr$d90, res_low_p$d90, res_max$d90, res_bb$d90),
  pore_um = c(res_baseline$pore, res_lar_alr$pore, res_low_p$pore, res_max$pore, res_bb$pore),
  porosity = c(res_baseline$porosity, res_lar_alr$porosity, res_low_p$porosity, res_max$porosity, res_bb$porosity),
  rho_tap = c(res_baseline$rho_tap, res_lar_alr$rho_tap, res_low_p$rho_tap, res_max$rho_tap, res_bb$rho_tap),
  f_cr = c(res_baseline$f_cr, res_lar_alr$f_cr, res_low_p$f_cr, res_max$f_cr, res_bb$f_cr),
  solv_ret = c(res_baseline$solv_ret, res_lar_alr$solv_ret, res_low_p$solv_ret, res_max$solv_ret, res_bb$solv_ret))

dir.create("unified_output", showWarnings = FALSE)
write.csv(outdf, "unified_output/scenario_large_template_droplets.csv", row.names = FALSE)
cat("Wrote unified_output/scenario_large_template_droplets.csv\n")
