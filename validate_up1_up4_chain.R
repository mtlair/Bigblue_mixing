# =============================================================================
# UP1 -> UP4 FULL-CHAIN validation harness (4 measured chain conditions)
# =============================================================================
# The 4 conditions in data/cond_up1234.csv were measured through the COMPLETE
# train UP1 (mixer) -> UP2 (foam wash) -> UP3 (centrifuge) -> UP4 (spray dryer),
# with SEM samples 114641-114644 = up3_1..up3_4 and a measured product PSD.
#
# This harness replays UP1 from its measured setpoints, then hands the dryer the
# CORRECT post-UP3 feed state (concentrated 25% -> ~40% solid, gas-free density,
# measured post-UP3 viscosity, residual cake gas) rather than the raw UP1 slurry.
# It is the "validate UP1 -> UP4" step: UP2/UP3 are represented by their MEASURED
# outputs (injection), not yet by their models -- that wiring is the next step
# (intermediate_stage_1/2 in unified/interface_stream.R).
#
# Compared against measurement:
#   dry d50  : model up2 D_particle_um  vs  measured d50   (product PSD)
#   dry d90  : model up2 Dp90_um        vs  measured d90
#
# Per-condition inputs from the data:
#   UP1 : up1_vtip, up1_exit_temp, up1_psig                (cond_up1234.csv)
#   UP3 : up3_solid_pct, cake gas holdup, post-UP3 viscosity(cond_up1234 + up3_viscometry.csv)
#   UP4 : up4_feed (mass-balance), up4_atom_scfm, atom psig, pump psig, Tin/Tout
#
# Run:  Rscript validate_up1_up4_chain.R
# =============================================================================

# --- load unified-model definitions (functions only; skip the Morris driver) -
lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))

# --- UP2 / UP3 stage models (for the predictive "wired" path) ----------------
source("foam_wash_module.R")      # UP2 foam-wash column: foam_wash_column(stream, pars)
source("up3_centrifuge_module.R") # UP3 centrifuge:       up3_centrifuge(stream, pars)

dat  <- read.csv("data/cond_up1234.csv")
visc <- read.csv("data/up3_viscometry.csv")
RHO_AIR <- 0.0765   # lb/ft3, standard air
RHO_S   <- 1700     # bare polymer colloid / dried-product skeletal density

# atomizer + feed-line pressures (from up1_2_3_4_visc_sem.xlsx `up1_2_3_4_psd`,
# not carried in cond_up1234.csv). psig gauge.
up4_atom_psig <- c(up3_1 = 9.051, up3_2 = 5.716, up3_3 = 7.289, up3_4 = 10.957)
up4_pump_psig <- c(up3_1 = 5.515, up3_2 = 6.817, up3_3 = 6.821, up3_4 = 6.805)

# post-UP3 slurry viscosity at the ~12.7 1/s reference used for UP1 mu_exit_PaS
# (log-log interpolation of the measured flow curve).
visc_ref <- function(cond, target = 12.7) {
  s <- visc[visc$cond == cond, ]
  s <- s[order(s$shear_rate_1s), ]
  exp(approx(log(s$shear_rate_1s), log(s$viscosity_PaS), log(target), rule = 2)$y)
}

# Dryer air mass flow from the evaporative energy balance (same closure as the
# direct harness; feed water -> <=0.5% residual moisture at measured Tin/Tout).
dryer_airflow <- function(feed_lbhr, s, Tin, Tout, Tf) {
  h_fg <- 2.30e6; cp_gas <- 1005; resid <- 0.005; LB <- 0.4536; HR <- 3600
  feed <- feed_lbhr * LB / HR
  mevap <- feed * (1 - s) - resid / (1 - resid) * feed * s
  cp_feed <- (1 - s) * 4180 + s * 1500
  Q <- mevap * h_fg + feed * cp_feed * (Tout - Tf)
  Q / (cp_gas * (Tin - Tout))
}

# --- one condition, UP1 -> [measured post-UP3 feed state] -> UP4 -------------
run_chain <- function(row, size_template = 0) {
  cond <- row$cond
  x <- nominal_x
  # UP1 mixer setpoints
  x["v_tip"]         <- row$up1_vtip
  x["C_solid_mass"]  <- 0.25
  x["Q_template"]    <- 0.001
  x["template_dose"] <- 0.02
  x["T_system"]      <- row$up1_exit_temp_C + 273.15
  x["P_system"]      <- (row$up1_psig + 14.696) / 14.696
  # UP4 atomizer / dryer setpoints
  air_lbhr  <- row$up4_atom_scfm * RHO_AIR * 60
  feed_lbhr <- row$up4_feed_mb_lbhr
  x["ALR"]        <- air_lbhr / feed_lbhr
  x["mdot_L"]     <- feed_lbhr * 0.4536 / 3600
  x["P_atom_air"] <- (up4_atom_psig[[cond]] + 14.696) * 6894.76
  x["P_feed"]     <- (up4_pump_psig[[cond]] + 14.696) * 6894.76
  x["T_dryer_in"] <- row$up4_Tin_C + 273.15
  s_solid <- row$up3_solid_pct / 100
  x["mdot_gas_dry"] <- dryer_airflow(feed_lbhr, s_solid, row$up4_Tin_C,
                                     row$up4_Tout_C, row$up1_exit_temp_C)
  x["size_template"] <- size_template

  eq <- equipment; eq$template_type <- TEMPLATE_TYPE
  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)

  # --- inject the MEASURED post-UP3 feed state (UP2/UP3 by measurement) ------
  alpha_cake <- if (!is.null(row$up3_alphag_cake_pct) &&
                    is.finite(row$up3_alphag_cake_pct)) row$up3_alphag_cake_pct / 100 else 0.076
  s$C_solid     <- s_solid
  s$rho_slurry  <- 1 / (s_solid / RHO_S + (1 - s_solid) / 1000)  # gas-free condensed
  s$alpha_g     <- alpha_cake
  s$mu_exit_PaS <- visc_ref(cond)                                # post-UP3 low-shear eta

  x_up2 <- as.list(x[up2_names])
  x_up2[["size_template"]] <- size_template
  r2 <- up2_run_dryer(s, x_up2)

  list(d50 = unname(r2[["D_particle_um"]]),
       d90 = unname(r2[["Dp90_um"]]),
       theta_skin = unname(r2[["theta_skin_z"]]),
       phi_por    = unname(r2[["phi_porosity_z"]]),
       omega      = unname(r2[["Omega_struct_z"]]),
       Xmoist     = unname(r2[["X_moisture"]]),
       mu_ref     = s$mu_exit_PaS,
       Csol       = s$C_solid)
}

# --- predictive WIRED path: UP1 -> UP2(foam wash) -> UP3(centrifuge) -> UP4 ---
# Same UP1/UP4 setpoints, but UP2/UP3 are the MODELS (not measured injection).
run_wired <- function(row, size_template = 1) {
  cond <- row$cond
  x <- nominal_x
  x["v_tip"]         <- row$up1_vtip
  x["C_solid_mass"]  <- 0.25
  x["Q_template"]    <- 0.001
  x["template_dose"] <- 0.02
  x["T_system"]      <- row$up1_exit_temp_C + 273.15
  x["P_system"]      <- (row$up1_psig + 14.696) / 14.696
  air_lbhr  <- row$up4_atom_scfm * RHO_AIR * 60
  feed_lbhr <- row$up4_feed_mb_lbhr
  x["ALR"]        <- air_lbhr / feed_lbhr
  x["mdot_L"]     <- feed_lbhr * 0.4536 / 3600
  x["P_atom_air"] <- (up4_atom_psig[[cond]] + 14.696) * 6894.76
  x["P_feed"]     <- (up4_pump_psig[[cond]] + 14.696) * 6894.76
  x["T_dryer_in"] <- row$up4_Tin_C + 273.15
  s_solid <- row$up3_solid_pct / 100
  x["mdot_gas_dry"] <- dryer_airflow(feed_lbhr, s_solid, row$up4_Tin_C,
                                     row$up4_Tout_C, row$up1_exit_temp_C)
  x["size_template"] <- size_template

  eq <- equipment; eq$template_type <- TEMPLATE_TYPE
  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)

  # UP2 foam wash (gas removal, pressure); UP3 centrifuge (concentrate to cake)
  s <- foam_wash_column(s)
  s <- up3_centrifuge(s, pars = list(
         Cs_target = row$up3_solid_pct / 100,       # per-condition measured cake solids
         Fg        = if (is.finite(row$up3_Fg)) row$up3_Fg else 430))

  x_up2 <- as.list(x[up2_names]); x_up2[["size_template"]] <- size_template
  r2 <- up2_run_dryer(s, x_up2)
  list(d50 = unname(r2[["D_particle_um"]]), d90 = unname(r2[["Dp90_um"]]),
       Csol = s$C_solid, mu = s$mu_exit_PaS, alpha = s$alpha_g, rho = s$rho_slurry)
}

res0 <- lapply(seq_len(nrow(dat)), function(i) run_chain(dat[i, ], size_template = 0))
res1 <- lapply(seq_len(nrow(dat)), function(i) run_chain(dat[i, ], size_template = 1))
resW <- lapply(seq_len(nrow(dat)), function(i) run_wired(dat[i, ], size_template = 1))

cat("\n============ UP1 -> UP4 FULL CHAIN: model vs measured (4 conditions) ============\n")
cat("(UP2/UP3 represented by measured post-UP3 feed state: solid%, gas-free rho, eta)\n\n")

cat("--- feed state handed to the dryer ---\n")
cat(sprintf("%6s %10s %12s %14s\n", "cond", "up3_solid%", "eta_ref[Pa.s]", "rho_nogas[kg/m3]"))
for (i in seq_len(nrow(dat))) {
  r <- dat[i, ]; s_solid <- r$up3_solid_pct / 100
  cat(sprintf("%6s %10.1f %12.2f %14.0f\n", r$cond, r$up3_solid_pct,
              res0[[i]]$mu_ref, 1 / (s_solid / RHO_S + (1 - s_solid) / 1000)))
}

cat("\n--- product d50 [um] ---\n")
cat(sprintf("%6s %9s %12s %14s\n", "cond", "meas", "model(mute)", "model(tmpl=1)"))
for (i in seq_len(nrow(dat)))
  cat(sprintf("%6s %9.2f %12.2f %14.2f\n", dat$cond[i], dat$d50[i],
              res0[[i]]$d50, res1[[i]]$d50))

cat("\n--- product d90 [um] ---\n")
cat(sprintf("%6s %9s %12s %14s\n", "cond", "meas", "model(mute)", "model(tmpl=1)"))
for (i in seq_len(nrow(dat)))
  cat(sprintf("%6s %9.2f %12.2f %14.2f\n", dat$cond[i], dat$d90[i],
              res0[[i]]$d90, res1[[i]]$d90))

cat("\n--- morphology (model; for SEM cross-check once wired) ---\n")
cat(sprintf("%6s %11s %11s %11s %10s\n", "cond", "theta_skin", "phi_poros", "Omega_sph", "X_moist"))
for (i in seq_len(nrow(dat)))
  cat(sprintf("%6s %11.3f %11.4f %11.3f %10.5f\n", dat$cond[i],
      res0[[i]]$theta_skin, res0[[i]]$phi_por, res0[[i]]$omega, res0[[i]]$Xmoist))

# --- error metrics ------------------------------------------------------------
logacc <- function(meas, pred) {
  ok <- is.finite(meas) & is.finite(pred) & meas > 0 & pred > 0
  exp(sqrt(mean((log(pred[ok]) - log(meas[ok]))^2)))
}
d50_meas <- dat$d50; d90_meas <- dat$d90
d50_p0 <- sapply(res0, `[[`, "d50"); d50_p1 <- sapply(res1, `[[`, "d50")
d90_p0 <- sapply(res0, `[[`, "d90"); d90_p1 <- sapply(res1, `[[`, "d90")

cat("\n--- geometric RMS ratio (1.0 = perfect) ---\n")
cat(sprintf("  d50, muted   : %.2f\n", logacc(d50_meas, d50_p0)))
cat(sprintf("  d50, tmpl=1  : %.2f\n", logacc(d50_meas, d50_p1)))
cat(sprintf("  d90, muted   : %.2f\n", logacc(d90_meas, d90_p0)))
cat(sprintf("  d90, tmpl=1  : %.2f\n", logacc(d90_meas, d90_p1)))

cat("\n--- PREDICTIVE WIRED chain: UP1 -> UP2(foam wash) -> UP3(centrifuge) -> UP4 ---\n")
cat("(UP2/UP3 are MODELS here; compare product d50 to measured and to injection)\n")
cat(sprintf("%6s %8s %10s %9s %9s   %8s %10s %10s\n",
            "cond", "Cs_out", "eta[Pa.s]", "alpha_g", "meas d50", "wired", "inject", "meas d90"))
for (i in seq_len(nrow(dat))) {
  w <- resW[[i]]
  cat(sprintf("%6s %8.3f %10.2f %9.3f %9.2f   %8.2f %10.2f %10.2f\n",
      dat$cond[i], w$Csol, w$mu, w$alpha, dat$d50[i], w$d50, res1[[i]]$d50, w$d90))
}
d50_w <- sapply(resW, `[[`, "d50"); d90_w <- sapply(resW, `[[`, "d90")
cat(sprintf("\n  wired d50 geometric RMS ratio : %.2f\n", logacc(d50_meas, d50_w)))
cat(sprintf("  wired d90 geometric RMS ratio : %.2f\n", logacc(d90_meas, d90_w)))

out <- data.frame(cond = dat$cond, up3_solid_pct = dat$up3_solid_pct,
  eta_ref = round(sapply(res0, `[[`, "mu_ref"), 2),
  d50_meas = d50_meas, d50_muted = round(d50_p0, 2), d50_tmpl = round(d50_p1, 2),
  d90_meas = d90_meas, d90_muted = round(d90_p0, 2), d90_tmpl = round(d90_p1, 2),
  theta_skin = round(sapply(res0, `[[`, "theta_skin"), 3),
  phi_porosity = round(sapply(res0, `[[`, "phi_por"), 4),
  Omega = round(sapply(res0, `[[`, "omega"), 3))
dir.create("unified_output", showWarnings = FALSE)
write.csv(out, "unified_output/up1_up4_chain_validation.csv", row.names = FALSE)
cat("\nWrote unified_output/up1_up4_chain_validation.csv\n")
