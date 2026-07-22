# =============================================================================
# UP1 -> UP4 DIRECT validation harness
# =============================================================================
# The visc.xlsx calibration runs did NOT pass through UP2 (foam-wash) or UP3
# (centrifuge) -- see DATA_REVIEW.md. This harness replays each of the nine
# measured conditions UP1 (mixer) -> UP4 (dryer) *directly*, bypassing UP2/UP3,
# and compares the calibrated model against the measurements:
#
#   wet d50   (post-UP1)  : model stream$D_agg_um      vs  measured d50 (sample m)
#   low-shear viscosity   : model stream$mu_exit_PaS   vs  measured eta @ ~12.7 1/s
#   dry d50   (post-UP4)  : model up2 D_particle_um    vs  measured d50 (sample p)
#
# Per-condition inputs taken from the data (data/cond_process.csv):
#   UP1 : v_tip, feed solids, exit temperature, exit gauge pressure
#   UP4 : ALR (up4atom_scfm air / up4_feed), atomizing-air pressure (up4_psig),
#         liquid feed mass flow (up4_feed).
# Everything else (formulation chemistry: Delta_pH, ionic strength, surfactant,
# ...) is held at the model nominal -- the data does not resolve it -- so the
# aggregation onset v_tip_crit is the SAME for every condition here.
#
# Run:  Rscript validate_up1_up4_direct.R
# =============================================================================

# --- load unified-model definitions (functions only; skip the Morris driver) -
lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))

dat     <- read.csv("data/cond_process.csv")
RHO_AIR <- 0.0765   # lb/ft3, standard air

# --- one condition, UP1 -> (bypass) -> UP4 -----------------------------------
run_cond <- function(row, size_template = 0) {
  x <- nominal_x
  # UP1 mixer setpoints from the data
  x["v_tip"]         <- row$up1_vtip
  x["C_solid_mass"]  <- row$solid_pct / 100
  x["Q_template"]    <- 0.001                                  # ~no template dilution
  x["template_dose"] <- 0.02
  x["T_system"]      <- row$up1_exit_temp_C + 273.15
  x["P_system"]      <- (row$up1_psig + 14.696) / 14.696       # atm, absolute
  # UP4 atomizer setpoints from the measured flows
  air_lbhr        <- row$up4atom_scfm * RHO_AIR * 60
  x["ALR"]        <- air_lbhr / row$up4_feed_lbhr
  x["P_atom_air"] <- (row$up4_psig + 14.696) * 6894.76         # Pa, absolute
  x["mdot_L"]     <- row$up4_feed_lbhr * 0.4536 / 3600         # kg/s
  x["size_template"] <- size_template

  eq <- equipment; eq$template_type <- TEMPLATE_TYPE
  p1 <- up1_pars_from_x(x)
  r1 <- up1_run_mixer(p1, eq)
  s  <- stream_from_up1(r1, p1, eq)

  # --- BYPASS UP2 + UP3: hand the UP1 stream straight to the dryer ---
  x_up2 <- as.list(x[up2_names])
  x_up2[["size_template"]] <- size_template
  r2 <- up2_run_dryer(s, x_up2)

  list(wet  = s$D_agg_um,
       visc = s$mu_exit_PaS,
       dry  = unname(r2[["D_particle_um"]]),
       crit = s$v_tip_crit,
       Csol = s$C_solid)
}

# --- run all conditions -------------------------------------------------------
res0 <- lapply(seq_len(nrow(dat)), function(i) run_cond(dat[i, ], size_template = 0))
res1 <- lapply(seq_len(nrow(dat)), function(i) run_cond(dat[i, ], size_template = 1))

fmt <- function(v, d = 2) formatC(v, format = "f", digits = d)
cat("\n================ UP1 -> UP4 DIRECT: model vs measured ================\n")
cat(sprintf("v_tip_crit (nominal chemistry) = %.2f m/s  [aggregation onset]\n\n", res0[[1]]$crit))

cat("--- UP1 wet d50 [um] : aggregation closure ---\n")
cat(sprintf("%5s %7s %9s %9s   %s\n", "cond", "v_tip", "meas", "model", "regime"))
for (i in seq_len(nrow(dat))) {
  d <- dat[i, ]; m <- res0[[i]]
  reg <- if (d$up1_vtip < m$crit) "atomizer-control" else "UP1-control"
  cat(sprintf("%5d %7.2f %9.2f %9.2f   %s\n", d$cond, d$up1_vtip, d$wet_d50_um, m$wet, reg))
}

cat("\n--- UP1 low-shear viscosity [Pa.s] ---\n")
cat(sprintf("%5s %10s %10s\n", "cond", "meas@12.7", "model"))
for (i in seq_len(nrow(dat))) cat(sprintf("%5d %10.4f %10.4f\n", dat$cond[i], dat$visc12_PaS[i], res0[[i]]$visc))

cat("\n--- UP4 dry d50 [um] : droplet-shell (muted) vs size-template ---\n")
cat(sprintf("%5s %9s %12s %14s\n", "cond", "meas", "model(mute)", "model(tmpl=1)"))
for (i in seq_len(nrow(dat))) cat(sprintf("%5d %9.2f %12.2f %14.2f\n",
    dat$cond[i], dat$dry_d50_um[i], res0[[i]]$dry, res1[[i]]$dry))

# --- error metrics ------------------------------------------------------------
logacc <- function(meas, pred) {
  ok <- is.finite(meas) & is.finite(pred) & meas > 0 & pred > 0
  exp(sqrt(mean((log(pred[ok]) - log(meas[ok]))^2)))   # geometric RMS ratio
}
wet_meas <- dat$wet_d50_um;  wet_pred <- sapply(res0, `[[`, "wet")
vis_meas <- dat$visc12_PaS;  vis_pred <- sapply(res0, `[[`, "visc")
dry_meas <- dat$dry_d50_um;  dry_p0 <- sapply(res0, `[[`, "dry"); dry_p1 <- sapply(res1, `[[`, "dry")
agg <- dat$cond != 2   # UP1-control subset

cat("\n--- geometric RMS ratio (1.0 = perfect; 1.5 = ~50% typical error) ---\n")
cat(sprintf("  wet d50 (UP1-control 8)    : %.2f\n", logacc(wet_meas[agg], wet_pred[agg])))
cat(sprintf("  wet d50 (all 9, incl cond2): %.2f  [cond2 dispersed reports ODE primary, not the 0.2um bead]\n", logacc(wet_meas, wet_pred)))
cat(sprintf("  viscosity (UP1-control 8)  : %.2f\n", logacc(vis_meas[agg], vis_pred[agg])))
cat(sprintf("  dry d50, muted (UP1-ctrl 8): %.2f\n", logacc(dry_meas[agg], dry_p0[agg])))
cat(sprintf("  dry d50, tmpl=1 (UP1-ctrl 8): %.2f\n", logacc(dry_meas[agg], dry_p1[agg])))
cat(sprintf("  dry d50, muted (cond2 atomizer): meas %.1f vs model %.1f um\n",
            dry_meas[!agg], dry_p0[!agg]))

# --- write a tidy CSV ---------------------------------------------------------
out <- data.frame(cond = dat$cond, v_tip = dat$up1_vtip, solid_pct = dat$solid_pct,
  wet_meas = wet_meas, wet_model = round(wet_pred, 3),
  visc_meas = vis_meas, visc_model = round(vis_pred, 4),
  dry_meas = dry_meas, dry_model_muted = round(dry_p0, 2), dry_model_tmpl = round(dry_p1, 2))
dir.create("unified_output", showWarnings = FALSE)
write.csv(out, "unified_output/up1_up4_direct_validation.csv", row.names = FALSE)
cat("\nWrote unified_output/up1_up4_direct_validation.csv\n")
