#!/usr/bin/env Rscript
# =============================================================================
# FULL-TRAIN MORRIS SCREEN
#   UP1 mixer -> UP2 foam-wash (491-line ODE) -> UP3 separator -> reslurry
#   -> UP4 spray dryer
# =============================================================================
# One Morris elementary-effects screen over the WHOLE train. The per-stage
# factor tables are MERGED into a single dictionary, with three rules:
#
#   1. UP1 output bounds the downstream process. Every downstream input the
#      stream physically determines (solids, gas holdup, densities, surface
#      chemistry, temperature, pressure, bubble size, floc strength, ...) is
#      OVERWRITTEN by the mixer-exit stream inside the stage adapters, so those
#      are NOT independent Morris factors and are dropped from the dictionary.
#   2. Equipment / geometry parameters are FIXED (held at nominal): UP2 column
#      geometry (H_total, H_pool, U_up), UP3 bowl geometry (L_cyl, r_bowl,
#      beach_angle_deg, length_feed). Hard to change -> not screened.
#   3. UP2 is included through its COALESCENCE knobs only (K_coal, K_break,
#      K_burst, d_b_burst, frac_gas_coarse_ref). Coalescence in the wash column
#      resets the bubble population and the gas state, which flows straight into
#      the UP3 / UP4 feed, so it belongs in the screen even with geometry fixed.
#
# Factors are prefixed by stage (up1./up2./up3./up4.) so names cannot collide
# (e.g. D_template exists in both UP1 and UP4).
#
# Every evaluation is NA-GUARDED (tryCatch + finiteness check): a corner that
# fails either ODE returns an all-NA output row instead of aborting the screen.
#
# Screened FEATURES OF INTEREST (after UP4):
#   porosity (phi_porosity_z), sphericity (Omega_struct_z),
#   size (D_particle_um), skin (theta_skin_z), tapped density (rho_tapped).
#
# Run:  Rscript full_train_morris.R           # default r_traj
#       R_TRAJ=20 Rscript full_train_morris.R # override trajectory count
# Deps: deSolve (mixer + column ODEs). Everything else is base R.
# =============================================================================

# ---- load the full-train DEFINITIONS without running its nominal demo -------
.ftL   <- readLines("full_train_mixer_to_dryer.R")
.ftCut <- grep("^# NOMINAL END-TO-END RUN", .ftL)[1]
eval(parse(text = paste(.ftL[seq_len(.ftCut - 1)], collapse = "\n")), envir = globalenv())

set.seed(42)
r_traj  <- as.integer(Sys.getenv("R_TRAJ", "12"))
levels  <- 8
out_dir <- "full_train_output"
dir.create(out_dir, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. MERGED FACTOR DICTIONARY
# -----------------------------------------------------------------------------
# UP3 operating factors kept free (drop the 13 the stream overwrites, fix the
# 4 equipment factors).
up3_free <- c("rpm", "delta_rpm", "r_pool", "flow_rate_lpm", "pop_frac",
              "gas_sat_frac", "t_recover_sec", "contact_angle_deg",
              "chi_parameter", "S_base", "S_ceiling")
# UP4 operating factors kept free. Dropped: the 12 the reslurry handoff
# overwrites, PLUS T_feed / phi_emulsion / D_template, now bound from the mixer
# stream (feed T = mixer-exit T; template emulsion fraction & droplet size come
# from UP1) — so UP1 output bounds them, they are not independent dryer knobs.
up4_free <- c("ALR", "P_system", "P_feed", "mdot_L", "T_system", "mdot_gas_dry",
              "Y_in", "t_hold", "Tg_polymer", "n_flow",
              "k_perm_mono", "k_perm_plast", "k_perm_bind", "T_bp_solv")

# UP2 coalescence factors (equipment/geometry fixed). Ranges bracket the
# calibrated column defaults; broken into decades where they span them (log).
up2_tbl <- data.frame(
  base = c("K_coal", "K_break", "K_burst", "d_b_burst", "frac_gas_coarse_ref"),
  min  = c( 3.0e-5,   3.0e-5,    3.0e-4,    1.5e-3,      0.15),
  max  = c( 6.0e-4,   6.0e-4,    6.0e-3,    6.0e-3,      0.65),
  log  = c( TRUE,     TRUE,      TRUE,      FALSE,       FALSE),
  stringsAsFactors = FALSE)

mk <- function(base, min, max, log, stage)
  data.frame(base = base, min = min, max = max, log = as.logical(log),
             stage = stage, stringsAsFactors = FALSE)

d_up1 <- with(mx_factors[mx_factors$module == "up1", ],
              mk(name, min, max, log, "up1"))
d_up2 <- mk(up2_tbl$base, up2_tbl$min, up2_tbl$max, up2_tbl$log, "up2")
d_up3 <- with(cen_factors[match(up3_free, cen_factors$name), ],
              mk(name, lo, hi, log, "up3"))
d_up4 <- with(spray_factors[match(up4_free, spray_factors$name), ],
              mk(name, min, max, log, "up4"))

factors <- rbind(d_up1, d_up2, d_up3, d_up4)
factors$name <- paste(factors$stage, factors$base, sep = ".")

# Screening ranges: a handful of the widest stage ranges are narrowed toward
# the dynamic core so the sampled hypercube stays in a physically valid regime
# and multi-factor corner failures (NA) drop. These are the DoE screening
# bounds, not the full operating envelope; widen once the sweep is trusted.
narrow <- list(
  "up1.v_tip"         = c(2,    22),   "up1.P_mix"         = c(0.5,  3.5),
  "up1.C_surfactant"  = c(2e-4, 8e-3), "up1.MW_surfactant" = c(400,  6000),
  "up3.pop_frac"      = c(0.25, 0.85), "up4.ALR"           = c(1.5,  8),
  "up4.mdot_gas_dry"  = c(0.15, 0.8),  "up4.T_system"      = c(350,  460))
for (nm in names(narrow)) {
  i <- match(nm, factors$name)
  if (!is.na(i)) { factors$min[i] <- narrow[[nm]][1]; factors$max[i] <- narrow[[nm]][2] }
}
k <- nrow(factors)

# Emit the screening-range table (factor, stage/UP, min, max, log) as a CSV.
UP <- c(up1 = "UP1 mixer", up2 = "UP2 foam-wash", up3 = "UP3 separator",
        up4 = "UP4 dryer")
range_tbl <- data.frame(factor = factors$base, UP = UP[factors$stage],
                        min = factors$min, max = factors$max, log = factors$log,
                        row.names = NULL)
write.csv(range_tbl, file.path(out_dir, "full_train_factor_ranges.csv"),
          row.names = FALSE)
cat(sprintf("Merged dictionary: %d factors (UP1 %d + UP2 %d + UP3 %d + UP4 %d)\n",
            k, nrow(d_up1), nrow(d_up2), nrow(d_up3), nrow(d_up4)))

# -----------------------------------------------------------------------------
# 2. FLAT, NA-GUARDED TRAIN MODEL:  model(x named by factors$name) -> features
# -----------------------------------------------------------------------------
feature_names <- c("porosity", "sphericity", "size_um", "skin", "tapped")
na_row <- setNames(rep(NA_real_, length(feature_names)), feature_names)

split_by <- function(x, stage) {
  idx <- factors$stage == stage
  setNames(as.numeric(x[factors$name[idx]]), factors$base[idx])
}

run_train_x <- function(x) {
  tryCatch({
    mixer_x <- mixer_nominal_x
    u1 <- split_by(x, "up1"); for (nm in names(u1)) mixer_x[nm] <- u1[[nm]]
    cen_op <- cen_nominal
    u3 <- split_by(x, "up3"); for (nm in names(u3)) cen_op[nm] <- u3[[nm]]
    spray_op <- sp_mid
    u4 <- split_by(x, "up4"); for (nm in names(u4)) spray_op[nm] <- u4[[nm]]
    wash_pars <- as.list(split_by(x, "up2"))

    r <- NULL
    # swallow deSolve/lsoda solver chatter (DINTDY etc.) on corner failures;
    # the NA-guard below is what actually records them
    suppressWarnings(suppressMessages(invisible(capture.output(
      r <- run_full_train(mixer_x = mixer_x, cen_op = cen_op, spray_op = spray_op,
                          wash_pars = wash_pars, up2 = "ode")))))
    if (is.null(r)) return(na_row)
    sp <- r$spray
    v <- c(porosity   = sp[["phi_porosity_z"]],
           sphericity = sp[["Omega_struct_z"]],
           size_um    = sp[["D_particle_um"]],
           skin       = sp[["theta_skin_z"]],
           tapped     = sp[["rho_tapped"]])
    if (any(!is.finite(v))) return(na_row)
    v
  }, error = function(e) na_row)
}

# -----------------------------------------------------------------------------
# 3. MORRIS OAT DESIGN + ELEMENTARY EFFECTS  (base-R, from unified_model.R)
# -----------------------------------------------------------------------------
delta <- levels / (2 * (levels - 1))

scale_design <- function(X01) {
  X <- X01
  for (j in seq_len(k)) {
    if (factors$log[j]) {
      X[, j] <- exp(log(factors$min[j]) +
                    X01[, j] * (log(factors$max[j]) - log(factors$min[j])))
    } else {
      X[, j] <- factors$min[j] + X01[, j] * (factors$max[j] - factors$min[j])
    }
  }
  colnames(X) <- factors$name
  X
}

morris_oat_design <- function(k, r, levels, delta) {
  grid <- seq(0, 1 - delta, length.out = levels - levels / 2)
  X <- matrix(NA_real_, nrow = r * (k + 1), ncol = k)
  info <- vector("list", r)
  row <- 1
  for (t in seq_len(r)) {
    base <- sample(grid, k, replace = TRUE)
    ord  <- sample.int(k)
    dirs <- numeric(k)
    x <- base
    X[row, ] <- x
    for (s in seq_len(k)) {
      j <- ord[s]
      d <- sample(c(-delta, delta), 1)
      if (x[j] + d < 0 || x[j] + d > 1) d <- -d
      x[j] <- x[j] + d
      dirs[s] <- d
      X[row + s, ] <- x
    }
    info[[t]] <- list(order = ord, dirs = dirs)
    row <- row + k + 1
  }
  list(X = X, info = info)
}

elementary_effects <- function(design, Y, k, r) {
  n_out <- ncol(Y)
  ee <- lapply(seq_len(n_out), function(i) matrix(NA_real_, r, k))
  for (t in seq_len(r)) {
    off <- (t - 1) * (k + 1)
    ord <- design$info[[t]]$order
    dirs <- design$info[[t]]$dirs
    for (s in seq_len(k)) {
      j <- ord[s]
      dy <- Y[off + s + 1, ] - Y[off + s, ]
      for (i in seq_len(n_out)) ee[[i]][t, j] <- dy[i] / dirs[s]
    }
  }
  ee
}

# -----------------------------------------------------------------------------
# 4. RUN THE SCREEN
# -----------------------------------------------------------------------------
cache_file <- file.path(out_dir, sprintf("full_train_morris_r%d.rds", r_traj))
if (file.exists(cache_file)) {
  cat("Loading cached runs from", cache_file, "\n")
  cached <- readRDS(cache_file); design <- cached$design; Y <- cached$Y
} else {
  design <- morris_oat_design(k, r_traj, levels, delta)
  Xphys  <- scale_design(design$X)
  n_runs <- nrow(Xphys)
  cat(sprintf("Screen: %d factors, %d trajectories -> %d full-train runs...\n",
              k, r_traj, n_runs))
  rows <- lapply(seq_len(n_runs), function(i) Xphys[i, ])
  t0 <- Sys.time()
  if (.Platform$OS.type == "unix") {
    cores <- max(1, parallel::detectCores() - 1)
    ylist <- parallel::mclapply(rows, run_train_x, mc.cores = cores)
  } else {
    ylist <- lapply(rows, run_train_x)
  }
  Y <- do.call(rbind, ylist)
  cat(sprintf("  done in %.1f s;  %d/%d runs valid (%.1f%% NA)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(stats::complete.cases(Y)), n_runs,
              100 * mean(!stats::complete.cases(Y))))
  saveRDS(list(design = design, Y = Y), cache_file)
}

# -----------------------------------------------------------------------------
# 5. INDICES + NA DIAGNOSTICS
# -----------------------------------------------------------------------------
ee_list <- elementary_effects(design, Y, k, r_traj)
names(ee_list) <- feature_names

morris_stats <- function(ee)
  data.frame(factor  = factors$name,
             stage   = factors$stage,
             mu.star = colMeans(abs(ee), na.rm = TRUE),
             sigma   = apply(ee, 2, sd, na.rm = TRUE),
             n_na    = colSums(is.na(ee)))
stats_list <- setNames(lapply(ee_list, morris_stats), feature_names)

all_stats <- do.call(rbind, lapply(feature_names, function(o)
  cbind(feature = o, stats_list[[o]])))
write.csv(all_stats, file.path(out_dir, "full_train_morris_indices.csv"),
          row.names = FALSE)

# per-factor NA rate (share of elementary-effect pairs that came back NA)
na_by_factor <- colMeans(is.na(ee_list[[1]]))
na_tbl <- data.frame(factor = factors$name, na_rate = round(na_by_factor, 3))
na_tbl <- na_tbl[order(-na_tbl$na_rate), ]
write.csv(na_tbl, file.path(out_dir, "full_train_na_by_factor.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. CONSOLE SUMMARY
# -----------------------------------------------------------------------------
cat("\nTop-8 drivers per feature (mu*, [stage]):\n")
for (o in feature_names) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat(sprintf("  %-11s : %s\n", o,
      paste(sprintf("%s", sub("^(up[0-9])\\.", "[\\1]", head(st$factor, 8))),
            collapse = ", ")))
}
# flag genuinely NA-prone factors: well above the run-level rate AND above the
# r-quantization floor (1/r), so a single unlucky trajectory doesn't get flagged
overall_na <- mean(!stats::complete.cases(Y))
na_thresh  <- max(1.5 / r_traj, 3 * overall_na)
hi_na <- na_tbl[na_tbl$na_rate > na_thresh, ]
if (nrow(hi_na) > 0) {
  cat(sprintf("\nNA-prone factors (>%.0f%%; narrow these ranges next):\n", 100*na_thresh))
  for (i in seq_len(nrow(hi_na)))
    cat(sprintf("  %-24s %.0f%%\n", hi_na$factor[i], 100 * hi_na$na_rate[i]))
} else cat(sprintf("\nNo factor exceeds the %.0f%% NA threshold (run-level NA %.1f%%).\n",
                   100*na_thresh, 100*overall_na))

# -----------------------------------------------------------------------------
# 7. FEATURE LENS PLOT (mu*-sigma, one panel per feature of interest)
# -----------------------------------------------------------------------------
stage_col <- c(up1 = "#2c7fb8", up2 = "#7b3294", up3 = "#d95f0e", up4 = "#1a9850")
png(file.path(out_dir, "full_train_morris_features.png"),
    width = 2600, height = 1700, res = 150)
op <- par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.6, 1))
for (o in feature_names) {
  st <- stats_list[[o]]
  plot(st$mu.star, st$sigma, pch = 19, col = stage_col[st$stage], cex = 1.2,
       xlab = expression(mu * "*"), ylab = expression(sigma), main = o,
       xlim = c(0, max(st$mu.star, na.rm = TRUE) * 1.15 + 1e-9))
  abline(0, 1, lty = 2, col = "grey60")
  top <- order(-st$mu.star)[seq_len(min(6, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], sub("^up[0-9]\\.", "", st$factor[top]),
       pos = 3, cex = 0.62, xpd = NA)
}
plot.new(); legend("center", legend = names(stage_col), pch = 19,
                   col = stage_col, title = "stage", bty = "n", cex = 1.2)
par(op); invisible(dev.off())

cat("\nWrote to", out_dir, ":\n")
cat("  full_train_morris_indices.csv, full_train_na_by_factor.csv\n")
cat("  full_train_morris_features.png\n")
