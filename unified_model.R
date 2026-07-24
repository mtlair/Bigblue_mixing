#!/usr/bin/env Rscript
# =============================================================================
# UNIFIED PROCESS MODEL: mixer (UP1) -> ... -> spray dryer (UP2)
# =============================================================================
# Connects the two existing modules end-to-end through a common feed-stream
# interface, with two identity placeholder stages reserved for the unit
# operations to be added later:
#
#   UP1 gassed/templated mixer      unified/up1_mixer_module.R
#     -> intermediate stage 1       unified/interface_stream.R  (placeholder)
#     -> intermediate stage 2       unified/interface_stream.R  (placeholder)
#       -> UP2 spray dryer          unified/up2_spray_dryer_module.R
#
# The stream carries composition, physical state, particulate/gas/template
# state and structure history from the mixer exit into the dryer feed, so
# the factors the two standalone screens duplicated (solids, additives,
# surface chemistry, template) are now defined ONCE, at the mixer, and the
# dryer sees their transformed values (dilution, gas holdup, wet skin,
# residual template fraction) instead of independent knobs.
#
# What this script does:
#   1. Nominal end-to-end run for each templating strategy (rigid / gas /
#      surface_weld / capillary_bridge), printing the stream at each
#      interface -> unified_output/nominal_chain_summary.txt
#   2. Morris elementary-effects screen over the FULL chain (41 factors:
#      24 mixer-side + 17 dryer-side) for the configured template type.
#      Base-R OAT implementation (no 'sensitivity' dependency).
#   3. Zone-classified mu*-sigma lens plots (Process / Surface Chemistry /
#      Polymer) and a CSV of all indices.
#
# Run from the repo root:  Rscript unified_model.R
# Hard dependency: deSolve (UP1 ODE integration). Everything else is base R.
# =============================================================================

source("unified/up1_mixer_module.R")
source("unified/interface_stream.R")
source("unified/up2_spray_dryer_module.R")

set.seed(42)

# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
TEMPLATE_TYPE <- 4      # Morris screen template strategy: 4 = capillary_bridge
                        # (immiscible liquid template == UP2's pore-templating
                        # emulsion; 3 = surface_weld maximizes the RTF risk path)
r_traj  <- 30           # Morris trajectories
levels  <- 8            # Morris grid levels
out_dir <- "unified_output"
dir.create(out_dir, showWarnings = FALSE)

equipment <- up1_default_equipment()

# -----------------------------------------------------------------------------
# 2. UNIFIED FACTOR DICTIONARY
# -----------------------------------------------------------------------------
# module: which stage consumes the factor directly.
# group : lens for the grouped plots (Process / Surface / Polymer).
# log   : sampled log-uniformly.
# Name collisions between the standalone screens are resolved with explicit
# names: P_mix/T_mix/tau_mix (mixer vessel) vs P_atom_air/T_dryer_in (dryer).
# Factors the dryer used to own that now COME FROM THE STREAM (alpha_g_0,
# C_solid_mass, rho_L, T_feed, phi_emulsion, D_template, additive and
# surface-chemistry concentrations) are mixer-side factors here.
fac <- function(name, min, max, log, unit, module, group, desc)
  data.frame(name = name, min = min, max = max, log = log, unit = unit,
             module = module, group = group, desc = desc,
             stringsAsFactors = FALSE)

factors <- rbind(
  # --- UP1 mixer side --------------------------------------------------------
  fac("v_tip",          0.5,    32.0,  FALSE, "m/s",   "up1", "Process", "impeller tip speed"),
  fac("tau_mix",        5.0,    40.0,  FALSE, "min",   "up1", "Process", "mixer residence time"),
  fac("P_mix",          5/14.7, 64.7/14.7, FALSE, "atm", "up1", "Process", "mixer headspace pressure"),
  fac("T_mix",          290,    360,   FALSE, "K",     "up1", "Process", "mixer temperature (-> dryer feed T)"),
  fac("C_gas_diss_temp",0.0,    0.25,  FALSE, "-",     "up1", "Process", "dissolved gas Henry-normalised saturation index (model units; C_eq ~ 0.2 at 1 atm)"),
  fac("Q_gas",          0.4,    2.0,   FALSE, "-",     "up1", "Process", "sparge gas flow"),
  fac("Q_colloid",      40/60,  200/60,FALSE, "-",     "up1", "Process", "colloid feed flow"),
  fac("Q_template",     0.1,    2.0,   FALSE, "-",     "up1", "Process", "template feed flow"),
  fac("C_solid_mass",   0.10,   0.60,  FALSE, "wt/wt", "up1", "Polymer", "feed solids (diluted in-tank)"),
  fac("C_temp_mass",    0.01,   0.10,  FALSE, "wt/wt", "up1", "Polymer", "rigid template mass fraction"),
  fac("D_particle",     0.5,    2.0,   FALSE, "um",    "up1", "Polymer", "primary colloid particle diameter"),
  fac("D_template",     0.2,    1.5,   FALSE, "um",    "up1", "Polymer", "template droplet/seed diameter"),
  fac("C_binder",       0.01,   0.15,  FALSE, "wt/wt", "up1", "Polymer", "binder concentration"),
  fac("C_monomer",      0.001,  0.03,  FALSE, "wt/wt", "up1", "Polymer", "residual monomer concentration"),
  fac("C_plasticizer",  0.005,  0.04,  FALSE, "wt/wt", "up1", "Polymer", "plasticizer concentration"),
  fac("C_surfactant",   5e-5,   0.015, TRUE,  "wt/wt", "up1", "Surface", "surfactant concentration"),
  fac("CMC",            0.001,  0.015, FALSE, "wt/wt", "up1", "Surface", "critical micelle concentration"),
  fac("HLB",            8.0,    16.0,  FALSE, "-",     "up1", "Surface", "surfactant HLB"),
  fac("MW_surfactant",  200,    10000, TRUE,  "g/mol", "up1", "Surface", "surfactant molecular weight"),
  fac("A_molecule",     0.2,    12.5,  FALSE, "nm2",   "up1", "Surface", "surfactant molecule area"),
  fac("D_surf",         1e-11,  1e-9,  TRUE,  "m2/s",  "up1", "Surface", "surfactant diffusivity"),
  fac("ionic_strength", 0.01,   0.50,  FALSE, "M",     "up1", "Surface", "ionic strength"),
  fac("Delta_pH",       0.0,    5.0,   FALSE, "pH",    "up1", "Surface", "delta pH vs isoelectric point"),
  fac("template_dose",  0.02,   1.10,  FALSE, "-",     "up1", "Polymer", "template liquid / interstitial void"),
  fac("chi_template",   0.30,   2.50,  FALSE, "-",     "up1", "Polymer", "Flory-Huggins chi (template solvent vs polymer; <0.5 good, 0.5 theta, >0.5 poor)"),
  # --- UP2 dryer side --------------------------------------------------------
  fac("ALR",            0.9,    1.8,   FALSE, "-",     "up2", "Process", "atomizing air-liquid mass ratio (measured up4atom_scfm/up4_feed, visc.xlsx: 0.90-1.78)"),
  fac("P_atom_air",     1.2e5,  1.8e5, FALSE, "Pa",    "up2", "Process", "atomizing air supply pressure (measured max_up4atom_psig, visc.xlsx: 2.7-11.4 psig abs)"),
  fac("P_feed",         1.5e5,  1.0e6, FALSE, "Pa",    "up2", "Process", "feed pump / hold-line pressure"),
  fac("mdot_L",         0.013,  0.020, FALSE, "kg/s",  "up2", "Process", "liquid feed mass flow to nozzle (measured up4_feed, visc.xlsx: 0.0135-0.0194 kg/s)"),
  fac("sigma",          0.030,  0.070, FALSE, "N/m",   "up2", "Surface", "liquid surface tension"),
  fac("mu_L",           0.0012, 0.0560,TRUE,  "Pa s",  "up2", "Polymer", "serum (continuous phase) viscosity"),
  fac("T_dryer_in",     410,    440,   FALSE, "K",     "up2", "Process", "dryer gas inlet temperature (measured up4_Tin, visc.xlsx: 151.6-153.5 C = 424.8-426.7 K)"),
  fac("mdot_gas_dry",   0.40,   0.62,  TRUE,  "kg/s",  "up2", "Process", "dryer gas mass flow (evaporative estimate from up4_feed + solids to <=0.5% moisture, Tin/Tout balance: 0.40-0.62 kg/s)"),
  fac("Y_in",           0.001,  0.020, FALSE, "kg/kg", "up2", "Process", "dryer inlet absolute humidity"),
  fac("t_hold",         5,      600,   TRUE,  "s",     "up2", "Process", "pressurized hold before nozzle"),
  fac("D_b",            2.0e-5, 2.0e-4,TRUE,  "m",     "up2", "Process", "feed bubble diameter (until the transfer-line stage sets it)"),
  fac("Tg_polymer",     280,    380,   FALSE, "K",     "up2", "Polymer", "dry-polymer glass transition"),
  fac("n_flow",         0.20,   0.67,  FALSE, "-",     "up2", "Polymer", "power-law flow index (measured post-UP1 slurry envelope, visc.xlsx: 0.20-0.67, mean 0.44)"),
  fac("k_perm_mono",    5.0,    60.0,  FALSE, "-",     "up2", "Polymer", "monomer permeation coefficient"),
  fac("k_perm_plast",   5.0,    60.0,  FALSE, "-",     "up2", "Polymer", "plasticizer permeation coefficient"),
  fac("k_perm_bind",    2.0,    40.0,  FALSE, "-",     "up2", "Polymer", "binder pore-blocking coefficient"),
  fac("T_bp_solv",      300,    360,   FALSE, "K",     "up2", "Polymer", "template solvent boiling point")
)
k <- nrow(factors)

up1_names <- factors$name[factors$module == "up1"]
up2_names <- factors$name[factors$module == "up2"]

# Unified name -> UP1 module's native name
up1_rename <- c(tau_mix = "tau", P_mix = "P_system", T_mix = "T_system")

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

nominal_x <- setNames(ifelse(factors$log,
                             sqrt(factors$min * factors$max),
                             (factors$min + factors$max) / 2),
                      factors$name)

# Override nominal to 25% solid (the normal operating point; factor midpoint
# would be 35% but calibration data and process target are at 20-25%).
nominal_x["C_solid_mass"] <- 0.25

# -----------------------------------------------------------------------------
# 3. THE UNIFIED CHAIN
# -----------------------------------------------------------------------------
up1_pars_from_x <- function(x) {
  p <- as.list(x[up1_names])
  for (un in names(up1_rename)) {
    p[[up1_rename[[un]]]] <- p[[un]]
    p[[un]] <- NULL
  }
  p
}

# Screened chain outputs: mixer diagnostics (up1_) + product properties (up2_)
up1_keep <- c("Mixing_Potential", "Blended_Porosity", "Blended_Size_um",
              "Blended_WetSkin", "Template_Target_Score",
              "Residual_Template_Fraction")
unified_output_names <- c(paste0("up1_", up1_keep),
                          paste0("up2_", up2_output_names))
n_out <- length(unified_output_names)

run_unified <- function(x, template_type = TEMPLATE_TYPE, verbose = FALSE) {
  res <- tryCatch({
    eq <- equipment; eq$template_type <- template_type
    p1 <- up1_pars_from_x(x)
    r1 <- up1_run_mixer(p1, eq)
    if (any(is.na(r1$outputs))) return(setNames(rep(NA_real_, n_out), unified_output_names))

    stream <- stream_from_up1(r1, p1, eq)
    if (verbose) print_stream(stream, "stream: mixer exit")
    stream <- intermediate_stage_1(stream)   # placeholder (transfer line / pump)
    stream <- intermediate_stage_2(stream)   # placeholder (hold / pre-conditioning)
    if (verbose) print_stream(stream, "stream: dryer feed")

    x_up2 <- as.list(x[up2_names])
    # Forward the (muted-by-default) overlay knobs if the caller set them;
    # they are not screened factors, so they are not in up2_names.
    if ("size_template" %in% names(x))    x_up2[["size_template"]]    <- x[["size_template"]]
    if ("morphology_recal" %in% names(x)) x_up2[["morphology_recal"]] <- x[["morphology_recal"]]
    r2 <- up2_run_dryer(stream, x_up2)
    setNames(c(r1$outputs[up1_keep], r2), unified_output_names)
  }, error = function(e) setNames(rep(NA_real_, n_out), unified_output_names))
  res
}

# -----------------------------------------------------------------------------
# 4. NOMINAL END-TO-END DEMONSTRATION (all four templating strategies)
# -----------------------------------------------------------------------------
template_modes <- c(rigid = 1, gas = 2, surface_weld = 3, capillary_bridge = 4)

sink(file.path(out_dir, "nominal_chain_summary.txt"))
cat("UNIFIED CHAIN - nominal factor values, per templating strategy\n")
cat(sprintf("(mixer -> 2 placeholder stages -> spray dryer; %d factors)\n\n", k))
nom_tbl <- NULL
for (mn in names(template_modes)) {
  cat(sprintf("================ template_type = %d (%s) ================\n",
              template_modes[[mn]], mn))
  y <- run_unified(nominal_x, template_type = template_modes[[mn]],
                   verbose = (mn == "capillary_bridge"))
  for (o in unified_output_names)
    cat(sprintf("  %-32s %s\n", o, format(signif(y[[o]], 4))))
  cat("\n")
  nom_tbl <- rbind(nom_tbl, data.frame(template = mn, t(y)))
}
sink()
write.csv(nom_tbl, file.path(out_dir, "nominal_chain_outputs.csv"), row.names = FALSE)
cat("Wrote", file.path(out_dir, "nominal_chain_summary.txt"), "\n")

# -----------------------------------------------------------------------------
# 5. MORRIS SCREEN OVER THE FULL CHAIN (base-R OAT design)
# -----------------------------------------------------------------------------
delta <- levels / (2 * (levels - 1))

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

cache_file <- file.path(out_dir, sprintf("unified_morris_t%d_r%d.rds",
                                         TEMPLATE_TYPE, r_traj))
cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
# Invalidate a stale cache whose output schema no longer matches the current
# module (e.g. after up2_output_names grew with Module 0f): the cache key only
# encodes TEMPLATE_TYPE and r_traj, not the output columns, so guard on the
# recorded factor count (k) and output width before trusting it.
if (!is.null(cached) &&
    (ncol(cached$Y) != n_out || nrow(cached$design$X) != r_traj * (k + 1) ||
     ncol(cached$design$X) != k)) {
  cat("Cached Morris runs in", cache_file,
      "are stale (schema changed); recomputing.\n")
  cached <- NULL
}
if (!is.null(cached)) {
  cat("Loading cached Morris runs from", cache_file, "\n")
  design <- cached$design; Y <- cached$Y
} else {
  design <- morris_oat_design(k, r_traj, levels, delta)
  Xphys  <- scale_design(design$X)
  n_runs <- nrow(Xphys)
  cat(sprintf("Morris screen: %d factors, %d trajectories -> %d chain runs...\n",
              k, r_traj, n_runs))
  rows <- lapply(seq_len(n_runs), function(i) Xphys[i, ])
  if (.Platform$OS.type == "unix") {
    cores <- max(1, parallel::detectCores() - 1)
    ylist <- parallel::mclapply(rows, run_unified, mc.cores = cores)
  } else {
    ylist <- lapply(rows, run_unified)
  }
  Y <- do.call(rbind, ylist)
  saveRDS(list(design = design, Y = Y), cache_file)
  cat(sprintf("  %d/%d runs valid; cached to %s\n",
              sum(stats::complete.cases(Y)), n_runs, cache_file))
}

ee_list <- elementary_effects(design, Y, k, r_traj)
names(ee_list) <- unified_output_names

morris_stats <- function(ee) {
  data.frame(factor  = factors$name,
             module  = factors$module,
             group   = factors$group,
             mu      = colMeans(ee, na.rm = TRUE),
             mu.star = colMeans(abs(ee), na.rm = TRUE),
             sigma   = apply(ee, 2, sd, na.rm = TRUE))
}
stats_list <- setNames(lapply(ee_list, morris_stats), unified_output_names)

all_stats <- do.call(rbind, lapply(unified_output_names, function(o)
  cbind(output = o, stats_list[[o]])))
all_stats$sigma_mu_ratio <- all_stats$sigma / pmax(all_stats$mu.star, 1e-12)
all_stats$interaction_type <- cut(all_stats$sigma_mu_ratio,
  breaks = c(0, 0.05, 0.3, 1.0, Inf),
  labels = c("No_interaction", "Linear", "Complex", "Volatile"),
  right = FALSE)
write.csv(all_stats, file.path(out_dir, "unified_morris_indices.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. LENS PLOTS (zone-classified mu*-sigma, one figure per variable group)
# -----------------------------------------------------------------------------
plot_targets <- c("up1_Mixing_Potential", "up1_Blended_Porosity",
                  "up1_Template_Target_Score", "up1_Residual_Template_Fraction",
                  "up2_d_droplet_um", "up2_span", "up2_D_particle_um",
                  "up2_theta_skin_z", "up2_Omega_struct_z", "up2_phi_porosity_z",
                  "up2_rho_tapped", "up2_Tg_product_K", "up2_X_moisture",
                  "up2_D_pore_um", "up2_f_burst_solv", "up2_sigma_y_cake_MPa")

titles <- c(up1_Mixing_Potential      = "UP1 Mixing potential [-]",
            up1_Blended_Porosity      = "UP1 exit gas holdup [-]",
            up1_Template_Target_Score = "UP1 template target score [-]",
            up1_Residual_Template_Fraction = "UP1 residual template fraction [-]",
            up2_d_droplet_um   = "Spray droplet size Dv50 [um]",
            up2_span           = "Droplet span (d90-d10)/d50 [-]",
            up2_D_particle_um  = "Final particle size Dp50 [um]",
            up2_theta_skin_z   = "Skin formation theta_skin [-]",
            up2_Omega_struct_z = "Sphericity Omega_struct [-]",
            up2_phi_porosity_z = "Porosity phi_porosity [-]",
            up2_rho_tapped     = "Tapped density [kg/m3]",
            up2_Tg_product_K   = "Product Tg_eff [K]",
            up2_X_moisture     = "Residual moisture [-]",
            up2_D_pore_um      = "Templated pore size [um]",
            up2_f_burst_solv   = "Micro-explosion severity f_burst [-]",
            up2_sigma_y_cake_MPa = "Cake yield strength [MPa]")

zone_col <- function(st) {
  ratio <- st$sigma / pmax(st$mu.star, 1e-12)
  cls <- findInterval(ratio, c(0.05, 0.3, 1.0)) + 1  # 1..4
  c("#999999", "#2ecc71", "#ffa500", "#e74c3c")[cls]
}

plot_panel <- function(st, title, grp, n_label = 6) {
  cols <- zone_col(st)
  idx  <- st$group == grp
  st   <- st[idx, ]; cols <- cols[idx]
  if (max(st$mu.star) <= 0) { plot.new(); title(main = title, cex.main = 0.85); return(invisible()) }
  pch_mod <- ifelse(st$module == "up1", 21, 24)   # circle=mixer, triangle=dryer
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25),
       ylim = c(0, max(st$sigma) * 1.30 + 1e-12),
       pch = pch_mod, bg = cols, cex = 1.4, col = "black", lwd = 0.7,
       xlab = expression(mu * "*"), ylab = expression(sigma),
       main = title, cex.main = 0.85)
  abline(0, 1, lty = 2, col = "grey55")
  abline(0, 0.3, lty = 3, col = "grey60")
  abline(0, 0.05, lty = 3, col = "grey75")
  top <- order(-st$mu.star)[seq_len(min(n_label, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], labels = st$factor[top], pos = 3,
       cex = 0.65, offset = 0.3, xpd = NA)
}

groups <- c(Process = "Process", Surface = "Surface", Polymer = "Polymer")
for (g in names(groups)) {
  png(file.path(out_dir, sprintf("unified_morris_%s_variables.png", tolower(g))),
      width = 3200, height = 3200, res = 160)
  op <- par(mfrow = c(4, 4), mar = c(4.2, 4.2, 2.6, 1), oma = c(4, 0, 2, 0))
  for (o in plot_targets) plot_panel(stats_list[[o]], titles[[o]], groups[[g]])
  mtext(sprintf("UNIFIED chain Morris screen - %s variables | template_type=%d | r=%d, %d factors, %d runs",
                g, TEMPLATE_TYPE, r_traj, k, nrow(Y)),
        side = 3, outer = TRUE, cex = 0.85, line = 0.3)
  mtext(paste("Marker: circle = mixer-side factor, triangle = dryer-side factor |",
              "Colour: grey = no interaction, green = linear, orange = complex, red = volatile"),
        side = 1, outer = TRUE, cex = 0.70, line = 1.5, font = 3)
  par(op); dev.off()
}

# -----------------------------------------------------------------------------
# 7. CONSOLE SUMMARY
# -----------------------------------------------------------------------------
cat("\nTop-5 chain drivers per key product property (mu*):\n")
for (o in c("up2_D_particle_um", "up2_phi_porosity_z", "up2_rho_tapped",
            "up2_Tg_product_K", "up2_X_moisture", "up2_f_burst_solv")) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat(sprintf("  %-22s : %s\n", o,
      paste(sprintf("%s[%s]", head(st$factor, 5), head(st$module, 5)),
            collapse = ", ")))
}
cat("\nWrote to", out_dir, ":\n")
cat("  nominal_chain_summary.txt / nominal_chain_outputs.csv\n")
cat("  unified_morris_indices.csv\n")
cat("  unified_morris_{process,surface,polymer}_variables.png\n")
