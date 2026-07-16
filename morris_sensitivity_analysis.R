#!/usr/bin/env Rscript
# =============================================================================
# Morris (Elementary Effects) sensitivity analysis of the reduced-order
# two-fluid (internal-mix) atomization + drying model described in
# "deepresearchreport.md" (Feed Properties -> Two-Phase Conditioning ->
# Nozzle Hydraulics -> Primary Atomization -> Secondary Breakup ->
# Drying / Particle Formation).
#
# Outputs screened (naming follows the "Latex Coagulation Engine: Master
# Nomenclature" sheet where a symbol exists):
#   1. D_particle     - final (dry) particle size                    [um]
#   2. theta_skin_z   - skin network fraction (skin formation)       [-]
#   3. Omega_struct_z - structural memory / sphericity state         [-]
#   4. phi_porosity_z - total porosity (void fraction state)         [-]
#   5. rho_tapped     - powder tapped density (bulk analogue of
#                       rho_colloid_out / SG_out)                    [kg/m3]
#
# Method: Morris one-at-a-time screening. Uses the `sensitivity` package
# (morris()) when installed; otherwise falls back to an internal base-R
# implementation of the same OAT trajectory design, so the script has no
# hard package dependencies. The Morris plot (mu* vs sigma) is drawn with
# base graphics.
#
# Run:  Rscript morris_sensitivity_analysis.R
# =============================================================================

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Input factors and ranges
# -----------------------------------------------------------------------------
# Symbols follow the nomenclature sheet where available:
#   C_solid_mass : solid mass fraction in feed (wt/wt)
#   alpha_g_0    : total gas holdup of the feed (foam quality)      [-]
#   T_system     : dryer gas temperature                            [K]
#   P_system     : atomizing air (gas) supply pressure              [Pa]
# Others use the report's spray notation (sigma, mu_L, ALR, D_b, ...).
# Ranges are taken from the report's stated validity ranges / worked example.
factors <- data.frame(
  name = c("ALR",          # air-liquid mass ratio  m_G/m_L            [-]
           "P_system",     # atomizing air pressure                    [Pa]
           "mdot_L",       # liquid feed mass flow                     [kg/s]
           "sigma",        # liquid surface tension                    [N/m]
           "mu_L",         # liquid (continuous phase) viscosity       [Pa s]
           "rho_L",        # liquid (slurry) density                   [kg/m3]
           "alpha_g_0",    # feed foam quality / gas holdup            [-]
           "D_b",          # feed bubble diameter                      [m]
           "C_solid_mass", # solid mass fraction in feed               [-]
           "T_system",     # dryer gas temperature                     [K]
           "RH_gas"),      # dryer gas relative humidity               [-]
  min  = c( 1.0, 2.0e5, 0.002, 0.030, 0.0012, 1000, 0.05, 2.0e-5, 0.05, 330, 0.01),
  max  = c(10.0, 7.0e5, 0.020, 0.070, 0.0560, 1300, 0.60, 2.0e-4, 0.40, 470, 0.30),
  # sample wide-ranging positive factors log-uniformly
  log  = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
k <- nrow(factors)

# Map the unit-hypercube Morris design onto physical values
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

# -----------------------------------------------------------------------------
# 2. Built-in reduced-order model (one row of physical inputs -> 5 outputs)
# -----------------------------------------------------------------------------
# Fixed constants / nozzle geometry (report worked example)
P_atm   <- 1.013e5   # ambient pressure                     [Pa]
T_feed  <- 300       # feed / ambient temperature           [K]
R_air   <- 287       # gas constant, air                    [J/kg K]
gamma_a <- 1.4       # heat capacity ratio, air             [-]
A_L     <- 1.0e-6    # liquid passage area                  [m2]
A_G     <- 5.0e-6    # gas passage area                     [m2]
D_h     <- 1.0e-3    # nozzle hydraulic diameter            [m]
rho_s   <- 1400      # dry solid (polymer) density          [kg/m3]
a_prim  <- 1.0e-7    # primary colloid particle radius      [m]
kB      <- 1.380649e-23

spray_dry_model <- function(x) {
  ALR    <- x[["ALR"]];       P_G   <- x[["P_system"]]
  mdot_L <- x[["mdot_L"]];    sigma <- x[["sigma"]]
  mu_L   <- x[["mu_L"]];      rho_L <- x[["rho_L"]]
  alpha0 <- x[["alpha_g_0"]]; D_b   <- x[["D_b"]]
  C_sol  <- x[["C_solid_mass"]]
  T_gas  <- x[["T_system"]];  RH    <- x[["RH_gas"]]

  ## --- Module 1: Feed properties -------------------------------------------
  rho_G0   <- P_G / (R_air * T_feed)              # gas density at feed P
  rho_eff0 <- (1 - alpha0) * rho_L + alpha0 * rho_G0   # rho_eff,z (HEM)
  mu_eff   <- mu_L * (1 + 2.5 * alpha0)           # mu_eff,z (dilute-bubble)
  sigma_eff<- sigma * (1 - 0.5 * alpha0)          # foam-reduced surface tension

  ## --- Module 2: Two-phase conditioning (Boyle expansion to exit) ----------
  r_exp   <- P_G / P_atm                          # expansion ratio
  alpha_e <- alpha0 * r_exp / (1 - alpha0 + alpha0 * r_exp)  # exact volume rule
  alpha_e <- min(alpha_e, 0.97)
  D_b_e   <- D_b * r_exp^(1/3)                    # bubble growth D_b ~ P^-1/3
  rho_Ge  <- P_atm / (R_air * T_feed)
  rho_eff <- (1 - alpha_e) * rho_L + alpha_e * rho_Ge

  ## --- Module 3: Nozzle hydraulics ------------------------------------------
  U_L <- mdot_L / (rho_L * A_L)                   # liquid exit velocity
  # gas: fully-expanded isentropic velocity (choked upstream if P_G/P_atm>1.89)
  T_exit <- T_gas * 0 + T_feed                    # air enters at feed temp
  U_G <- 0.9 * sqrt(pmax(2 * gamma_a / (gamma_a - 1) * R_air * T_exit *
                         (1 - (P_atm / P_G)^((gamma_a - 1) / gamma_a)), 0))
  U_rel <- max(U_G - U_L, 5)
  rho_A <- P_atm / (R_air * T_feed)               # ambient air density

  ## --- Module 4: Primary atomization (Rizk & Lefebvre plain-jet airblast) --
  SMD <- 0.48 * D_h * (sigma_eff / (rho_A * U_rel^2 * D_h))^0.4 *
           (1 + 1 / ALR)^0.4 +
         0.15 * D_h * sqrt(mu_eff^2 / (sigma_eff * rho_L * D_h)) *
           (1 + 1 / ALR)

  ## --- Module 5: Secondary breakup (Weber-number correction) ---------------
  We_g <- rho_A * U_rel^2 * SMD / sigma_eff       # gas Weber number of drop
  f_sec <- 0.55 + 0.45 / (1 + We_g / 80)          # ~20-45 % diameter reduction
  d_drop <- SMD * f_sec                           # droplet Dv50 entering dryer

  ## --- Module 6: Drying / particle formation -------------------------------
  # d2-law evaporation coefficient, scaled by dryer driving force
  kappa <- 5e-8 * pmax(T_gas - T_feed, 1) / 100 * (1 - RH)      # [m2/s]
  # Stokes-Einstein diffusivity of primary colloid particles
  D_diff <- kB * T_feed / (6 * pi * mu_L * a_prim)              # [m2/s]
  Pe <- kappa / (8 * D_diff)                      # drying Peclet number

  # theta_skin,z : skin network fraction, driven by Pe x solids loading
  S_skin <- Pe * C_sol
  theta_skin <- S_skin / (S_skin + 500)

  # alpha_trap,z : gas retained inside droplets - bubbles comparable to or
  # larger than the droplet are shed during breakup, small ones survive
  alpha_trap <- alpha_e / (1 + D_b_e / d_drop)
  # drying-induced (skin-locked) porosity
  phi_dry <- 0.25 * theta_skin * (1 - exp(-Pe / 2000))
  # phi_porosity,z : total porosity (void fraction state)
  phi_porosity <- alpha_trap + (1 - alpha_trap) * phi_dry

  # Omega_struct,z : sphericity; high sigma & mu resist skin buckling/collapse
  resist <- (sigma_eff / 0.05) * (mu_eff / 0.01)^0.3
  Omega_struct <- 1 - 0.55 * theta_skin / (1 + resist)

  # Final particle size: solids mass balance on one droplet
  # solid volume fraction of droplet liquid phase (phi_s in nomenclature)
  D_particle <- d_drop * ((1 - alpha_trap) * rho_L * C_sol /
                          (rho_s * (1 - phi_porosity)))^(1/3)

  # Envelope density and tapped density (bulk analogue of rho_colloid,out)
  rho_env  <- (1 - phi_porosity) * rho_s
  f_pack   <- 0.64 * (0.55 + 0.45 * Omega_struct) *
              (D_particle / (D_particle + 8e-6))   # cohesion penalty for fines
  rho_tapped <- rho_env * f_pack

  c(D_particle_um = D_particle * 1e6,
    theta_skin_z  = theta_skin,
    Omega_struct_z = Omega_struct,
    phi_porosity_z = phi_porosity,
    rho_tapped    = rho_tapped)
}

run_model <- function(X01) {
  Xphys <- scale_design(as.matrix(X01))
  t(apply(Xphys, 1, spray_dry_model))
}

# -----------------------------------------------------------------------------
# 3. Morris design: `sensitivity` package if available, else base-R fallback
# -----------------------------------------------------------------------------
r_traj <- 60   # number of Morris trajectories
levels <- 8    # grid levels
delta  <- levels / (2 * (levels - 1))   # standard Morris step

morris_oat_design <- function(k, r, levels, delta) {
  grid <- seq(0, 1 - delta, length.out = levels - levels / 2)  # feasible bases
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
  # Y: matrix, one column per output; returns list of (r x k) EE matrices
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

use_pkg <- requireNamespace("sensitivity", quietly = TRUE)

if (use_pkg) {
  message("Using sensitivity::morris()")
  mor <- sensitivity::morris(model = NULL, factors = factors$name, r = r_traj,
                             design = list(type = "oat", levels = levels,
                                           grid.jump = levels / 2),
                             binf = 0, bsup = 1)
  Y <- run_model(mor$X)
  ee_list <- lapply(seq_len(ncol(Y)), function(i) {
    m <- mor
    sensitivity::tell(m, Y[, i])
    m$ee
  })
} else {
  message("Package 'sensitivity' not found - using built-in Morris OAT design")
  design <- morris_oat_design(k, r_traj, levels, delta)
  Y <- run_model(design$X)
  ee_list <- elementary_effects(design, Y, k, r_traj)
}

outputs <- colnames(Y)

morris_stats <- function(ee) {
  data.frame(factor  = factors$name,
             mu      = colMeans(ee),
             mu.star = colMeans(abs(ee)),
             sigma   = apply(ee, 2, sd))
}
stats_list <- setNames(lapply(ee_list, morris_stats), outputs)

# -----------------------------------------------------------------------------
# 4. Morris plots (mu* vs sigma), one panel per output
# -----------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)

titles <- c(D_particle_um  = "Final particle size  D_particle [um]",
            theta_skin_z   = "Skin formation  theta_skin,z [-]",
            Omega_struct_z = "Sphericity  Omega_struct,z [-]",
            phi_porosity_z = "Porosity  phi_porosity,z [-]",
            rho_tapped     = "Tapped density  rho_tapped [kg/m3]")

plot_morris_panel <- function(st, title) {
  xr <- range(st$mu.star); pad <- 0.12 * diff(xr); if (pad == 0) pad <- 1
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25 + pad * 0.01),
       ylim = c(0, max(st$sigma) * 1.30),
       pch = 21, bg = "steelblue", cex = 1.3,
       xlab = expression(mu * "*  (mean |elementary effect|)"),
       ylab = expression(sigma * "  (sd of elementary effects)"),
       main = title, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey55")           # sigma = mu* (non-linear)
  abline(0, 0.1, lty = 3, col = "grey70")         # sigma = 0.1 mu* (~linear)
  text(st$mu.star, st$sigma, labels = st$factor, pos = 3, cex = 0.72,
       offset = 0.35, xpd = NA)
}

png(file.path("output", "morris_sensitivity_plots.png"),
    width = 2200, height = 1500, res = 160)
op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2.5, 1))
for (o in outputs) plot_morris_panel(stats_list[[o]], titles[[o]])
plot.new()
legend("center", bty = "n", cex = 1.0, title = "Morris screening",
       legend = c(sprintf("%d trajectories, %d levels", r_traj, levels),
                  sprintf("%d model runs", nrow(Y)),
                  "dashed: sigma = mu* (interaction / non-linear)",
                  "dotted: sigma = 0.1 mu* (near-linear)"))
par(op); dev.off()

# -----------------------------------------------------------------------------
# 5. Console summary + CSV export
# -----------------------------------------------------------------------------
all_stats <- do.call(rbind, lapply(outputs, function(o)
  cbind(output = o, stats_list[[o]])))
write.csv(all_stats, file.path("output", "morris_indices.csv"),
          row.names = FALSE)

for (o in outputs) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat("\n==", titles[[o]], "==\n")
  print(st, row.names = FALSE, digits = 3)
}
cat("\nWrote output/morris_sensitivity_plots.png and output/morris_indices.csv\n")
