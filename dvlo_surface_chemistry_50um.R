#!/usr/bin/env Rscript
# =============================================================================
# DLVO + SURFACE CHEMISTRY + TEMPLATING PATH TO 50 µm HIGH-POROSITY PARTICLES
# =============================================================================
# Explores the UP1-2-3 (mixer, foam-wash, separator) space to find routes to
# 50 µm particles with limited surface fusion and high porosity WITHOUT relying
# on UP4 atomizer adjustments (those are the known path).
#
# Three independent levers identified from mechanistic analysis:
#
#   LEVER 1 — SIZE via large-seed heteroaggregation (UP1)
#     Cationic 35-42 µm seed particles + anionic colloids -> DLVO-driven
#     selective attachment -> D_agg ≈ seed + shell ≈ 45-50 µm.
#     In the model: stream$D_agg_um override after UP1, size_template = 1
#     to make UP4 track the aggregate rather than the atomizer droplet.
#
#   LEVER 2 — SURFACE FUSION suppression via DLVO + plasticizer stripping (UP1)
#     Lower ionic_strength + higher Delta_pH -> higher electrostatic barrier
#     -> higher E_stability -> lower theta_skin_pe contribution.
#     Remove plasticizer (C_plasticizer -> 0) -> Tg_plas closer to Tg_polymer
#     -> lower theta_skin_fus.
#     NOTE: hard floor at ~0.46 from Tg_polymer (330 K) < T_wet_bulb (373 K).
#     Large rigid seed physically prevents surface polymer from coalescing
#     even when polymer Tg is below the drying surface temperature.
#
#   LEVER 3 — POROSITY via dual gas+template strategy (UP1 + UP3)
#     a) Gas holdup: Q_gas up -> alpha_g up, then UP3 at low Fg
#        (alpha_floor = alpha_gas_cake * Fg_ref/Fg) preserves gas into dryer.
#     b) Capillary-bridge template (type 4, butyl butyrate):
#        phi_templ_free drives phi_porosity_z via the escape/inflation closure.
#     c) Structural porosity from large d_ratio inversion:
#        d_ratio = D_primary/D_agg; for D_agg=42 µm, d_ratio=0.82/42=0.019
#        -> phi_struct = 0.30*theta_skin*(1-d_ratio^(3-Df)) ~ 0.24.
#
# Scenario outputs -> unified_output/dvlo_50um_*.csv
# =============================================================================

library(deSolve)

source("unified/up1_mixer_module.R")
source("unified/interface_stream.R")
source("unified/up2_atomizer_dryer_module.R")
source("up3_separator_module.R")

out_dir <- "unified_output"
dir.create(out_dir, showWarnings = FALSE)

# =============================================================================
# NOMINAL FACTOR BASELINE  (capillary_bridge template type 4)
# =============================================================================
TEMPLATE_TYPE <- 4

nominal_x <- list(
  # UP1 mixer side
  v_tip           = 8.5,        # m/s  (just above aggregation onset)
  tau_mix         = 20.0,       # min
  P_mix           = 1.5,        # atm
  T_mix           = 325.0,      # K
  C_gas_diss_temp = 0.05,       # -  Henry-norm dissolved gas
  Q_gas           = 1.0,        # sparge flow
  Q_colloid       = 1.5,        # colloid feed flow
  Q_template      = 0.5,        # template feed flow
  C_solid_mass    = 0.25,       # wt/wt  nominal 25% solids
  C_temp_mass     = 0.04,       # wt/wt  rigid template mass fraction (not used in type 4)
  D_particle      = 1.25,       # µm  primary colloid input
  D_template      = 0.85,       # µm  template droplet
  C_binder        = 0.08,       # wt/wt
  C_monomer       = 0.015,      # wt/wt  residual monomer
  C_plasticizer   = 0.022,      # wt/wt
  C_surfactant    = 3e-3,       # wt/wt
  CMC             = 0.008,      # wt/wt
  HLB             = 12.0,       # -
  MW_surfactant   = 1414,       # g/mol
  A_molecule      = 6.35,       # nm2
  D_surf          = 5e-11,      # m2/s
  ionic_strength  = 0.25,       # M
  Delta_pH        = 2.5,        # pH above isoelectric point
  template_dose   = 0.5,        # template / interstitial void
  chi_template    = 1.4,        # Flory-Huggins chi (BB = poor solvent, chi>0.5)
  # UP2 dryer side (used for the UP4 stage — intentionally kept at nominal
  # because UP4 is the KNOWN path, not the study variable here)
  ALR             = 1.2,        # -
  P_atom_air      = 1.5e5,      # Pa
  P_feed          = 2.5e5,      # Pa
  mdot_L          = 0.016,      # kg/s
  sigma           = 0.05,       # N/m
  mu_L            = 0.005,      # Pa.s
  T_dryer_in      = 425.0,      # K
  mdot_gas_dry    = 0.50,       # kg/s
  Y_in            = 0.005,      # kg/kg
  t_hold          = 60,         # s
  D_b             = 8e-5,       # m
  Tg_polymer      = 330.0,      # K  ~ 57°C nominal
  n_flow          = 0.44,       # power-law index
  k_perm_mono     = 20.0,       # -
  k_perm_plast    = 15.0,       # -
  k_perm_bind     = 10.0,       # -
  T_bp_solv       = 439.0       # K  butyl butyrate bp = 166°C = 439 K
)

equipment_base <- up1_default_equipment()
equipment_base$template_type <- TEMPLATE_TYPE

up1_rename <- c(tau_mix = "tau", P_mix = "P_system", T_mix = "T_system")
up1_names  <- c("v_tip","tau_mix","P_mix","T_mix","C_gas_diss_temp","Q_gas",
                "Q_colloid","Q_template","C_solid_mass","C_temp_mass",
                "D_particle","D_template","C_binder","C_monomer","C_plasticizer",
                "C_surfactant","CMC","HLB","MW_surfactant","A_molecule",
                "D_surf","ionic_strength","Delta_pH","template_dose","chi_template")
up2_names  <- c("ALR","P_atom_air","P_feed","mdot_L","sigma","mu_L",
                "T_dryer_in","mdot_gas_dry","Y_in","t_hold","D_b",
                "Tg_polymer","n_flow","k_perm_mono","k_perm_plast",
                "k_perm_bind","T_bp_solv")

# Rename UP1 params to match ODE internal names
up1_pars_from_x <- function(x) {
  p <- as.list(x[up1_names])
  for (un in names(up1_rename)) {
    p[[up1_rename[[un]]]] <- p[[un]]
    p[[un]] <- NULL
  }
  p
}

# =============================================================================
# CHAIN RUNNER (base nominal chain, UP4 fixed at nominal)
# =============================================================================
run_chain <- function(x, ttype = TEMPLATE_TYPE,
                      up3_pars = list(),
                      verbose = FALSE) {
  eq <- equipment_base; eq$template_type <- ttype
  p1 <- up1_pars_from_x(x)
  r1 <- tryCatch(up1_run_mixer(p1, eq), error = function(e) NULL)
  if (is.null(r1) || any(is.na(r1$outputs))) return(NULL)

  stream <- stream_from_up1(r1, p1, eq)
  if (verbose) print_stream(stream, "mixer exit")

  # UP2 foam-wash column (placeholder — wired if available)
  options(unified.wire_up2 = FALSE, unified.wire_up3 = FALSE)
  stream <- intermediate_stage_1(stream)

  # UP3 decanting separator (always wired here)
  stream <- up3_separator(stream, up3_pars)
  if (verbose) print_stream(stream, "UP3 exit -> dryer feed")

  x_up2 <- as.list(x[up2_names])
  r2 <- tryCatch(up2_run_dryer(stream, x_up2), error = function(e) NULL)
  if (is.null(r2)) return(NULL)

  list(up1 = r1$outputs, stream = stream, up2 = r2,
       D_agg_um = stream$D_agg_um, alpha_g = stream$alpha_g)
}

# =============================================================================
# LARGE-SEED CHAIN RUNNER
# Override D_agg_um after UP1 to simulate heteroaggregation onto a large seed
# particle (e.g. 35-42 µm cationic PMMA/silica). Then enable size_template=1
# so UP4 tracks the aggregate size instead of the atomizer droplet.
#
# D_seed_um   : seed particle diameter (e.g. 40 µm)
# SHELL_UM    : thickness of the colloid coating shell per side (µm)
# D_primary_um: physical colloid size to use for d_ratio (default 0.82 µm)
# =============================================================================
run_chain_large_seed <- function(x, D_seed_um, SHELL_UM = 5.0,
                                  D_primary_um = 0.82,
                                  ttype = TEMPLATE_TYPE,
                                  up3_pars = list(),
                                  verbose = FALSE) {
  eq <- equipment_base; eq$template_type <- ttype
  p1 <- up1_pars_from_x(x)
  r1 <- tryCatch(up1_run_mixer(p1, eq), error = function(e) NULL)
  if (is.null(r1) || any(is.na(r1$outputs))) return(NULL)

  stream <- stream_from_up1(r1, p1, eq)

  # === HETEROAGGREGATION OVERRIDE ===
  # D_agg_um is now set by the seed + coating shell geometry, not the
  # calibrated 10.2 µm homogeneous aggregation plateau.
  D_agg_seed <- D_seed_um + 2 * SHELL_UM       # seed + two shell radii
  stream$D_agg_um           <- D_agg_seed
  # Set D_primary_phys_um to a value consistent with the DLVO heteroaggregation
  # mechanism: the colloid primaries coating the seed are 0.82 µm calibration
  # beads. This gives d_ratio = 0.82/D_agg_seed which drives phi_struct.
  stream$D_primary_phys_um  <- D_primary_um
  stream$D_primary_exit_um  <- D_primary_um
  if (verbose) cat(sprintf("  [seed heteroagg] D_agg=%.1f µm, d_ratio=%.4f\n",
                           D_agg_seed, D_primary_um / D_agg_seed))

  options(unified.wire_up2 = FALSE, unified.wire_up3 = FALSE)
  stream <- intermediate_stage_1(stream)
  stream <- up3_separator(stream, up3_pars)
  if (verbose) print_stream(stream, "UP3 exit (large-seed)")

  x_up2 <- as.list(x[up2_names])
  x_up2[["size_template"]] <- 1    # particle size tracks aggregate, not droplet

  r2 <- tryCatch(up2_run_dryer(stream, x_up2), error = function(e) NULL)
  if (is.null(r2)) return(NULL)

  list(up1 = r1$outputs, stream = stream, up2 = r2,
       D_agg_um = D_agg_seed, alpha_g = stream$alpha_g,
       D_primary_um = D_primary_um)
}

# Helper: extract key product scalars from a run result
extract_kpis <- function(res) {
  if (is.null(res)) return(NULL)
  list(
    D_agg_um       = res$D_agg_um,
    alpha_g_feed   = res$alpha_g,
    Dp50_um        = unname(res$up2["D_particle_um"]),   # already µm
    Dp90_um        = unname(res$up2["Dp90_um"]),          # already µm
    phi_porosity   = unname(res$up2["phi_porosity_z"]),
    theta_skin     = unname(res$up2["theta_skin_z"]),
    rho_tapped     = unname(res$up2["rho_tapped"]),
    Tg_product_K   = unname(res$up2["Tg_product_K"]),
    X_moisture     = unname(res$up2["X_moisture"]),
    solv_retained  = unname(res$up2["solv_retained"]),
    sigma_y_MPa    = unname(res$up2["sigma_y_cake_MPa"]),
    Blended_Porosity = unname(res$up1["Blended_Porosity"])
  )
}

cat("=== DLVO + Surface Chemistry Path to 50 µm High-Porosity Particles ===\n\n")

# =============================================================================
# SWEEP 0 — NOMINAL BASELINE
# =============================================================================
cat("--- Sweep 0: Nominal baseline ---\n")
nom <- run_chain(nominal_x)
nom_kpi <- extract_kpis(nom)
cat(sprintf("  D_agg=%.1f µm, Dp50=%.1f µm, phi=%.3f, theta_skin=%.3f, alpha_g=%.4f\n",
            nom_kpi$D_agg_um, nom_kpi$Dp50_um, nom_kpi$phi_porosity,
            nom_kpi$theta_skin, nom_kpi$alpha_g_feed))

# =============================================================================
# SWEEP 1 — DLVO CHEMISTRY: Delta_pH × ionic_strength -> theta_skin + phi
# =============================================================================
cat("\n--- Sweep 1: DLVO chemistry (Delta_pH × ionic_strength) ---\n")

dpH_vals <- c(0.5, 1.0, 1.5, 2.5, 3.5, 5.0)
ion_vals  <- c(0.01, 0.05, 0.10, 0.25, 0.50)

sw1_rows <- list()
for (dpH in dpH_vals) {
  for (ion in ion_vals) {
    x <- nominal_x
    x$Delta_pH       <- dpH
    x$ionic_strength <- ion
    res <- run_chain(x)
    kpi <- extract_kpis(res)
    if (!is.null(kpi)) {
      sw1_rows[[length(sw1_rows) + 1]] <- data.frame(
        Delta_pH = dpH, ionic_strength = ion,
        theta_skin = kpi$theta_skin, phi_porosity = kpi$phi_porosity,
        D_agg_um = kpi$D_agg_um, alpha_g = kpi$alpha_g_feed,
        rho_tapped = kpi$rho_tapped
      )
    }
  }
}
sw1 <- do.call(rbind, sw1_rows)
sw1_file <- file.path(out_dir, "dvlo_50um_sweep1_dlvo_chemistry.csv")
write.csv(sw1, sw1_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw1), sw1_file))
cat(sprintf("  theta_skin range: [%.3f, %.3f]\n", min(sw1$theta_skin), max(sw1$theta_skin)))
cat(sprintf("  phi_porosity range: [%.3f, %.3f]\n", min(sw1$phi_porosity), max(sw1$phi_porosity)))
best1 <- sw1[which.min(sw1$theta_skin), ]
cat(sprintf("  Best theta_skin=%.3f at Delta_pH=%.1f, ionic_strength=%.2f\n",
            best1$theta_skin, best1$Delta_pH, best1$ionic_strength))

# =============================================================================
# SWEEP 2 — GAS + TEMPLATE POROSITY: Q_gas × template_dose -> alpha_g + phi
# =============================================================================
cat("\n--- Sweep 2: Gas holdup + template dose ---\n")

qgas_vals  <- c(0.4, 0.8, 1.2, 1.6, 2.0)
tdose_vals <- c(0.1, 0.3, 0.5, 0.8, 1.1)

sw2_rows <- list()
for (qg in qgas_vals) {
  for (td in tdose_vals) {
    x <- nominal_x
    x$Q_gas          <- qg
    x$template_dose  <- td
    res <- run_chain(x)
    kpi <- extract_kpis(res)
    if (!is.null(kpi)) {
      sw2_rows[[length(sw2_rows) + 1]] <- data.frame(
        Q_gas = qg, template_dose = td,
        alpha_g_feed = kpi$alpha_g_feed,
        phi_porosity = kpi$phi_porosity,
        theta_skin = kpi$theta_skin,
        Blended_Porosity = kpi$Blended_Porosity
      )
    }
  }
}
sw2 <- do.call(rbind, sw2_rows)
sw2_file <- file.path(out_dir, "dvlo_50um_sweep2_gas_template.csv")
write.csv(sw2, sw2_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw2), sw2_file))
cat(sprintf("  phi_porosity range: [%.3f, %.3f]\n", min(sw2$phi_porosity), max(sw2$phi_porosity)))
cat(sprintf("  alpha_g range: [%.4f, %.4f]\n", min(sw2$alpha_g_feed), max(sw2$alpha_g_feed)))

# =============================================================================
# SWEEP 3 — UP3 G-FORCE: Fg -> gas retention -> porosity in dryer feed
# =============================================================================
cat("\n--- Sweep 3: UP3 centrifugal force (Fg) ---\n")

Fg_vals <- c(50, 100, 150, 200, 300, 430)

sw3_rows <- list()
for (Fg in Fg_vals) {
  up3_p <- list(Fg = Fg)
  # Use elevated Q_gas to have gas available to preserve
  x <- nominal_x; x$Q_gas <- 1.6
  res <- run_chain(x, up3_pars = up3_p)
  kpi <- extract_kpis(res)
  if (!is.null(res) && !is.null(kpi)) {
    # UP3 module physics: alpha_floor = alpha_gas_cake * Fg_ref/Fg
    alpha_floor <- 0.076 * (430 / Fg)
    sw3_rows[[length(sw3_rows) + 1]] <- data.frame(
      Fg_g           = Fg,
      alpha_floor_th = min(alpha_floor, 0.99),   # theoretical maximum gas floor
      alpha_g_feed   = kpi$alpha_g_feed,
      phi_porosity   = kpi$phi_porosity,
      theta_skin     = kpi$theta_skin,
      Cs_cake        = res$stream$C_solid
    )
  }
}
sw3 <- do.call(rbind, sw3_rows)
sw3_file <- file.path(out_dir, "dvlo_50um_sweep3_up3_gforce.csv")
write.csv(sw3, sw3_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw3), sw3_file))
cat("  Fg [g] | alpha_floor | alpha_g_feed | phi_porosity\n")
for (i in seq_len(nrow(sw3))) {
  cat(sprintf("  %5.0f  | %10.4f  | %12.4f | %.4f\n",
              sw3$Fg_g[i], sw3$alpha_floor_th[i],
              sw3$alpha_g_feed[i], sw3$phi_porosity[i]))
}

# =============================================================================
# SWEEP 4 — LARGE SEED HETEROAGGREGATION: D_seed -> Dp50 + phi_struct
# =============================================================================
cat("\n--- Sweep 4: Large-seed heteroaggregation (D_seed sweep) ---\n")

D_seed_vals <- c(10, 15, 20, 25, 30, 35, 40, 45)
shell_um    <- 5.0     # colloid coating shell thickness per side [µm]
D_prim_um   <- 0.82    # physical colloid bead size for d_ratio [µm]

sw4_rows <- list()
for (Ds in D_seed_vals) {
  res <- run_chain_large_seed(nominal_x, D_seed_um = Ds,
                              SHELL_UM = shell_um, D_primary_um = D_prim_um)
  kpi <- extract_kpis(res)
  if (!is.null(kpi)) {
    D_agg  <- Ds + 2 * shell_um
    d_ratio <- D_prim_um / D_agg
    # Analytical phi_struct estimate (Df~2.5 for moderate stability)
    Df_est  <- 2.5
    phi_str_est <- 0.30 * kpi$theta_skin * (1 - d_ratio^(3 - Df_est))
    sw4_rows[[length(sw4_rows) + 1]] <- data.frame(
      D_seed_um    = Ds,
      D_agg_um     = D_agg,
      d_ratio      = d_ratio,
      phi_struct_est = phi_str_est,
      Dp50_um      = kpi$Dp50_um,
      Dp90_um      = kpi$Dp90_um,
      phi_porosity = kpi$phi_porosity,
      theta_skin   = kpi$theta_skin,
      rho_tapped   = kpi$rho_tapped
    )
  }
}
sw4 <- do.call(rbind, sw4_rows)
sw4_file <- file.path(out_dir, "dvlo_50um_sweep4_large_seed.csv")
write.csv(sw4, sw4_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw4), sw4_file))
cat("  D_seed | D_agg | Dp50 | phi_porosity | theta_skin | rho_tapped\n")
for (i in seq_len(nrow(sw4))) {
  cat(sprintf("  %6.0f | %5.0f | %4.1f | %12.4f | %10.4f | %.1f\n",
              sw4$D_seed_um[i], sw4$D_agg_um[i], sw4$Dp50_um[i],
              sw4$phi_porosity[i], sw4$theta_skin[i], sw4$rho_tapped[i]))
}

# =============================================================================
# SWEEP 5 — PLASTICIZER STRIPPING: C_plasticizer -> theta_skin (DLVO + Tg effect)
# =============================================================================
cat("\n--- Sweep 5: Plasticizer stripping (C_plasticizer) ---\n")

cplas_vals <- c(0.0, 0.005, 0.010, 0.015, 0.020, 0.030, 0.040)
sw5_rows   <- list()
for (cp in cplas_vals) {
  x <- nominal_x; x$C_plasticizer <- cp
  res <- run_chain(x)
  kpi <- extract_kpis(res)
  if (!is.null(kpi)) {
    sw5_rows[[length(sw5_rows) + 1]] <- data.frame(
      C_plasticizer = cp,
      theta_skin    = kpi$theta_skin,
      phi_porosity  = kpi$phi_porosity,
      Tg_product_K  = kpi$Tg_product_K,
      Dp50_um       = kpi$Dp50_um
    )
  }
}
sw5 <- do.call(rbind, sw5_rows)
sw5_file <- file.path(out_dir, "dvlo_50um_sweep5_plasticizer.csv")
write.csv(sw5, sw5_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw5), sw5_file))
cat("  C_plas | theta_skin | phi_porosity | Tg_product_K\n")
for (i in seq_len(nrow(sw5))) {
  cat(sprintf("  %.4f | %10.4f | %12.4f | %.1f\n",
              sw5$C_plasticizer[i], sw5$theta_skin[i],
              sw5$phi_porosity[i], sw5$Tg_product_K[i]))
}

# =============================================================================
# SWEEP 6 — Tg_polymer sweep: effect on surface fusion floor
# =============================================================================
cat("\n--- Sweep 6: Tg_polymer sweep ---\n")

tg_vals  <- c(280, 300, 320, 340, 360, 380)
sw6_rows <- list()
for (tg in tg_vals) {
  x <- nominal_x; x$Tg_polymer <- tg; x$C_plasticizer <- 0.005
  res <- run_chain(x)
  kpi <- extract_kpis(res)
  if (!is.null(kpi)) {
    # Estimated Tg_plas using Fox equation (without UP4's w_core)
    phi_solvent_est <- x$C_monomer + x$C_plasticizer
    inv_tg <- (1 - phi_solvent_est) / tg + phi_solvent_est / 150
    Tg_plas_est <- 1 / inv_tg
    theta_fus_est <- 1 / (1 + exp(-(373.15 - Tg_plas_est) / 15))
    sw6_rows[[length(sw6_rows) + 1]] <- data.frame(
      Tg_polymer_K   = tg,
      Tg_plas_est_K  = Tg_plas_est,
      theta_fus_est  = theta_fus_est,
      theta_skin     = kpi$theta_skin,
      phi_porosity   = kpi$phi_porosity
    )
  }
}
sw6 <- do.call(rbind, sw6_rows)
sw6_file <- file.path(out_dir, "dvlo_50um_sweep6_tg_polymer.csv")
write.csv(sw6, sw6_file, row.names = FALSE)
cat(sprintf("  Wrote %d rows -> %s\n", nrow(sw6), sw6_file))
cat("  Tg_poly | Tg_plas | theta_fus | theta_skin_actual\n")
for (i in seq_len(nrow(sw6))) {
  cat(sprintf("  %7.0f | %7.1f | %9.4f | %.4f\n",
              sw6$Tg_polymer_K[i], sw6$Tg_plas_est_K[i],
              sw6$theta_fus_est[i], sw6$theta_skin[i]))
}

# =============================================================================
# COMBINED BEST-CASE: 50 µm, high porosity, low surface fusion
# Apply all three levers simultaneously:
#   1. Large seed: D_seed = 40 µm -> D_agg = 50 µm, size_template = 1
#   2. DLVO: Delta_pH = 3.5, ionic_strength = 0.01 -> low coagulation rate
#   3. Gas + template: Q_gas = 1.6, template_dose = 0.8
#   4. Plasticizer stripped: C_plasticizer = 0.005
#   5. High Tg polymer: Tg_polymer = 360 K (87°C)
#   6. UP3 at low Fg: Fg = 100 g -> alpha_floor ~ 33% gas preserved
# =============================================================================
cat("\n=== COMBINED BEST-CASE RECIPE ===\n")

x_best <- nominal_x
x_best$Delta_pH        <- 3.5
x_best$ionic_strength  <- 0.03
x_best$Q_gas           <- 1.6
x_best$template_dose   <- 0.8
x_best$chi_template    <- 1.4    # butyl butyrate: RED=1.33, chi~1.4 (poor solvent)
x_best$T_bp_solv       <- 439    # butyl butyrate bp = 166°C
x_best$C_plasticizer   <- 0.005  # minimal plasticizer
x_best$C_monomer       <- 0.005  # reduce residual monomer
x_best$Tg_polymer      <- 355    # high-Tg variant polymer (82°C)

up3_best <- list(Fg = 100, alpha_gas_cake = 0.076, Fg_ref = 430,
                 Cs_target = 0.405, reject_solid_frac = 0.003)

D_seed_best <- 40.0

cat("  Conditions:\n")
cat(sprintf("    Seed diameter: %.0f µm (D_agg_target = %.0f µm)\n",
            D_seed_best, D_seed_best + 2 * shell_um))
cat(sprintf("    Delta_pH=%.1f, ionic_strength=%.3f M\n",
            x_best$Delta_pH, x_best$ionic_strength))
cat(sprintf("    Q_gas=%.1f, template_dose=%.1f\n",
            x_best$Q_gas, x_best$template_dose))
cat(sprintf("    C_plasticizer=%.4f, Tg_polymer=%.0f K\n",
            x_best$C_plasticizer, x_best$Tg_polymer))
cat(sprintf("    UP3 Fg=%.0f g (alpha_floor_th=%.3f)\n",
            up3_best$Fg,
            up3_best$alpha_gas_cake * up3_best$Fg_ref / up3_best$Fg))

res_best <- run_chain_large_seed(x_best,
                                 D_seed_um     = D_seed_best,
                                 SHELL_UM      = shell_um,
                                 D_primary_um  = D_prim_um,
                                 up3_pars      = up3_best,
                                 verbose       = FALSE)

kpi_best <- extract_kpis(res_best)

if (!is.null(kpi_best)) {
  cat("\n  RESULTS:\n")
  cat(sprintf("    D_agg_um (stream)    : %.1f µm\n",  kpi_best$D_agg_um))
  cat(sprintf("    Dp50_um (product)    : %.1f µm\n",  kpi_best$Dp50_um))
  cat(sprintf("    Dp90_um              : %.1f µm\n",  kpi_best$Dp90_um))
  cat(sprintf("    phi_porosity         : %.4f (%.1f%%)\n",
              kpi_best$phi_porosity, 100 * kpi_best$phi_porosity))
  cat(sprintf("    theta_skin           : %.4f\n",     kpi_best$theta_skin))
  cat(sprintf("    rho_tapped kg/m3     : %.1f\n",     kpi_best$rho_tapped))
  cat(sprintf("    Tg_product_K         : %.1f K (%.0f°C)\n",
              kpi_best$Tg_product_K, kpi_best$Tg_product_K - 273.15))
  cat(sprintf("    X_moisture           : %.5f\n",     kpi_best$X_moisture))
  cat(sprintf("    solv_retained        : %.5f\n",     kpi_best$solv_retained))
  cat(sprintf("    sigma_y_cake_MPa     : %.4f\n",     kpi_best$sigma_y_MPa))
  cat(sprintf("    alpha_g_feed (UP3)   : %.4f\n",     kpi_best$alpha_g_feed))

  # Save best-case stream for inspection
  best_stream_file <- file.path(out_dir, "dvlo_50um_best_case_stream.csv")
  s <- res_best$stream
  stream_df <- data.frame(
    field = names(s)[vapply(s, function(v) is.numeric(v) && length(v) == 1, logical(1))],
    stringsAsFactors = FALSE
  )
  stream_df$value <- unlist(s[stream_df$field])
  write.csv(stream_df, best_stream_file, row.names = FALSE)
  cat(sprintf("\n  Stream state saved -> %s\n", best_stream_file))

  # Save best-case product properties
  kpi_df <- data.frame(field = names(kpi_best), value = unlist(kpi_best),
                       stringsAsFactors = FALSE)
  kpi_file <- file.path(out_dir, "dvlo_50um_best_case_kpis.csv")
  write.csv(kpi_df, kpi_file, row.names = FALSE)
  cat(sprintf("  KPIs saved           -> %s\n", kpi_file))
} else {
  cat("  WARNING: best-case run returned NULL (model error)\n")
}

# =============================================================================
# COMPARISON TABLE: nominal vs best-case
# =============================================================================
cat("\n=== COMPARISON: Nominal vs 50 µm Best-Case ===\n")
cat(sprintf("  %-24s %12s %12s\n", "KPI", "Nominal", "Best-Case (50µm)"))
cat(paste0(rep("-", 52), collapse = ""), "\n")

kpi_fields <- c("D_agg_um", "Dp50_um", "Dp90_um", "phi_porosity",
                "theta_skin", "rho_tapped", "Tg_product_K", "alpha_g_feed")
for (f in kpi_fields) {
  nom_v  <- if (!is.null(nom_kpi[[f]])) nom_kpi[[f]] else NA
  best_v <- if (!is.null(kpi_best[[f]])) kpi_best[[f]] else NA
  cat(sprintf("  %-24s %12.4g %12.4g\n", f, nom_v, best_v))
}

cat("\n=== DONE ===\n")
cat("Output files in", out_dir, ":\n")
for (f in list.files(out_dir, pattern = "dvlo_50um.*\\.csv")) {
  cat(sprintf("  %s\n", f))
}
