# =============================================================================
# UP1 MODULE: gassed / templated mixing vessel (Latex Coagulation Engine)
# =============================================================================
# Physics extracted verbatim from up1_module_rev38_dryer_risk.r (v52 / rev38):
# two-compartment (cavern A / dead-wall B) ODE model with C-infinity smoothed
# algebra, surface chemistry (CMC cap, HLB efficiency, DLVO), templating
# closure (rigid / gas / surface-weld / capillary-bridge) and the rev38
# spray-dryer risk outputs (Swelling_Softness_exit, Residual_Template_Fraction).
#
# This file contains NO study driver (no Morris, no sweeps, no plotting) —
# only the physics, packaged as:
#
#   up1_run_mixer(pars, equipment) -> list(outputs = <named numeric>,
#                                          extras  = <named list>)
#
# `outputs` are the 13 rev38 screening outputs; `extras` carries the exit
# state the downstream stream interface needs (diluted solids, phi_s, slurry
# density, template fill, exit T/P, Theta_surf).
#
# Keep this file in sync with up1_module_rev38_dryer_risk.r if the upstream
# module evolves — the RHS and post-processing below are a 1:1 copy.
# =============================================================================

library(deSolve)

# Saturation regime boundaries shared with the templating closure
UP1_PEND_MAX <- 0.25
UP1_CAP_MAX  <- 0.80
UP1_EPS_RCP  <- 0.36   # random-close-packing void fraction of the packed solids

up1_default_equipment <- function() {
  list(chi_npgas    = 0.0,
       D_imp2D_tank = 0.71,
       N_impellers  = 4.0,
       template_type = 1)   # 1=rigid, 2=gas, 3=surface_weld, 4=capillary_bridge
}

# --- ODE right-hand side (verbatim from rev38) -------------------------------
up1_kinetics_rhs <- function(t, state, pars) {
  with(as.list(c(state, pars)), {

    # --- C-INFINITY SMOOTHING ALGEBRA ---
    s_pos <- function(x) 0.5 * (x + sqrt(x^2 + 1e-8))
    s_min <- function(a, b) 0.5 * ((a + b) - sqrt((a - b)^2 + 1e-8))
    s_max <- function(a, b) 0.5 * ((a + b) + sqrt((a - b)^2 + 1e-8))
    cbrt  <- function(x) sign(x) * abs(x)^(1/3)

    # --- 1. GLOBAL GEOMETRY & RHEOLOGY ---
    D_tank_ref <- 1.0
    D_impeller <- D_imp2D_tank * D_tank_ref
    V_tank_total <- (pi / 4) * (D_tank_ref^3)

    denom <- (C_solid_mass / rho_polymer) + ((1.0 - C_solid_mass) / rho_water)
    phi_s <- (C_solid_mass / rho_polymer) / denom

    # Cavitation & Compressible Expansion
    P_local <- s_pos(P_system - (0.001 * (v_tip^2))) + 0.05
    T_scale <- (T_system / 298.15)
    rho_gas_local <- 1.2 * (P_local / T_scale)

    rho_eff <- rho_water * (1.0 - a_free_A) + rho_gas_local * a_free_A
    safe_rho <- s_pos(rho_eff) + 1.0

    # Base Fluid & Swelling Plasticization
    phi_solvent <- C_monomer + C_plasticizer
    Swelling_Softness <- exp(25.0 * phi_solvent)

    K_fluid <- 1.5 + (10.0 * C_binder)
    n_fluid <- 0.65
    shear_rate <- v_tip / D_impeller

    # --- 2. ZONE A PHYSICS (THE CAVERN) ---
    eps_0D <- N_impellers * (0.005 * (v_tip^3) / D_impeller)
    epsilon_A_raw <- eps_0D * (V_tank_total / (s_pos((pi/6)*D_impeller^3) + 1e-4)) * s_pos(1.0 - 1.65 * a_free_A) + 1e-4
    epsilon_A <- 1000.0 * tanh(epsilon_A_raw / 1000.0)

    softness_A <- (1.0 + (5.0 * C_binder)) * Swelling_Softness

    # --- 1.5 SURFACE CHEMISTRY, DLVO ELECTROSTATICS & COVERAGE ---
    C_surf_eff <- s_min(C_surfactant, CMC)
    Moles_Surf <- C_surf_eff / s_pos(MW_surfactant + 1e-6)

    HLB_efficiency <- s_max(0.1, 1.0 - 0.01 * (HLB - 12.0)^2)
    Coverage_Capacity <- Moles_Surf * A_molecule * 1e6 * HLB_efficiency
    SA_Gas_Dynamic <- s_pos(a_free_A + a_trap_A) * (s_pos(epsilon_A)^0.5) + 1e-6

    Theta_surf_raw <- Coverage_Capacity / SA_Gas_Dynamic
    Theta_surf <- s_min(1.0, Theta_surf_raw)
    Foam_Stability <- 1.0 + (5.0 * Theta_surf * (HLB / 10.0))

    # DLVO Exponential Energy Barrier
    dpH <- s_pos(Delta_pH)
    W_barrier <- exp(1.0 * dpH - 5.0 * sqrt(ionic_strength))
    E_stability <- s_max(0.01, s_min(100.0, W_barrier)) * s_pos(1.0 + 10.0 * Theta_surf)

    Debye_Screening <- 1.0 / s_pos(1.0 + 10.0 * ionic_strength)
    E_repulsion <- dpH * Debye_Screening * s_pos(1.0 - Theta_surf)

    Mobility_Factor <- s_max(0.1, 1000.0 / s_pos(MW_surfactant))

    # --- 1.6 YIELD STRESS (wired to additive solvents) ---
    Yield_Stress <- ((1.0 + (K_fluid * (shear_rate^n_fluid))) / Swelling_Softness) * (1.0 + 100.0 * s_pos(a_trap_A) * phi_s * Foam_Stability)

    # Volume Mapping
    N_rps <- v_tip / (pi * D_impeller)
    Power_Number <- 0.5
    Dc_cubed_ratio <- (Power_Number * safe_rho * (N_rps^2) * (D_impeller^2)) / Yield_Stress
    D_cavern_single <- D_impeller * cbrt(Dc_cubed_ratio)
    D_cavern_actual <- s_min(D_tank_ref, s_max(D_impeller, D_cavern_single))
    V_cavern_single <- (pi / 6) * (D_cavern_actual^3)

    V_stage_max <- V_tank_total / s_pos(N_impellers)
    V_cavern_isolated <- s_min(V_stage_max, V_cavern_single)
    V_A <- s_min(V_tank_total, V_cavern_isolated * N_impellers)
    V_B <- s_max(0.01 * V_tank_total, V_tank_total - V_A)

    Q_ideal <- 0.7 * N_rps * (D_impeller^3) * N_impellers
    cavern_fraction <- V_A / V_tank_total
    Q_exchange <- Q_ideal * (cavern_fraction^3)
    ex_rate_A <- 10.0 * tanh((Q_exchange / V_A) / 10.0)
    ex_rate_B <- 10.0 * tanh((Q_exchange / V_B) / 10.0)

    # --- AGGREGATE BREAKAGE & FORMATION ---
    Rate_floc_A <- (0.01 * (0.2 / s_max(0.01, D_particle)) * (phi_s^2) * epsilon_A * softness_A) / E_stability

    A_template_A <- s_pos(n_temp_A) * pi * (D_template * 1e-6)^2
    A_latex_A <- (6.0 * phi_s) / s_pos(D_particle * 1e-6)
    Seed_Area_Frac <- A_template_A / s_pos(A_template_A + A_latex_A)
    S_adh_seed <- (1.0 / (1.0 + exp(1.0 * dpH))) * Debye_Screening * (1.0 + 5.0 * C_monomer)
    R_het_A <- 0.05 * Seed_Area_Frac * phi_s * epsilon_A * softness_A * S_adh_seed

    denom_break_A <- (1.0 + 0.5 * s_pos(a_trap_A)) * (1.0 + 0.5 * s_pos(t_skin_A)) * (1.0 + 100.0 * s_pos(t_surf_A)) * softness_A
    R_break_raw <- (0.001 * s_pos(D_particle - 0.05) * epsilon_A * s_pos(1.0 + E_repulsion)) / denom_break_A
    R_break_A <- 10.0 * tanh(R_break_raw / 10.0)

    Rate_N_A <- Rate_floc_A + R_het_A - (R_break_A * N_A)

    # --- GAS THERMODYNAMICS & FLASHING ---
    C_eq <- P_local / s_pos(5.0 * T_scale)
    S_flash_A <- 0.1 * s_pos(C_gas_A - C_eq)
    S_diss_A  <- 0.001 * s_pos(C_eq - C_gas_A) * epsilon_A
    Rate_C_gas_A <- S_diss_A - S_flash_A

    S_sparge_A <- k_sparge * (Q_gas / (s_pos(Q_colloid) + 0.01)) * (1.0 - a_free_A)
    capture_eff_A <- (1.0 * phi_s) / s_pos(1.0 + 0.1 * E_stability)
    J_free_trap_A <- (1.0 - chi_npgas) * capture_eff_A * a_free_A * (1.0 + (0.5 * Foam_Stability) * s_pos(a_trap_A)) * (1.0 - a_trap_A)

    Rate_a_free_A <- S_sparge_A + S_flash_A - J_free_trap_A
    Rate_a_trap_A <- J_free_trap_A

    # Morphological States
    k_omega <- 10.0 * tanh((0.01 * softness_A * epsilon_A) / 10.0)
    Rate_Omega_A <- k_omega * (1.0 - Omega_A) - 0.005 * Omega_A

    k_surf_base <- 10.0 * tanh((1.0 * sqrt(epsilon_A)) / 10.0)
    k_surf <- k_surf_base * Mobility_Factor * (1.0 + 1e9 * D_surf)
    Rate_t_surf_A <- k_surf * (1.0 - t_surf_A) - 0.01 * t_surf_A

    k_skin <- 10.0 * tanh((0.001 * epsilon_A * Swelling_Softness) / 10.0)
    Rate_t_skin_A <- (0.0001 * epsilon_A * Swelling_Softness) * (1.0 - t_skin_A) - k_skin * t_skin_A

    # --- 3. ZONE B PHYSICS (THE DEAD WALL) ---
    Rate_C_gas_B <- 0.0
    J_free_trap_B <- (1.0 - chi_npgas) * (0.1 * capture_eff_A) * a_free_B * (1.0 + (0.5 * Foam_Stability) * s_pos(a_trap_B)) * (1.0 - a_trap_B)

    Rate_Omega_B <- - 0.005 * Omega_B
    Rate_N_B <- 0.001 * phi_s^2
    Rate_a_free_B <- - J_free_trap_B
    Rate_a_trap_B <- J_free_trap_B
    Rate_t_surf_B <- 0.01 * Mobility_Factor * (1.0 - t_surf_B)
    Rate_t_skin_B <- 0.0

    # --- 4. COMPARTMENTAL MASS TRANSFER ---
    dN_A_dt <- Rate_N_A + ex_rate_A * (N_B - N_A)
    dN_B_dt <- Rate_N_B + ex_rate_B * (N_A - N_B)

    da_free_A_dt <- Rate_a_free_A + ex_rate_A * (a_free_B - a_free_A)
    da_free_B_dt <- Rate_a_free_B + ex_rate_B * (a_free_A - a_free_B)

    dC_gas_A_dt <- Rate_C_gas_A + ex_rate_A * (C_gas_B - C_gas_A)
    dC_gas_B_dt <- Rate_C_gas_B + ex_rate_B * (C_gas_A - C_gas_B)

    da_trap_A_dt <- Rate_a_trap_A + ex_rate_A * (a_trap_B - a_trap_A)
    da_trap_B_dt <- Rate_a_trap_B + ex_rate_B * (a_trap_A - a_trap_B)

    dt_skin_A_dt <- Rate_t_skin_A + ex_rate_A * (t_skin_B - t_skin_A)
    dt_skin_B_dt <- Rate_t_skin_B + ex_rate_B * (t_skin_A - t_skin_B)

    dOmega_A_dt <- Rate_Omega_A + ex_rate_A * (Omega_B - Omega_A)
    dOmega_B_dt <- Rate_Omega_B + ex_rate_B * (Omega_A - Omega_B)

    dt_surf_A_dt <- Rate_t_surf_A + ex_rate_A * (t_surf_B - t_surf_A)
    dt_surf_B_dt <- Rate_t_surf_B + ex_rate_B * (t_surf_A - t_surf_B)

    dC_chem_A_dt <- 0; dC_chem_B_dt <- 0
    dn_temp_A_dt <- ex_rate_A * (n_temp_B - n_temp_A)
    dn_temp_B_dt <- ex_rate_B * (n_temp_A - n_temp_B)

    return(list(c(dN_A_dt, da_free_A_dt, dC_gas_A_dt, dC_chem_A_dt, da_trap_A_dt, dt_skin_A_dt, dOmega_A_dt, dt_surf_A_dt, dn_temp_A_dt,
                  dN_B_dt, da_free_B_dt, dC_gas_B_dt, dC_chem_B_dt, da_trap_B_dt, dt_skin_B_dt, dOmega_B_dt, dt_surf_B_dt, dn_temp_B_dt)))
  })
}

up1_output_names <- c("Blended_Porosity", "Blended_Sphericity", "Mixing_Potential",
                      "Blended_WetSkin", "Blended_Size_um", "Blended_Viscosity_PaS",
                      "SG_Colloid", "Blended_SG_Gassed", "Bond_Strength",
                      "Retained_Porosity", "Template_Target_Score",
                      "Swelling_Softness_exit", "Residual_Template_Fraction")

# --- Single mixer run (rev38 post-processing, verbatim) ----------------------
# pars: named list/vector of the rev38 base_vars factors.
# equipment: list(chi_npgas, D_imp2D_tank, N_impellers, template_type).
up1_run_mixer <- function(pars, equipment = up1_default_equipment()) {
  p <- as.list(pars)
  for (name in names(equipment)) p[[name]] <- equipment[[name]]

  s_pos <- function(x) 0.5 * (x + sqrt(x^2 + 1e-8))
  s_clamp <- function(x, lower=0, upper=1) {
    x_low <- lower + 0.5 * ((x - lower) + sqrt((x - lower)^2 + 1e-8))
    upper - 0.5 * ((upper - x_low) + sqrt((upper - x_low)^2 + 1e-8))
  }
  s_min <- function(a, b) 0.5 * ((a + b) - sqrt((a - b)^2 + 1e-8))
  s_max <- function(a, b) 0.5 * ((a + b) + sqrt((a - b)^2 + 1e-8))
  cbrt <- function(x) sign(x) * abs(x)^(1/3)

  p$rho_water <- 1000.0; p$rho_polymer <- 1050.0; p$k_sparge <- 0.05

  t_type <- if (!is.null(p$template_type)) p$template_type else 1
  if (t_type == 2) p$chi_npgas <- 1.0

  Q_total <- p$Q_colloid + p$Q_template
  tau_effective <- p$tau * (p$Q_colloid / Q_total)

  C_solid_effective <- (p$Q_colloid * p$C_solid_mass) / Q_total

  denom_poly <- (C_solid_effective / p$rho_polymer) + ((1.0 - C_solid_effective) / p$rho_water)
  phi_s_effective <- (C_solid_effective / p$rho_polymer) / denom_poly

  denom_temp <- (p$C_temp_mass / 2000.0) + ((1.0 - p$C_temp_mass) / p$rho_water)
  phi_template_effective <- (((p$Q_template * p$C_temp_mass) / Q_total) / 2000.0) / denom_temp

  C_gas_diss_effective <- (p$Q_template * p$C_gas_diss_temp) / Q_total

  V_single_template <- (pi / 6) * (p$D_template * 1e-6)^3
  n_temp_initial <- phi_template_effective / V_single_template

  p$tau <- tau_effective
  p$C_solid_mass <- C_solid_effective

  state_init <- c(N_A=100.0, a_free_A=0.0, C_gas_A=C_gas_diss_effective, C_chem_A=0.1, a_trap_A=0.0, t_skin_A=0.0, Omega_A=0.2, t_surf_A=0.0, n_temp_A=n_temp_initial,
                  N_B=100.0, a_free_B=0.0, C_gas_B=C_gas_diss_effective, C_chem_B=0.1, a_trap_B=0.0, t_skin_B=0.0, Omega_B=0.2, t_surf_B=0.0, n_temp_B=n_temp_initial)

  ode_out <- tryCatch({
    deSolve::ode(y = state_init, times = c(0, p$tau), func = up1_kinetics_rhs, parms = p, method = "lsoda", atol = 1e-6, rtol = 1e-4, maxsteps = 100000)
  }, error = function(e) { return(NULL) })

  if (is.null(ode_out) || nrow(ode_out) < 2)
    return(list(outputs = setNames(rep(NA_real_, length(up1_output_names)), up1_output_names),
                extras = NULL))

  final_state <- as.list(ode_out[nrow(ode_out), ])

  final_state$a_trap_A <- s_clamp(final_state$a_trap_A, 0, 0.99)
  final_state$a_free_A <- s_clamp(final_state$a_free_A, 0, 0.99)
  final_state$a_trap_B <- s_clamp(final_state$a_trap_B, 0, 0.99)
  final_state$a_free_B <- s_clamp(final_state$a_free_B, 0, 0.99)
  final_state$Omega_A <- s_clamp(final_state$Omega_A, 0, 1.0)
  final_state$Omega_B <- s_clamp(final_state$Omega_B, 0, 1.0)
  final_state$t_skin_A <- s_clamp(final_state$t_skin_A, 0, 1.0)
  final_state$t_skin_B <- s_clamp(final_state$t_skin_B, 0, 1.0)
  final_state$N_A <- s_pos(final_state$N_A) + 1e-6
  final_state$N_B <- s_pos(final_state$N_B) + 1e-6

  denom <- (p$C_solid_mass / p$rho_polymer) + ((1.0 - p$C_solid_mass) / p$rho_water)
  phi_s <- (p$C_solid_mass / p$rho_polymer) / denom

  K_fluid <- 1.5 + (10.0 * p$C_binder)
  n_fluid <- 0.65
  D_impeller <- p$D_imp2D_tank * 1.0
  shear_rate <- p$v_tip / s_max(0.01, D_impeller)

  C_surf_eff <- s_min(p$C_surfactant, p$CMC)
  Moles_Surf <- C_surf_eff / s_pos(p$MW_surfactant + 1e-6)
  HLB_efficiency <- s_max(0.1, 1.0 - 0.01 * (p$HLB - 12.0)^2)
  Coverage_Capacity <- Moles_Surf * p$A_molecule * 1e6 * HLB_efficiency
  SA_Gas_Dynamic_Final <- s_pos(final_state$a_free_A + final_state$a_trap_A) * 10.0 + 1e-6
  Theta_surf <- s_min(1.0, Coverage_Capacity / SA_Gas_Dynamic_Final)
  Foam_Stability <- 1.0 + (5.0 * Theta_surf * (p$HLB / 10.0))

  phi_solvent <- p$C_monomer + p$C_plasticizer
  Swelling_Softness <- exp(25.0 * phi_solvent)

  Yield_Stress <- ((1.0 + (K_fluid * (shear_rate^n_fluid))) / Swelling_Softness) * (1.0 + 100.0 * s_pos(final_state$a_trap_A) * phi_s * Foam_Stability)

  N_rps <- p$v_tip / (pi * s_max(0.01, D_impeller))
  Dc_cubed_ratio <- (0.5 * 1000.0 * (N_rps^2) * (D_impeller^2)) / s_max(0.1, Yield_Stress)

  D_cavern_single <- D_impeller * cbrt(Dc_cubed_ratio)
  V_cavern_single_raw <- (pi / 6) * (D_cavern_single^3)

  D_cavern_actual <- s_min(1.0, s_max(D_impeller, D_cavern_single))
  V_cavern_single <- (pi / 6) * (D_cavern_actual^3)

  V_tank_total <- (pi / 4) * (1.0^3)
  V_stage_max <- V_tank_total / s_max(1.0, p$N_impellers)

  Mixing_Potential <- V_cavern_single_raw / V_stage_max

  V_cavern_isolated <- s_min(V_stage_max, V_cavern_single)
  V_A <- s_min(V_tank_total, V_cavern_isolated * p$N_impellers)
  V_B <- s_max(0.01 * V_tank_total, V_tank_total - V_A)

  Fraction_A <- V_A / V_tank_total
  Fraction_B <- V_B / V_tank_total

  Porosity_A <- s_min(1.0 - phi_s, final_state$a_trap_A)
  Porosity_B <- s_min(1.0 - phi_s, final_state$a_trap_B)

  Blended_Porosity <- (Porosity_A * Fraction_A) + (Porosity_B * Fraction_B)
  Blended_Sphericity <- (final_state$Omega_A * Fraction_A) + (final_state$Omega_B * Fraction_B)
  Blended_WetSkin <- (final_state$t_skin_A * Fraction_A) + (final_state$t_skin_B * Fraction_B)

  V_solid_A <- phi_s
  V_template_A <- final_state$n_temp_A * ((pi * (p$D_template * 1e-6)^3) / 6)

  V_env_A <- (V_solid_A + V_template_A) / s_max(1e-4, 1.0 - (0.05 * final_state$a_trap_A) - final_state$a_free_A)
  Size_A_um <- cbrt((6.0 * V_env_A) / (pi * final_state$N_A)) * 1e6

  V_solid_B <- phi_s
  V_template_B <- final_state$n_temp_B * ((pi * (p$D_template * 1e-6)^3) / 6)
  V_env_B <- (V_solid_B + V_template_B) / s_max(1e-4, 1.0 - (0.05 * final_state$a_trap_B) - final_state$a_free_B)
  Size_B_um <- cbrt((6.0 * V_env_B) / (pi * final_state$N_B)) * 1e6

  Blended_Size_um <- (Size_A_um * Fraction_A) + (Size_B_um * Fraction_B)

  phi_effective_exit <- phi_s + 0.1 * Blended_WetSkin + V_template_A + (p$chi_npgas * final_state$a_free_A)
  phi_m <- 0.64
  safe_phi_exit <- s_min(phi_effective_exit, phi_m - 0.01)
  mu_base_exit <- 0.001 * (1 - (safe_phi_exit / phi_m))^(-2.5 * phi_m)
  Blended_Viscosity_PaS <- mu_base_exit * (1 + (100.0 / 50.0))^(0.7 - 1)

  rho_colloid_kg_m3 <- (phi_s * p$rho_polymer) + ((1.0 - phi_s) * p$rho_water)
  SG_Colloid <- rho_colloid_kg_m3 / 997.0

  gas_expansion_factor <- (p$T_system / 298.15) / s_max(0.05, p$P_system)
  rho_gas_exit <- 1.2 * (1.0 / gas_expansion_factor)
  Blended_SG_Gassed <- (SG_Colloid * (1.0 - Blended_Porosity)) + ((rho_gas_exit / 997.0) * Blended_Porosity)

  # --- Templating closure (rev37) ---
  fill <- s_pos(p$template_dose)

  sat <- tanh(fill / (0.5 * UP1_PEND_MAX))

  gamma_proxy <- s_max(0.1, 1.0 - 0.8 * Theta_surf)
  size_factor <- s_min(2.0, 0.7 / s_max(0.1, p$D_particle))
  weld_soft <- s_min(3.0, Swelling_Softness * (p$T_system / 298.15) * (p$tau / 20.0))

  A_template_f <- s_pos(final_state$n_temp_A) * pi * (p$D_template * 1e-6)^2
  A_latex_f    <- (6.0 * phi_s) / s_pos(p$D_particle * 1e-6)
  Seed_Frac_f  <- A_template_f / s_pos(A_template_f + A_latex_f)

  bond_rigid <- s_clamp(0.15 * Seed_Frac_f + 0.05 * sat, 0, 0.4)
  bond_gas   <- 0.02 + 0.0 * fill
  bond_weld  <- s_clamp(sat * s_min(1.0, 0.6 * weld_soft * gamma_proxy), 0, 1)
  bond_cap   <- s_clamp(sat * s_min(1.0, 0.9 * gamma_proxy * size_factor), 0, 1)

  coalesce_pen <- s_clamp(s_pos(fill - UP1_PEND_MAX) * weld_soft, 0, 1)
  densify_pen  <- s_clamp(s_pos(fill - UP1_CAP_MAX) / (1.0 - UP1_CAP_MAX), 0, 1)

  ret_rigid <- s_clamp(UP1_EPS_RCP * (1.0 - 0.6 * fill), 0, UP1_EPS_RCP)
  ret_gas   <- s_clamp(Blended_Porosity, 0, 0.99)
  ret_weld  <- s_clamp(UP1_EPS_RCP * (1.0 - coalesce_pen), 0, UP1_EPS_RCP)
  ret_cap   <- s_clamp(UP1_EPS_RCP * (1.0 - 0.4 * densify_pen), 0, UP1_EPS_RCP)

  Bond_Strength <- switch(as.character(t_type),
                          "1" = bond_rigid, "2" = bond_gas,
                          "3" = bond_weld,  "4" = bond_cap, bond_rigid)
  Retained_Porosity <- switch(as.character(t_type),
                          "1" = ret_rigid, "2" = ret_gas,
                          "3" = ret_weld,  "4" = ret_cap, ret_rigid)
  Template_Target_Score <- s_clamp(Bond_Strength * (Retained_Porosity / UP1_EPS_RCP), 0, 1)

  # --- Spray-dryer collapse risk (rev38) ---
  r_particle_m  <- (p$D_particle * 1e-6) / 2.0
  D_polymer_ref <- 5e-17
  T_scale_38    <- (p$T_system / 298.15)^1.5
  D_polymer_eff <- D_polymer_ref * Swelling_Softness * T_scale_38
  tau_D_s       <- (r_particle_m^2) / D_polymer_eff
  f_diff        <- s_clamp(1.0 - exp(-3.0 * (p$tau * 60.0) / tau_D_s), 0.0, 1.0)

  Swelling_Softness_exit <- Swelling_Softness

  Residual_Template_Fraction <- switch(as.character(t_type),
    "1" = 0.0,
    "2" = 0.0,
    "3" = f_diff,
    "4" = s_clamp(0.05 * f_diff * s_clamp(fill, 0.0, 1.0), 0.0, 0.15),
    0.0)

  outputs <- c(Blended_Porosity = Blended_Porosity,
               Blended_Sphericity = Blended_Sphericity,
               Mixing_Potential = Mixing_Potential,
               Blended_WetSkin = Blended_WetSkin,
               Blended_Size_um = Blended_Size_um,
               Blended_Viscosity_PaS = Blended_Viscosity_PaS,
               SG_Colloid = SG_Colloid,
               Blended_SG_Gassed = Blended_SG_Gassed,
               Bond_Strength = Bond_Strength,
               Retained_Porosity = Retained_Porosity,
               Template_Target_Score = Template_Target_Score,
               Swelling_Softness_exit = Swelling_Softness_exit,
               Residual_Template_Fraction = Residual_Template_Fraction)

  # Exit-state extras consumed by the stream interface (not screened outputs)
  extras <- list(C_solid_exit   = p$C_solid_mass,   # diluted by template stream
                 phi_s_exit     = phi_s,
                 rho_slurry     = rho_colloid_kg_m3,
                 rho_polymer    = p$rho_polymer,
                 T_exit_K       = p$T_system,
                 P_exit_Pa      = p$P_system * 101325,
                 Theta_surf_exit = Theta_surf,
                 template_type  = t_type,
                 template_fill  = fill,             # liquid template / void volume
                 C_temp_rigid   = (p$Q_template * p$C_temp_mass) / (p$Q_colloid + p$Q_template))

  list(outputs = outputs, extras = extras)
}
