# ==============================================================================
# INTEGRATED SENSITIVITY ANALYSIS CONTROLLER & PLOTTER (v49 - Restored Factors)
# ==============================================================================
# rev35: restores factors specified in the Latex Coagulation Engine master
# documents (v28.0.0 Production Package / revG / v30.0.0 Armored Architecture)
# that were dropped or left inert in rev34 (v48):
#   1. CMC cap on interfacial surfactant:  C_eff = min(C_surfactant, CMC)
#   2. HLB packing-efficiency penalty:     max(0.1, 1 - 0.01*(HLB-12)^2)
#   3. DLVO repulsion in breakage kernel:  (1 + E_repulsion), with
#      E_repulsion = dpH * Debye_Screening * (1 - Theta_surf)
#   4. Wet-skin shielding on breakage:     (1 + 0.5*theta_skin)
#   5. Particle-size flocculation kernel:  (0.2 / max(0.01, D_particle))
#   6. Rigid template phase un-gated from chi_npgas (conserved solid inventory
#      per spec Section 4/9), with corrected volume-fraction conversion
#      (missing /rho_template) and a heterogeneous seed-deposition source
#      R_het = k_het * A_template * phi_s * eps_reg * S_adh_seed  (revG 4.2)
# ==============================================================================

# 1. LIBRARIES
library(deSolve)
if (requireNamespace("sensitivity", quietly = TRUE)) {
  library(sensitivity)
} else if (file.exists("morris_shim.R")) {
  # Offline fallback: bundled Morris OAT implementation extracted from the
  # GPL-2 'sensitivity' package (see morris_shim.R header for attribution).
  source("morris_shim.R")
} else {
  stop("Neither the 'sensitivity' package nor morris_shim.R is available.")
}
library(parallel)
library(pbapply)
library(ggplot2)
library(ggrepel)

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================
execution_mode    <- "Morris"  
target_feature    <- "Mixing_Potential" 
trajectories      <- 25        
n_samples         <- 10        
force_recalculate <- TRUE     
database_file     <- sprintf("Morris_Omnibus_rev35_r%d.rds", trajectories)

# ==============================================================================
# 3. VARIABLE DICTIONARY (Transition Zone Bounds)
# ==============================================================================
base_vars <- list(
  v_tip             = c(1.0, 32.0),      
  tau               = c(5.0, 40.0),     
  P_system          = c(5/14.7, 64.7/14.7), 
  T_system          = c(290, 360),       
  C_gas_diss_temp   = c(0.0, 1.0),       
  Q_gas             = c(0.4, 2.0),       
  Q_colloid         = c(40/60, 200/60),       
  Q_template        = c(0.1, 2.0),       
  C_solid_mass      = c(0.10, 0.60),         
  C_temp_mass       = c(0.01, 0.10),     
  D_particle        = c(0.5, 2.0),       
  D_template        = c(0.2, 1.5),       
  C_binder          = c(0.01, 0.15),     
  C_monomer         = c(0.001, 0.03),    
  C_plasticizer     = c(0.005, 0.04),    
  C_surfactant      = c(0.00005, 0.015), 
  CMC               = c(0.001, 0.015),   
  HLB               = c(8.0, 16.0),      
  MW_surfactant     = c(200, 10000),     
  A_molecule        = c(0.2, 12.5),      
  D_surf            = c(1e-11, 1e-9),    
  ionic_strength    = c(0.01, 0.50),      
  Delta_pH          = c(0.0, 5.0)        
)

static_equipment <- list(
  chi_npgas         = 0.0,   
  D_imp2D_tank      = 0.71,  
  N_impellers       = 4.0    
)

# ==============================================================================
# 4. CHUNKING CONFIGURATION 
# ==============================================================================
sweep_targets <- list(
  v_tip          = list(bounds = c(1.0, 32.0),  n_chunks = 3), 
  C_solid_mass   = list(bounds = c(0.10, 0.60), n_chunks = 3), 
  Delta_pH       = list(bounds = c(0.0, 5.0),   n_chunks = 3), 
  ionic_strength = list(bounds = c(0.01, 0.50), n_chunks = 2), 
  C_surfactant   = list(bounds = c(1e-5, 0.015), n_chunks = 3)  
)

# ==============================================================================
# 5. PHYSICS ENGINE (100% Continuous Jacobian + Chemistry)
# ==============================================================================
kinetics_system_rhs <- function(t, state, pars) {
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

    # --------------------------------------------------------------------------
    # --- 1.5 SURFACE CHEMISTRY, DLVO ELECTROSTATICS & COVERAGE ---
    # --------------------------------------------------------------------------
    # [RESTORED #1] Only surfactant below the CMC is interfacially active;
    # micellized excess does not contribute to coverage (spec revG 3.4).
    C_surf_eff <- s_min(C_surfactant, CMC)
    Moles_Surf <- C_surf_eff / s_pos(MW_surfactant + 1e-6)

    # [RESTORED #2] HLB packing-efficiency penalty, peak efficiency at HLB ~ 12
    HLB_efficiency <- s_max(0.1, 1.0 - 0.01 * (HLB - 12.0)^2)
    Coverage_Capacity <- Moles_Surf * A_molecule * 1e6 * HLB_efficiency
    SA_Gas_Dynamic <- s_pos(a_free_A + a_trap_A) * (s_pos(epsilon_A)^0.5) + 1e-6

    Theta_surf_raw <- Coverage_Capacity / SA_Gas_Dynamic
    Theta_surf <- s_min(1.0, Theta_surf_raw)
    Foam_Stability <- 1.0 + (5.0 * Theta_surf * (HLB / 10.0))

    # DLVO Exponential Energy Barrier (Restores Delta_pH impact)
    dpH <- s_pos(Delta_pH)
    W_barrier <- exp(1.0 * dpH - 5.0 * sqrt(ionic_strength))
    E_stability <- s_max(0.01, s_min(100.0, W_barrier)) * s_pos(1.0 + 10.0 * Theta_surf)

    # [RESTORED #3] DLVO double-layer repulsion (report Section 6, Eq. 14):
    # electrostatic repulsion of bare surface assists redispersion of flocs.
    Debye_Screening <- 1.0 / s_pos(1.0 + 10.0 * ionic_strength)
    E_repulsion <- dpH * Debye_Screening * s_pos(1.0 - Theta_surf)

    Mobility_Factor <- s_max(0.1, 1000.0 / s_pos(MW_surfactant))
    
    # --------------------------------------------------------------------------
    # --- 1.6 UPDATED YIELD STRESS (Wired to Additive Solvents) ---
    # --------------------------------------------------------------------------
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
    
    # --------------------------------------------------------------------------
    # --- AGGREGATE BREAKAGE & FORMATION (Now respects Chemistry!) ---
    # --------------------------------------------------------------------------
    # [RESTORED #5] Smaller primary particles collide/flocculate faster
    # (report Eq. 15: R_homo ~ 0.2 / max(0.01, D_particle)).
    Rate_floc_A <- (0.01 * (0.2 / s_max(0.01, D_particle)) * (phi_s^2) * epsilon_A * softness_A) / E_stability

    # [RESTORED #6] Heterogeneous seed deposition onto rigid templates
    # (revG 4.2): R_het = k_het * A_template * phi_s * eps * S_adh_seed.
    # Template area is expressed as the fraction of total collision area.
    A_template_A <- s_pos(n_temp_A) * pi * (D_template * 1e-6)^2
    A_latex_A <- (6.0 * phi_s) / s_pos(D_particle * 1e-6)
    Seed_Area_Frac <- A_template_A / s_pos(A_template_A + A_latex_A)
    S_adh_seed <- (1.0 / (1.0 + exp(1.0 * dpH))) * Debye_Screening * (1.0 + 5.0 * C_monomer)
    R_het_A <- 0.05 * Seed_Area_Frac * phi_s * epsilon_A * softness_A * S_adh_seed

    # [RESTORED #3, #4] Breakage kernel per report Eq. 16: adds wet-skin
    # shielding (1 + 0.5*theta_skin) and DLVO repulsion boost (1 + E_repulsion).
    denom_break_A <- (1.0 + 0.5 * s_pos(a_trap_A)) * (1.0 + 0.5 * s_pos(t_skin_A)) * (1.0 + 100.0 * s_pos(t_surf_A)) * softness_A
    R_break_raw <- (0.001 * s_pos(D_particle - 0.05) * epsilon_A * s_pos(1.0 + E_repulsion)) / denom_break_A
    R_break_A <- 10.0 * tanh(R_break_raw / 10.0)

    Rate_N_A <- Rate_floc_A + R_het_A - (R_break_A * N_A)
    
    # --------------------------------------------------------------------------
    # --- GAS THERMODYNAMICS & FLASHING (Restores P_system & T_system) ---
    # --------------------------------------------------------------------------
    # Calculate equilibrium dissolved gas using Henry's Law proxy
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
    Rate_C_gas_B <- 0.0 # Stagnant zone does not actively flash/dissolve strongly
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
    # [RESTORED #6] Rigid template number is a conserved transported state
    # (spec 4.1: R_loss,template = 0 by default; exchange transport only).
    dn_temp_A_dt <- ex_rate_A * (n_temp_B - n_temp_A)
    dn_temp_B_dt <- ex_rate_B * (n_temp_A - n_temp_B)
    
    return(list(c(dN_A_dt, da_free_A_dt, dC_gas_A_dt, dC_chem_A_dt, da_trap_A_dt, dt_skin_A_dt, dOmega_A_dt, dt_surf_A_dt, dn_temp_A_dt,
                  dN_B_dt, da_free_B_dt, dC_gas_B_dt, dC_chem_B_dt, da_trap_B_dt, dt_skin_B_dt, dOmega_B_dt, dt_surf_B_dt, dn_temp_B_dt)))
  })
}

# ==============================================================================
# 5. EXECUTION WRAPPERS & POST-PROCESSING
# ==============================================================================
run_single_simulation <- function(parameters) {
  p <- as.list(parameters)
  for (name in names(static_equipment)) { p[[name]] <- static_equipment[[name]] }
  
  s_pos <- function(x) 0.5 * (x + sqrt(x^2 + 1e-8))
  s_clamp <- function(x, lower=0, upper=1) {
    x_low <- lower + 0.5 * ((x - lower) + sqrt((x - lower)^2 + 1e-8))
    upper - 0.5 * ((upper - x_low) + sqrt((upper - x_low)^2 + 1e-8))
  }
  s_min <- function(a, b) 0.5 * ((a + b) - sqrt((a - b)^2 + 1e-8))
  s_max <- function(a, b) 0.5 * ((a + b) + sqrt((a - b)^2 + 1e-8))
  cbrt <- function(x) sign(x) * abs(x)^(1/3)
  
  p$rho_water <- 1000.0; p$rho_polymer <- 1050.0; p$k_sparge <- 0.05
  
  Q_total <- p$Q_colloid + p$Q_template
  tau_effective <- p$tau * (p$Q_colloid / Q_total) 
  
  C_solid_effective <- (p$Q_colloid * p$C_solid_mass) / Q_total
  
  denom_poly <- (C_solid_effective / p$rho_polymer) + ((1.0 - C_solid_effective) / p$rho_water)
  phi_s_effective <- (C_solid_effective / p$rho_polymer) / denom_poly
  
  # [RESTORED #6] Corrected mass->volume-fraction conversion (was missing the
  # /rho_template division, rho_template = 2000 kg/m^3).
  denom_temp <- (p$C_temp_mass / 2000.0) + ((1.0 - p$C_temp_mass) / p$rho_water)
  phi_template_effective <- (((p$Q_template * p$C_temp_mass) / Q_total) / 2000.0) / denom_temp

  C_gas_diss_effective <- (p$Q_template * p$C_gas_diss_temp) / Q_total

  V_single_template <- (pi / 6) * (p$D_template * 1e-6)^3
  # [RESTORED #6] Rigid templates are a conserved solid inventory, NOT gated
  # by the non-penetrating-gas switch (chi_npgas gates gas crossing only).
  n_temp_initial <- phi_template_effective / V_single_template
  
  p$tau <- tau_effective
  p$C_solid_mass <- C_solid_effective
  
  state_init <- c(N_A=100.0, a_free_A=0.0, C_gas_A=C_gas_diss_effective, C_chem_A=0.1, a_trap_A=0.0, t_skin_A=0.0, Omega_A=0.2, t_surf_A=0.0, n_temp_A=n_temp_initial,
                  N_B=100.0, a_free_B=0.0, C_gas_B=C_gas_diss_effective, C_chem_B=0.1, a_trap_B=0.0, t_skin_B=0.0, Omega_B=0.2, t_surf_B=0.0, n_temp_B=n_temp_initial)
  
  ode_out <- tryCatch({
    deSolve::ode(y = state_init, times = c(0, p$tau), func = kinetics_system_rhs, parms = p, method = "lsoda", atol = 1e-6, rtol = 1e-4, maxsteps = 100000)
  }, error = function(e) { return(NULL) })
  
  if (is.null(ode_out) || nrow(ode_out) < 2) return(rep(NA, 8))
  
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
  
  # Ensure Chemistry computes fully in Post Processing
  # [RESTORED #1, #2] Mirror the CMC cap and HLB efficiency used in the RHS
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
  
  return(c(Blended_Porosity = Blended_Porosity, 
           Blended_Sphericity = Blended_Sphericity, 
           Mixing_Potential = Mixing_Potential, 
           Blended_WetSkin = Blended_WetSkin,
           Blended_Size_um = Blended_Size_um,
           Blended_Viscosity_PaS = Blended_Viscosity_PaS,
           SG_Colloid = SG_Colloid,
           Blended_SG_Gassed = Blended_SG_Gassed))
}

run_batch <- function(param_matrix, v_names) {
  mat <- as.matrix(param_matrix); colnames(mat) <- v_names
  cores <- max(1, parallel::detectCores() - 1)
  input_list <- lapply(seq_len(nrow(mat)), function(i) mat[i, ])
  
  if (.Platform$OS.type == "unix") { raw_results <- parallel::mclapply(input_list, run_single_simulation, mc.cores = cores)
  } else {
    cl <- parallel::makeCluster(cores)
    parallel::clusterEvalQ(cl, { library(deSolve) })
    script_env <- environment(run_single_simulation)
    parallel::clusterExport(cl, c("run_single_simulation", "kinetics_system_rhs", "static_equipment"), envir = script_env)
    raw_results <- pbapply::pblapply(input_list, run_single_simulation, cl = cl)
    parallel::stopCluster(cl)
  }
  
  safe_results <- lapply(raw_results, function(x) { if (is.numeric(x) && length(x) == 8) return(x) else return(rep(NA, 8)) })
  res_df <- as.data.frame(do.call(rbind, safe_results), row.names = FALSE)
  colnames(res_df) <- c("Blended_Porosity", "Blended_Sphericity", "Mixing_Potential", "Blended_WetSkin", "Blended_Size_um", "Blended_Viscosity_PaS", "SG_Colloid", "Blended_SG_Gassed")
  return(res_df)
}

format_time <- function(s) {
  if (is.na(s) || is.infinite(s) || s < 0) return("Calculating...")
  sprintf("%02d:%02d:%02d", as.integer(s %/% 3600), as.integer((s %% 3600) %/% 60), as.integer(s %% 60))
}

# ==============================================================================
# 6. EXECUTION ENGINE: MORRIS OMNIBUS
# ==============================================================================
if (file.exists(database_file) && !force_recalculate) {
  cat(sprintf("\n>>> Master Database found! Loading pre-calculated Morris simulations from '%s'...\n", database_file))
  omnibus_data <- readRDS(database_file)
} else {
  cat("\n>>> Starting fresh MORRIS OMNIBUS execution...\n")
  cat(">>> Generating Auto-Chunked Scenarios...\n")
  
  chunk_intervals <- list()
  for (var in names(sweep_targets)) {
    if (var %in% names(base_vars)) {
      min_b <- sweep_targets[[var]]$bounds[1]; max_b <- sweep_targets[[var]]$bounds[2]
      n <- sweep_targets[[var]]$n_chunks; step <- (max_b - min_b) / n
      ints <- list()
      for (i in 1:n) ints[[i]] <- c(min_b + (i-1)*step, min_b + i*step)
      chunk_intervals[[var]] <- ints
    }
  }
  grid_indices <- expand.grid(lapply(chunk_intervals, seq_along))
  total_scenarios <- nrow(grid_indices)
  
  global_start_time <- Sys.time()
  omnibus_data <- list()
  
  for (row in 1:total_scenarios) {
    current_time <- Sys.time()
    if (row == 1) { eta_str <- "Calculating..." } else {
      elapsed <- as.numeric(difftime(current_time, global_start_time, units = "secs"))
      eta_str <- format_time((elapsed / (row - 1)) * (total_scenarios - row + 1))
    }
    
    progress_pct <- row / total_scenarios
    bar_str <- paste0("[", strrep("=", round(30 * progress_pct)), strrep(" ", 30 - round(30 * progress_pct)), "]")
    scenario_name <- paste0("Regime_", row, "_of_", total_scenarios)
    
    cat("\n====================================================================\n")
    cat(sprintf("STARTING: %s | PROGRESS: %s %d%%\n", scenario_name, bar_str, round(progress_pct * 100)))
    cat(sprintf("ESTIMATED TIME REMAINING: %s\n", eta_str))
    cat("====================================================================\n")
    
    current_bounds <- base_vars
    title_details <- c()
    for (var in names(chunk_intervals)) {
      chunk_index <- grid_indices[row, var]
      specific_chunk <- chunk_intervals[[var]][[chunk_index]]
      current_bounds[[var]] <- specific_chunk
      title_details <- c(title_details, sprintf("%s: [%.1f - %.1f]", var, specific_chunk[1], specific_chunk[2]))
    }
    
    var_names <- names(current_bounds)
    min_bounds <- sapply(current_bounds, function(x) x[1])
    max_bounds <- sapply(current_bounds, function(x) x[2])
    
    morris_obj <- morris(model = NULL, factors = var_names, r = trajectories, design = list(type = "oat", levels = 5, grid.jump = 1), binf = min_bounds, bsup = max_bounds)
    output_df <- run_batch(morris_obj$X, var_names)
    
    omnibus_data[[scenario_name]] <- list(
      obj = morris_obj,
      results = output_df,
      title_details = title_details
    )
    
    saveRDS(omnibus_data, database_file)
    cat(sprintf("   > Progress Auto-Saved to '%s'\n", database_file))
  }
  
  cat("\n>>> MORRIS OMNIBUS EXECUTION COMPLETE! <<<\n")
}

# ==============================================================================
# 7. POST-PROCESSING: MORRIS OMNIBUS AUTO-PLOTTER 
# ==============================================================================
all_targets <- c("Mixing_Potential", "Blended_Porosity", "Blended_Sphericity", 
                 "Blended_WetSkin", "Blended_Size_um", "Blended_Viscosity_PaS", 
                 "SG_Colloid", "Blended_SG_Gassed")

cat(sprintf("\n>>> Master Database loaded. Initiating Auto-Plot sequence for %d targets...\n", length(all_targets)))

for (target_feature in all_targets) {
  
  plot_df <- data.frame()
  valid_scenarios <- 0
  
  cat(sprintf("\n--- Processing Target: %s ---\n", target_feature))
  
  for (scen_name in names(omnibus_data)) {
    scen <- omnibus_data[[scen_name]]
    if (!target_feature %in% colnames(scen$results)) next
    
    y_data <- scen$results[[target_feature]]
    
    if (all(is.na(y_data))) {
      cat(sprintf("   - %s: Skipped (All NAs)\n", scen_name))
      next
    }
    
    valid_scenarios <- valid_scenarios + 1
    morris_temp <- scen$obj
    tell(morris_temp, y_data)
    
    mu_star <- apply(morris_temp$ee, 2, function(x) mean(abs(x), na.rm = TRUE))
    sigma <- apply(morris_temp$ee, 2, function(x) sd(x, na.rm = TRUE))
    
    max_mu <- max(mu_star, na.rm = TRUE)
    max_sig <- max(sigma, na.rm = TRUE)
    
    if (max_mu < 1e-5 && max_sig < 1e-5) {
      cat(sprintf("   - %s: Skipped (Completely static outcome)\n", scen_name))
      next
    }
    
    temp_df <- data.frame(
      Factor = names(mu_star),
      mu_star = mu_star,
      sigma = sigma,
      Regime = scen_name
    )
    
    for (detail in scen$title_details) {
      parts <- strsplit(detail, ": ")[[1]]
      var_name <- parts[1]
      val_range <- parts[2]
      temp_df[[paste0(var_name, "_Regime")]] <- val_range
    }
    
    plot_df <- rbind(plot_df, temp_df)
    cat(sprintf("   - %s: Appended to master map.\n", scen_name))
  }
  
  if (nrow(plot_df) > 0) {
    regime_cols <- grep("_Regime$", colnames(plot_df), value = TRUE)
    active_facets <- names(which(sapply(plot_df[regime_cols], function(col) length(unique(col)) > 1)))
    
    # FORCED LAYOUT LOGIC 
    # Prioritizing horizontal vs vertical layout configuration
    if ("v_tip_Regime" %in% active_facets && "C_solid_mass_Regime" %in% active_facets) {
      plot_facets <- c("v_tip_Regime", "C_solid_mass_Regime")
      facet_formula <- as.formula("C_solid_mass_Regime ~ v_tip_Regime")
    } else {
      active_facets <- active_facets[order(sapply(plot_df[active_facets], function(col) length(unique(col))), decreasing = TRUE)]
      plot_facets <- head(active_facets, 2)
      
      if (length(plot_facets) == 0) {
        agg_formula <- as.formula("cbind(mu_star, sigma) ~ Factor")
        facet_formula <- NULL
      } else if (length(plot_facets) == 2) {
        facet_formula <- as.formula(paste(plot_facets[2], "~", plot_facets[1]))
      } else {
        facet_formula <- as.formula(paste(". ~", plot_facets[1]))
      }
    }
    
    if (length(plot_facets) > 0) {
      agg_formula <- as.formula(paste("cbind(mu_star, sigma) ~ Factor +", paste(plot_facets, collapse = " + ")))
    }
    
    clean_df <- aggregate(agg_formula, data = plot_df, FUN = mean, na.rm = TRUE)
    
    # NUMERIC SORTING & TWO-LINE FACET LABELS
    for (col in plot_facets) {
      first_numbers <- as.numeric(gsub("\\[([0-9.]+).*", "\\1", clean_df[[col]]))
      clean_var_name <- gsub("_Regime", "", col)
      new_labels <- paste0(clean_var_name, "\n", clean_df[[col]])
      level_order <- unique(new_labels[order(first_numbers)])
      clean_df[[col]] <- factor(new_labels, levels = level_order)
    }
    
    # CLASSIFY SENSITIVITY ZONES 
    max_mu_panel <- max(clean_df$mu_star, na.rm = TRUE)
    max_sig_panel <- max(clean_df$sigma, na.rm = TRUE)
    
    clean_df$Zone <- "1. Negligible"
    active_mask <- (clean_df$mu_star > 0.025 * max_mu_panel) | (clean_df$sigma > 0.025 * max_sig_panel)
    
    clean_df$Zone[active_mask & (clean_df$sigma < clean_df$mu_star)] <- "2. Predictable / Additive"
    clean_df$Zone[active_mask & (clean_df$sigma >= clean_df$mu_star)] <- "3. Interacting / Complex"
    
    volatility_mask <- active_mask & (clean_df$mu_star < 0.10 * max_mu_panel) & (clean_df$sigma > 0.20 * max_sig_panel)
    clean_df$Zone[volatility_mask] <- "4. Pure Volatility"
    
    # THE MULTI-LENS GENERATOR 
    group_chem <- c("CMC", "HLB", "A_molecule", "MW_surfactant", "Delta_pH", "ionic_strength", "C_surfactant", "D_surf")
    group_add  <- c("C_binder", "C_plasticizer", "C_monomer")
    group_phys <- setdiff(unique(clean_df$Factor), c(group_chem, group_add)) 
    
    lenses <- list(
      "Chemistry" = group_chem,
      "Additives" = group_add,
      "Physical"  = group_phys
    )
    
    for (lens_name in names(lenses)) {
      lens_factors <- lenses[[lens_name]]
      
      lens_df <- clean_df
      lens_df$Zone[!(lens_df$Factor %in% lens_factors)] <- "1. Negligible"
      
      p <- ggplot(lens_df, aes(x = mu_star, y = sigma, color = Zone)) +
        geom_point(size = 3, alpha = 0.8) +
        geom_text_repel(data = subset(lens_df, Zone != "1. Negligible"),
                        aes(label = Factor), 
                        size = 3.5, 
                        fontface = "bold", 
                        box.padding = 0.5, 
                        max.overlaps = Inf,
                        min.segment.length = 0) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.6) +
        scale_color_manual(values = c("1. Negligible" = "grey85", "2. Predictable / Additive" = "#27ae60", 
                                      "3. Interacting / Complex" = "#f39c12", "4. Pure Volatility" = "#e74c3c")) +
        labs(title = paste("Regime Trellis Map:", target_feature, "| LENS:", lens_name), 
             subtitle = "Averaged Multi-Dimensional Sub-Regimes",
             x = expression("Main Effect (" * mu^"*" * ")"), y = expression("Interactions (" * sigma * ")")) +
        theme_bw() + theme(legend.position = "bottom", strip.background = element_rect(fill = "#b0c4de"),
                           strip.text = element_text(face = "bold", size = 10))
      
      if (!is.null(facet_formula)) {
        p <- p + facet_grid(facet_formula, scales = "free")
      }
      
      filename <- paste0("Lens_", lens_name, "_", target_feature, ".png")
      ggsave(filename, p, width = 14, height = 10, dpi = 300)
      cat(sprintf("   > Saved Lens Plot: %s\n", filename))
    }
    
  } else {
    cat(sprintf(">>> No valid variance recorded for %s. Skipping plot.\n", target_feature))
  }
}
cat("\n>>> ALL EXECUTIONS COMPLETE. <<<\n")
