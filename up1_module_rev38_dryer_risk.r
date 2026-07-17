# ==============================================================================
# INTEGRATED SENSITIVITY ANALYSIS CONTROLLER & PLOTTER (v52 - Dryer Risk)
# ==============================================================================
# rev38 (on top of rev37): spray-dryer collapse risk outputs.
#   * Swelling_Softness_exit — bulk matrix softness at mixer exit; the same
#     exp(25*(C_monomer+C_plasticizer)) computed in the post-processor but now
#     returned as an explicit output for the downstream spray-dryer thread.
#   * Residual_Template_Fraction — fraction of loaded template volume that has
#     diffused into particle cores by mixer exit (spray-dryer collapse risk).
#     Physics: Fickian short-time penetration into a sphere:
#       tau_D = r^2 / (D_ref * Swelling_Softness * T_scale^1.5)
#       f_diff = 1 - exp(-3 * tau_eff / tau_D)
#     Gated by template_type: only surface_weld (type 3) has significant core
#     penetration (polymer-soluble); capillary_bridge (type 4) stays mostly
#     interstitial; rigid/gas (1/2) carry no liquid template.
#   * Deterministic collapse-risk dose sweep: Swelling_Softness_exit and
#     Residual_Template_Fraction vs template_dose (rev38_Dryer_Collapse_Risk.png).
#   * Morris screen extended to the two new outputs.
#
# Design rationale (2026-07-17 discussion):
#   Spray-dryer failure = "hard shell around soft interior" with confined vapor.
#   If template stays interstitial (pendular dose, theta-solvent choice) the
#   pore network stays open and the failure geometry never forms. If template
#   has diffused into particle cores, cores are plasticized and the volatile
#   solvent flashes under confinement → inflation or rupture.
#   Theta solvent (chi ≈ 0.5): D_eff ~ 10-100× lower than a good solvent;
#   tau_D >> mixer residence time → cores intact. Monomer (chi ≈ 0): D_eff
#   maximized, sub-minute penetration → worst case.
#
# rev37 (on top of rev36): liquid/gas/rigid TEMPLATING strategy study.
#   * template_type switch: 1=rigid seed, 2=gas (non-penetrating), 3=surface-
#     weld liquid (polymer-soluble, volatile -> Frenkel neck welding),
#     4=capillary-bridge liquid (both-immiscible, wetting -> Rumpf pendular
#     bridges). See Liquid_Template_Summary.md for the physics rationale.
#   * template_dose factor = template-liquid volume as a fraction of the
#     interstitial VOID volume of the packed solids (pendular 0-0.25,
#     funicular 0.25-0.80, capillary 0.80-1.0).
#   * New product-property outputs (terminal algebraic, no new ODE states):
#       Bond_Strength         - neck (Frenkel) or bridge (Rumpf) bond, 0..1
#       Retained_Porosity     - dry-product void fraction after strip/extract
#       Template_Target_Score - Bond_Strength * (Retained_Porosity/eps_rcp),
#                               peaks in the pendular window (the design goal:
#                               particles stick a bit, voids preserved)
#   * Deterministic pendular DOSE SWEEP (Pendular_Dose_Window.csv / .png).
#   * Morris screen per template_type; plots grouped into Process / Surface-
#     Chemistry / Polymer lenses with zone colours (green=linear/additive,
#     light-orange=complex/interacting, red=volatile, grey=no interaction)
#     and FIXED per-feature axes so variables are directly comparable.
#
# rev36 (on top of rev35): sweep-design improvements driven by the pooled
# range audit of the rev35 Morris omnibus:
#   a. v_tip lower bound extended to 0.5 m/s; non-uniform chunk breaks
#      [0.5-4 | 4-12 | 12-32] to resolve the cavern-collapse transition that
#      lives almost entirely below ~4 m/s (92% of rev35 runs were saturated).
#   b. tau promoted to a sweep dimension (3 chunks) - it is the top-ranked
#      driver of Mixing_Potential / Porosity / SG_Gassed but was not chunked.
#   c. ionic_strength: 3 non-uniform chunks [0.01-0.05 | 0.05-0.15 | 0.15-0.5]
#      to resolve the steep Debye (sqrt-I) end.
#   d. C_surfactant: bounds aligned to the base_vars dictionary (5e-5, was
#      inconsistently 1e-5 in the sweep config), reduced to 2 chunks split at
#      1e-3 (the bottom of the CMC range), since the top of its range is
#      CMC-capped in ~46% of runs.
#   e. Log-uniform sampling for decade-spanning factors MW_surfactant, D_surf,
#      C_surfactant (linear sampling put ~90% of points in the top decade).
#      Their Morris elementary effects are therefore per log10-unit.
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
force_recalculate <- FALSE     
database_file     <- sprintf("Morris_Templating_rev38_r%d.rds", trajectories)
plot_prefix       <- "rev38_"
# Factors sampled log-uniformly (Morris design operates in log10 space):
log_vars          <- c("MW_surfactant", "D_surf", "C_surfactant")

# ==============================================================================
# 3. VARIABLE DICTIONARY (Transition Zone Bounds)
# ==============================================================================
base_vars <- list(
  v_tip             = c(0.5, 32.0),
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
  Delta_pH          = c(0.0, 5.0),
  template_dose     = c(0.02, 1.10)      # template liquid / interstitial void volume
)

static_equipment <- list(
  chi_npgas         = 0.0,
  D_imp2D_tank      = 0.71,
  N_impellers       = 4.0,
  template_type     = 1      # 1=rigid, 2=gas, 3=surface_weld, 4=capillary_bridge
)

# Templating strategies swept in the outer study loop (Section 6)
template_modes <- c(rigid = 1, gas = 2, surface_weld = 3, capillary_bridge = 4)

# Pendular / funicular / capillary saturation boundaries (fraction of void)
PEND_MAX <- 0.25
CAP_MAX  <- 0.80
EPS_RCP  <- 0.36    # random-close-packing void fraction of the packed solids

# ==============================================================================
# 4. CHUNKING CONFIGURATION 
# ==============================================================================
# Entries may specify either equal-width chunks (bounds + n_chunks) or
# explicit non-uniform chunk edges (breaks).
sweep_targets <- list(
  v_tip          = list(breaks = c(0.5, 4.0, 12.0, 32.0)),
  tau            = list(bounds = c(5.0, 40.0),  n_chunks = 3),
  C_solid_mass   = list(bounds = c(0.10, 0.60), n_chunks = 3),
  Delta_pH       = list(bounds = c(0.0, 5.0),   n_chunks = 3),
  ionic_strength = list(breaks = c(0.01, 0.05, 0.15, 0.50)),
  C_surfactant   = list(breaks = c(5e-5, 1e-3, 0.015))
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

  # Templating mode: gas template activates the non-penetrating gas path so
  # template loading is routed into trapped gas rather than rigid inventory.
  t_type <- if (!is.null(p$template_type)) p$template_type else 1
  if (t_type == 2) p$chi_npgas <- 1.0

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
  
  if (is.null(ode_out) || nrow(ode_out) < 2) return(rep(NA, 13))
  
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

  # ==========================================================================
  # TEMPLATING CLOSURE (rev37): bond strength, retained porosity, target score
  # See Liquid_Template_Summary.md for the physics. All terminal-algebraic.
  # ==========================================================================
  fill <- s_pos(p$template_dose)   # design var = template liquid / void volume

  # Pendular saturation shape: bonds rise with fill and saturate near PEND_MAX
  sat <- tanh(fill / (0.5 * PEND_MAX))

  # Surfactant lowers interfacial tension -> weaker capillary/weld driving force
  gamma_proxy <- s_max(0.1, 1.0 - 0.8 * Theta_surf)
  # Rumpf: agglomerate strength ~ 1/d_particle (smaller spheres bond stronger)
  size_factor <- s_min(2.0, 0.7 / s_max(0.1, p$D_particle))
  # Frenkel neck welding accelerates with local surface softening (solvent+heat)
  # and available contact time; surface-weld does NOT feed the bulk softness pool
  weld_soft <- s_min(3.0, Swelling_Softness * (p$T_system / 298.15) * (p$tau / 20.0))

  # Seed-deposition tack for the rigid case (low, from R_het surface coverage)
  A_template_f <- s_pos(final_state$n_temp_A) * pi * (p$D_template * 1e-6)^2
  A_latex_f    <- (6.0 * phi_s) / s_pos(p$D_particle * 1e-6)
  Seed_Frac_f  <- A_template_f / s_pos(A_template_f + A_latex_f)

  bond_rigid <- s_clamp(0.15 * Seed_Frac_f + 0.05 * sat, 0, 0.4)
  bond_gas   <- 0.02 + 0.0 * fill                       # gas templates do not glue
  bond_weld  <- s_clamp(sat * s_min(1.0, 0.6 * weld_soft * gamma_proxy), 0, 1)
  bond_cap   <- s_clamp(sat * s_min(1.0, 0.9 * gamma_proxy * size_factor), 0, 1)

  # Void retention (dry product). PEND window preserves voids; over-dose loses
  # them: surface-weld by coalescence, capillary by densification.
  coalesce_pen <- s_clamp(s_pos(fill - PEND_MAX) * weld_soft, 0, 1)
  densify_pen  <- s_clamp(s_pos(fill - CAP_MAX) / (1.0 - CAP_MAX), 0, 1)

  ret_rigid <- s_clamp(EPS_RCP * (1.0 - 0.6 * fill), 0, EPS_RCP)     # solid seed occupies void
  ret_gas   <- s_clamp(Blended_Porosity, 0, 0.99)                    # gas voids are the porosity
  ret_weld  <- s_clamp(EPS_RCP * (1.0 - coalesce_pen), 0, EPS_RCP)
  ret_cap   <- s_clamp(EPS_RCP * (1.0 - 0.4 * densify_pen), 0, EPS_RCP)

  Bond_Strength <- switch(as.character(t_type),
                          "1" = bond_rigid, "2" = bond_gas,
                          "3" = bond_weld,  "4" = bond_cap, bond_rigid)
  Retained_Porosity <- switch(as.character(t_type),
                          "1" = ret_rigid, "2" = ret_gas,
                          "3" = ret_weld,  "4" = ret_cap, ret_rigid)
  # Objective: strong-enough bonds AND preserved voids -> peaks in pendular window
  Template_Target_Score <- s_clamp(Bond_Strength * (Retained_Porosity / EPS_RCP), 0, 1)

  # ==========================================================================
  # SPRAY-DRYER COLLAPSE RISK (rev38)
  # Residual_Template_Fraction: fraction of loaded template that has diffused
  # into particle cores by mixer exit.  High value means cores are plasticised
  # and contain volatile solvent -> inflation / rupture in the dryer.
  #
  # Fickian penetration into a sphere (short-time):
  #   tau_D = r_particle^2 / D_polymer_eff      (diffusion timescale)
  #   f_diff = 1 - exp(-3 * tau / tau_D)        (fraction reached core)
  # D_polymer_eff = D_ref * Swelling_Softness * (T/298)^1.5
  #   D_ref = 2e-13 m^2/s  (theta-solvent diffusivity in dry polymer at 298 K;
  #   good-solvent D is 10-100x higher, captured by Swelling_Softness proxy)
  #
  # Type gating:
  #   rigid (1), gas (2) : no liquid template -> 0
  #   surface_weld  (3)  : polymer-soluble -> full f_diff
  #   capillary_bridge(4): polymer-immiscible, stays interstitial;
  #                        small residual wetting fraction (0.05 * f_diff)
  # ==========================================================================
  r_particle_m  <- (p$D_particle * 1e-6) / 2.0
  # D_polymer_ref: theta-solvent diffusivity in near-Tg latex at 298 K.
  # For a good solvent (monomer, chi~0): D ~ 1e-13 to 1e-12 m^2/s -> tau_D < 1 s
  #   -> 100% penetration at any mixer tau (captured via Swelling_Softness path).
  # For theta solvent (chi~0.5): D ~ 5e-17 m^2/s -> tau_D ~ 40 min at r=0.6 µm
  #   -> partial penetration; tau and D_particle are the key controlling variables.
  # Swelling_Softness (1.2 - 5.8x over model range) scales D_eff upward when
  # pre-existing monomer/plasticizer has already opened the polymer network.
  D_polymer_ref <- 5e-17                             # m^2/s, theta condition, 298 K
  T_scale_38    <- (p$T_system / 298.15)^1.5        # Arrhenius-like scaling
  D_polymer_eff <- D_polymer_ref * Swelling_Softness * T_scale_38  # always > 0
  tau_D_s       <- (r_particle_m^2) / D_polymer_eff               # always > 0; no s_pos (would swamp 1e-16)
  f_diff        <- s_clamp(1.0 - exp(-3.0 * (p$tau * 60.0) / tau_D_s), 0.0, 1.0)  # p$tau in min -> sec

  Swelling_Softness_exit <- Swelling_Softness        # alias for downstream thread

  Residual_Template_Fraction <- switch(as.character(t_type),
    "1" = 0.0,
    "2" = 0.0,
    "3" = f_diff,
    "4" = s_clamp(0.05 * f_diff * s_clamp(fill, 0.0, 1.0), 0.0, 0.15),
    0.0)

  return(c(Blended_Porosity = Blended_Porosity,
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
           Residual_Template_Fraction = Residual_Template_Fraction))
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
  
  out_names <- c("Blended_Porosity", "Blended_Sphericity", "Mixing_Potential", "Blended_WetSkin",
                 "Blended_Size_um", "Blended_Viscosity_PaS", "SG_Colloid", "Blended_SG_Gassed",
                 "Bond_Strength", "Retained_Porosity", "Template_Target_Score",
                 "Swelling_Softness_exit", "Residual_Template_Fraction")
  n_out <- length(out_names)
  safe_results <- lapply(raw_results, function(x) { if (is.numeric(x) && length(x) == n_out) return(x) else return(rep(NA, n_out)) })
  res_df <- as.data.frame(do.call(rbind, safe_results), row.names = FALSE)
  colnames(res_df) <- out_names
  return(res_df)
}

format_time <- function(s) {
  if (is.na(s) || is.infinite(s) || s < 0) return("Calculating...")
  sprintf("%02d:%02d:%02d", as.integer(s %/% 3600), as.integer((s %% 3600) %/% 60), as.integer(s %% 60))
}

# ==============================================================================
# 6. TEMPLATING STUDY: PENDULAR DOSE SWEEP + MORRIS PER TEMPLATE TYPE
# ==============================================================================

# --- Variable groupings for the lens plots (rev37 requested split) ---------
group_process <- c("v_tip", "tau", "P_system", "T_system", "C_gas_diss_temp",
                   "Q_gas", "Q_colloid", "Q_template")
group_surfchem <- c("C_surfactant", "CMC", "HLB", "MW_surfactant", "A_molecule",
                    "D_surf", "ionic_strength", "Delta_pH")
group_polymer <- c("C_solid_mass", "C_temp_mass", "D_particle", "D_template",
                   "C_binder", "C_monomer", "C_plasticizer", "template_dose")
var_groups <- list("Process" = group_process,
                   "Surface Chemistry" = group_surfchem,
                   "Polymer" = group_polymer)
group_of <- function(f) {
  for (g in names(var_groups)) if (f %in% var_groups[[g]]) return(g)
  "Polymer"
}

all_targets <- c("Template_Target_Score", "Bond_Strength", "Retained_Porosity",
                 "Blended_Porosity", "Blended_Sphericity", "Mixing_Potential",
                 "Blended_WetSkin", "Blended_Size_um", "Blended_Viscosity_PaS",
                 "SG_Colloid", "Blended_SG_Gassed",
                 "Swelling_Softness_exit", "Residual_Template_Fraction")

# --------------------------------------------------------------------------
# 6a. DETERMINISTIC PENDULAR DOSE SWEEP
#     Fix formulation/process at nominal midpoints; sweep template_dose for
#     each template type; record Bond_Strength / Retained_Porosity / Score.
# --------------------------------------------------------------------------
cat("\n>>> Pendular dose-window sweep...\n")
nominal <- sapply(base_vars, function(b) mean(b))
dose_grid <- seq(0.02, 1.10, length.out = 60)
dose_rows <- list()
for (mode_name in names(template_modes)) {
  static_equipment$template_type <- as.numeric(template_modes[[mode_name]])
  for (d in dose_grid) {
    pars <- nominal; pars["template_dose"] <- d
    r <- run_single_simulation(pars)
    dose_rows[[length(dose_rows) + 1]] <- data.frame(
      template = mode_name, dose = d,
      Bond_Strength = r[["Bond_Strength"]],
      Retained_Porosity = r[["Retained_Porosity"]],
      Template_Target_Score = r[["Template_Target_Score"]],
      Swelling_Softness_exit = r[["Swelling_Softness_exit"]],
      Residual_Template_Fraction = r[["Residual_Template_Fraction"]])
  }
}
static_equipment$template_type <- 1
dose_df <- do.call(rbind, dose_rows)
write.csv(dose_df, "Pendular_Dose_Window.csv", row.names = FALSE)
cat("   > Saved Pendular_Dose_Window.csv\n")

# Locate the pendular optimum (max score) per template type
opt_tbl <- do.call(rbind, lapply(split(dose_df, dose_df$template), function(s) {
  i <- which.max(s$Template_Target_Score)
  data.frame(template = s$template[1], opt_dose = round(s$dose[i], 3),
             max_score = round(s$Template_Target_Score[i], 3))
}))
cat("   Pendular optima (fraction of void):\n"); print(opt_tbl, row.names = FALSE)

dose_long <- rbind(
  data.frame(template = dose_df$template, dose = dose_df$dose,
             metric = "Bond_Strength", value = dose_df$Bond_Strength),
  data.frame(template = dose_df$template, dose = dose_df$dose,
             metric = "Retained_Porosity", value = dose_df$Retained_Porosity),
  data.frame(template = dose_df$template, dose = dose_df$dose,
             metric = "Template_Target_Score", value = dose_df$Template_Target_Score))
dose_long$metric <- factor(dose_long$metric,
  levels = c("Bond_Strength", "Retained_Porosity", "Template_Target_Score"))

p_dose <- ggplot(dose_long, aes(x = dose, y = value, color = template)) +
  annotate("rect", xmin = 0, xmax = PEND_MAX, ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "#2ecc71") +
  annotate("rect", xmin = PEND_MAX, xmax = CAP_MAX, ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "#f39c12") +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(title = "Pendular Dose Window: template loading vs bond / void retention / target score",
       subtitle = "Green band = pendular (target); orange band = funicular; beyond = capillary/densified",
       x = "template_dose (template liquid / interstitial void volume)", y = NULL,
       color = "Template type") +
  theme_bw() + theme(legend.position = "bottom",
                     strip.background = element_rect(fill = "#b0c4de"),
                     strip.text = element_text(face = "bold"))
ggsave("rev38_Pendular_Dose_Window.png", p_dose, width = 11, height = 10, dpi = 200)
cat("   > Saved rev38_Pendular_Dose_Window.png\n")

# Collapse-risk plot: Swelling_Softness_exit and Residual_Template_Fraction vs dose
collapse_long <- rbind(
  data.frame(template = dose_df$template, dose = dose_df$dose,
             metric = "Swelling_Softness_exit",
             value = dose_df$Swelling_Softness_exit),
  data.frame(template = dose_df$template, dose = dose_df$dose,
             metric = "Residual_Template_Fraction",
             value = dose_df$Residual_Template_Fraction))
collapse_long$metric <- factor(collapse_long$metric,
  levels = c("Swelling_Softness_exit", "Residual_Template_Fraction"))

p_collapse <- ggplot(collapse_long, aes(x = dose, y = value, color = template)) +
  annotate("rect", xmin = 0,        xmax = PEND_MAX, ymin = -Inf, ymax = Inf,
           alpha = 0.12, fill = "#2ecc71") +
  annotate("rect", xmin = PEND_MAX, xmax = CAP_MAX,  ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "#f39c12") +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  labs(title = "Spray-dryer collapse risk vs template dose",
       subtitle = paste0("Green band = pendular (target, low risk); orange = funicular.\n",
         "Residual_Template_Fraction > 0 only for surface_weld (polymer-soluble template)."),
       x = "template_dose (template liquid / interstitial void volume)", y = NULL,
       color = "Template type") +
  theme_bw() + theme(legend.position = "bottom",
                     strip.background = element_rect(fill = "#f9c784"),
                     strip.text = element_text(face = "bold"))
ggsave("rev38_Dryer_Collapse_Risk.png", p_collapse, width = 11, height = 8, dpi = 200)
cat("   > Saved rev38_Dryer_Collapse_Risk.png\n")

# --------------------------------------------------------------------------
# 6b. MORRIS SCREEN PER TEMPLATE TYPE (full-range single regime per type)
# --------------------------------------------------------------------------
if (file.exists(database_file) && !force_recalculate) {
  cat(sprintf("\n>>> Loading cached templating Morris DB '%s'...\n", database_file))
  omnibus_data <- readRDS(database_file)
} else {
  cat("\n>>> Morris screen per template type...\n")
  var_names <- names(base_vars)
  min_bounds <- sapply(base_vars, function(x) x[1])
  max_bounds <- sapply(base_vars, function(x) x[2])
  # log-uniform factors sampled in log10 space
  dmin <- min_bounds; dmax <- max_bounds
  for (lv in intersect(log_vars, var_names)) { dmin[lv] <- log10(dmin[lv]); dmax[lv] <- log10(dmax[lv]) }

  omnibus_data <- list()
  for (mode_name in names(template_modes)) {
    cat(sprintf("   - Template: %s\n", mode_name))
    static_equipment$template_type <- as.numeric(template_modes[[mode_name]])
    set.seed(2024)
    mobj <- morris(model = NULL, factors = var_names, r = trajectories,
                   design = list(type = "oat", levels = 5, grid.jump = 1),
                   binf = dmin, bsup = dmax)
    Xmodel <- mobj$X
    for (lv in intersect(log_vars, var_names)) Xmodel[, lv] <- 10^Xmodel[, lv]
    out_df <- run_batch(Xmodel, var_names)
    omnibus_data[[mode_name]] <- list(obj = mobj, results = out_df, template = mode_name)
  }
  static_equipment$template_type <- 1
  saveRDS(omnibus_data, database_file)
  cat(sprintf("   > Saved %s\n", database_file))
}

# --------------------------------------------------------------------------
# 6c. AGGREGATE + ZONE CLASSIFY (per target, per template type, per factor)
# --------------------------------------------------------------------------
# Zone colours (rev37 request):
#   grey  = no interaction (negligible)
#   green = linear (additive, sigma < mu*)
#   light orange = complex (interacting, sigma >= mu*)
#   red   = volatile (low mu*, high sigma)
zone_levels <- c("No interaction", "Linear", "Complex", "Volatile")
zone_cols <- c("No interaction" = "grey80", "Linear" = "#27ae60",
               "Complex" = "#f5b041", "Volatile" = "#e74c3c")

classify_zone <- function(mu, sg, mu_max, sg_max) {
  active <- (mu > 0.05 * mu_max) | (sg > 0.05 * sg_max)
  z <- ifelse(!active, "No interaction",
       ifelse(sg < mu, "Linear", "Complex"))
  vol <- active & (mu < 0.15 * mu_max) & (sg > 0.30 * sg_max)
  z[vol] <- "Volatile"
  factor(z, levels = zone_levels)
}

build_target_df <- function(tg) {
  rows <- list()
  for (mode_name in names(omnibus_data)) {
    s <- omnibus_data[[mode_name]]
    y <- s$results[[tg]]
    if (is.null(y) || all(is.na(y))) next
    m <- s$obj; tell(m, y)
    mu <- apply(m$ee, 2, function(x) mean(abs(x), na.rm = TRUE))
    sg <- apply(m$ee, 2, function(x) sd(x, na.rm = TRUE))
    rows[[mode_name]] <- data.frame(template = mode_name, Factor = names(mu),
                                    mu_star = as.numeric(mu), sigma = as.numeric(sg))
  }
  if (length(rows) == 0) return(NULL)
  df <- do.call(rbind, rows)
  df$Group <- factor(sapply(as.character(df$Factor), group_of), levels = names(var_groups))
  df$template <- factor(df$template, levels = names(template_modes))
  df
}

# --------------------------------------------------------------------------
# 6d. PLOTS: one figure per target; rows = template type, cols = variable
#     group; FIXED axes across the whole figure so factors are comparable.
# --------------------------------------------------------------------------
cat("\n>>> Rendering grouped zone plots (fixed per-feature axes)...\n")
for (tg in all_targets) {
  df <- build_target_df(tg)
  if (is.null(df) || nrow(df) == 0) { cat(sprintf("   - %s: no data\n", tg)); next }
  mu_max <- max(df$mu_star, na.rm = TRUE); sg_max <- max(df$sigma, na.rm = TRUE)
  if (mu_max < 1e-12 && sg_max < 1e-12) { cat(sprintf("   - %s: static\n", tg)); next }
  df$Zone <- classify_zone(df$mu_star, df$sigma, mu_max, sg_max)
  xlim_hi <- mu_max * 1.15 + 1e-12
  ylim_hi <- sg_max * 1.15 + 1e-12

  p <- ggplot(df, aes(x = mu_star, y = sigma, color = Zone)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
    geom_point(size = 2.6, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = subset(df, Zone != "No interaction"),
      aes(label = Factor), size = 2.7, fontface = "bold",
      max.overlaps = Inf, box.padding = 0.4, min.segment.length = 0, show.legend = FALSE) +
    scale_color_manual(values = zone_cols, drop = FALSE) +
    facet_grid(template ~ Group) +
    coord_cartesian(xlim = c(0, xlim_hi), ylim = c(0, ylim_hi)) +
    labs(title = paste("Morris templating screen:", tg),
         subtitle = "Fixed axes across all panels. Colour: green=linear, orange=complex, red=volatile, grey=no interaction",
         x = expression("Main effect  " * mu^"*"), y = expression("Interactions  " * sigma),
         color = "Zone") +
    theme_bw() + theme(legend.position = "bottom",
                       strip.background = element_rect(fill = "#b0c4de"),
                       strip.text = element_text(face = "bold", size = 9))
  fn <- paste0(plot_prefix, "Morris_", tg, ".png")
  ggsave(fn, p, width = 13, height = 12, dpi = 180)
  cat(sprintf("   > %s\n", fn))
}

# --------------------------------------------------------------------------
# 6e. TEXT SUMMARY of rankings per template type for the target score
# --------------------------------------------------------------------------
sink("rev38_templating_summary.md")
cat("# rev38 Templating + Dryer Risk Morris Summary\n\n")
cat(sprintf("Trajectories r=%d; factors=%d; template types=%d.\n\n",
            trajectories, length(base_vars), length(template_modes)))
cat("## Pendular dose optima (Template_Target_Score vs template_dose)\n\n")
cat("| Template | Optimal dose (frac of void) | Max score |\n|---|---|---|\n")
for (i in seq_len(nrow(opt_tbl)))
  cat(sprintf("| %s | %.3f | %.3f |\n", opt_tbl$template[i], opt_tbl$opt_dose[i], opt_tbl$max_score[i]))
cat("\n## Top factors for Template_Target_Score by template type\n\n")
df_sc <- build_target_df("Template_Target_Score")
for (mode_name in names(template_modes)) {
  s <- df_sc[df_sc$template == mode_name, ]
  if (nrow(s) == 0) next
  s <- s[order(-s$mu_star), ]
  mm <- max(s$mu_star, 1e-30)
  cat(sprintf("**%s**: ", mode_name))
  cat(paste(sprintf("%s(%.2f,%s)", head(s$Factor, 6), head(s$mu_star / mm, 6),
                    ifelse(head(s$sigma, 6) >= head(s$mu_star, 6), "cplx", "lin")), collapse = ", "), "\n\n")
}

cat("\n## Top factors for Residual_Template_Fraction by template type (spray-dryer collapse risk)\n\n")
cat("_Non-zero only for surface_weld (type 3) and minimally for capillary_bridge (type 4)._\n\n")
df_rtf <- build_target_df("Residual_Template_Fraction")
for (mode_name in names(template_modes)) {
  s <- df_rtf[df_rtf$template == mode_name, ]
  if (nrow(s) == 0) next
  s <- s[order(-s$mu_star), ]
  mm <- max(s$mu_star, 1e-30)
  if (mm < 1e-20) { cat(sprintf("**%s**: zero (no liquid template or immiscible).\n\n", mode_name)); next }
  cat(sprintf("**%s**: ", mode_name))
  cat(paste(sprintf("%s(%.2f,%s)", head(s$Factor, 6), head(s$mu_star / mm, 6),
                    ifelse(head(s$sigma, 6) >= head(s$mu_star, 6), "cplx", "lin")), collapse = ", "), "\n\n")
}
sink()
cat("\n>>> rev38 DRYER RISK STUDY COMPLETE.\n")
