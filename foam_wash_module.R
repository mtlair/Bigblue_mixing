# =============================================================================
# UP2 — FOAM-WASH COLUMN MODULE
# Algebraic closure derived from the standalone ODE (foam_wash_column_psd.R)
# Captures hindered settling, bimodal bubble population, surfactant film
# elasticity, burst logic, and per-class particle loss in steady-state form.
# =============================================================================

# Constants (match standalone model + interface_stream.R)
GRAV <- 9.81

# Stokes settling and aggregate density (reuse from standalone)
stokes  <- function(dRho, d, mu) dRho * GRAV * d^2 / (18 * mu)
rho_agg <- function(dp, db, rp, rg) (rp*dp^3 + rg*db^3) / (dp^3 + db^3)

# --- THERMODYNAMIC + SURFACTANT DERIVATION (steps 1-3) ----------------------
# Langmuir adsorption + Gibbs elasticity + film drainage/holdup closure
# All dimensionless factors are 1.0 at the baseline state
derive_foam_wash_state <- function(p) {
  Tk <- p[["T_col"]]

  # T/P thermodynamic state
  mu_T      <- p[["mu_ref"]] * exp(p[["E_visc"]] * (1/Tk - 1/p[["T_ref"]]))
  rho_gas_T <- p[["P_col_Pa"]] * p[["MW_gas"]] / (p[["R_gas"]] * Tk)
  sigma0_T  <- p[["sigma_ref"]] + p[["dsigma_dT"]] * (Tk - p[["T_ref"]])

  # Langmuir adsorption (monomer activity capped at the CMC)
  c_eff  <- min(p[["c_surf"]], p[["cmc"]])
  theta  <- p[["K_ads"]] * c_eff / (1 + p[["K_ads"]] * c_eff)
  theta  <- min(theta, 0.995)
  Gamma  <- p[["Gamma_inf"]] * theta
  sigma  <- sigma0_T + p[["R_gas"]] * Tk * p[["Gamma_inf"]] * log(1 - theta)
  # Gibbs elasticity with micelle buffering above the CMC
  micelle_buffer <- 1 / (1 + p[["K_mic"]] * max(0, p[["c_surf"]] - p[["cmc"]]) / p[["cmc"]])
  E_gibbs <- p[["R_gas"]] * Tk * p[["Gamma_inf"]] * theta / (1 - theta) * micelle_buffer

  # Film stability (normalized ~1 at baseline)
  film_stab <- E_gibbs / p[["E_stab_ref"]]

  # MIXER (up1) birth law: Hinze turbulent breakup
  d_b_in_eff <- p[["d_b_in"]] * (p[["V_tip_ref"]] / p[["V_tip"]])^p[["n_hinze"]] *
                (sigma / p[["sigma_ref_film"]])^0.6
  p[["d_b_in"]] <- d_b_in_eff

  # Drainage: surface mobility + effective tau
  mobility      <- 0.5 * (1 + p[["mu_surf"]] / p[["mu_surf_ref"]])
  tau_drain_eff <- p[["tau_drain"]] * (mu_T / p[["mu_ref"]]) * mobility

  # Film thickness -> equilibrium holdup
  r_pb_dry  <- p[["c_pb"]] * p[["d_b_in"]] * sqrt(p[["eps_l_dry"]])
  P_c_dry   <- sigma / r_pb_dry
  Pi_e      <- p[["Pi_charge"]] * theta
  h_eq_dry  <- p[["lambda_D"]] * log(max(Pi_e / P_c_dry, 1 + 1e-9))
  film_wet  <- h_eq_dry / p[["h_eq_ref"]]
  eps_l_dry_eff <- min(p[["eps_l_in"]],
                       p[["eps_l_dry"]] * film_stab^p[["hold_exp"]] * film_wet^p[["q_film"]])

  # Store derived values
  p[["mu_cont"]]        <- mu_T
  p[["rho_gas"]]        <- rho_gas_T
  p[["film_stability"]] <- film_stab
  p[["tau_drain"]]      <- tau_drain_eff
  p[["eps_l_dry"]]      <- eps_l_dry_eff
  p[["sigma_film"]]     <- sigma

  p
}

# --- PARTICLE RETENTION (algebraic from ODE steady-state logic) ----------------
# Returns per-class retention fractions: c(fine, mid, coarse)
retention_fractions <- function(stream, p) {
  # Unpack stream
  rho_liquid <- stream$rho_slurry
  rho_p      <- stream$rho_polymer
  d_fine     <- p[["d_fine"]]
  d_mid      <- p[["d_mid"]]
  d_crs      <- p[["d_crs"]]

  # Load and local state (nominal)
  load_rel <- 1.0  # assume nominal loading for retention estimate
  eps_s_in <- p[["eps_s_in"]] * load_rel
  eps_l_in <- p[["eps_l_in"]]
  eps_g_in <- p[["eps_g_in"]]

  # Solve for the pool depth (where alpha_g transitions from pool to foam)
  # The ODE has H_pool as a parameter; we approximate the transition as a
  # volume balance: pool fraction ~ H_pool / H_total.
  H_pool   <- p[["H_pool"]]
  H_total  <- p[["H_total"]]
  pool_frac <- H_pool / H_total

  # Krieger-Dougherty local viscosity (at nominal state)
  phi_cond <- eps_s_in / (eps_s_in + eps_l_in)
  phi_smax <- p[["phi_smax"]]
  mu_eff   <- p[["mu_cont"]] * (1 - min(phi_cond, 0.999*phi_smax)/phi_smax)^(-2.5*phi_smax)

  # Hindered settling (Richardson-Zaki) in pool
  phi <- 1 - eps_g_in
  n_RZ <- p[["n_RZ"]]
  U_settle <- function(dp) {
    d_agg <- (dp^3 + p[["d_b_in"]]^3)^(1/3)
    U_stokes <- stokes(rho_agg(dp, p[["d_b_in"]], rho_p, p[["rho_gas"]]) - rho_liquid, d_agg, mu_eff)
    U_stokes * (1 - phi)^n_RZ / p[["H_pool"]]
  }

  # Pool settling loss (residence time ~ H_pool / U_settle)
  settle_loss <- function(dp) pmin(1.0, U_settle(dp))
  sl_fine <- settle_loss(d_fine)
  sl_mid  <- settle_loss(d_mid)
  sl_crs  <- settle_loss(d_crs)

  # Foam-zone particle loss (detachment + burst-driven collapse)
  # Weighted by foam residence time; default detachment rates from ODE
  lf_fine <- p[["k_det_fine"]]
  lf_mid  <- p[["k_det_mid"]]
  lf_crs  <- p[["k_det_crs"]]

  # Burst contribution (scales with bubble coalescence; Gast model:
  # burst rate ~ K_burst * max(0, d_b - d_b_crit))
  d_b_crit <- p[["d_b_burst"]] * p[["film_stability"]]^p[["a_fs"]]
  burst_rate <- max(0, p[["K_burst"]] * (p[["d_b_in"]] - d_b_crit) / p[["d_b_in"]])
  lf_fine <- lf_fine + p[["k_burstloss_fine"]] * burst_rate
  lf_mid  <- lf_mid  + p[["k_burstloss_mid"]]  * burst_rate
  lf_crs  <- lf_crs  + p[["k_burstloss_crs"]]  * burst_rate

  # Total loss: pool fraction settled + foam fraction detached/burst
  # Normalized per residence time
  U_up <- p[["U_up"]]
  t_foam <- (1 - pool_frac) * H_total / U_up
  t_pool <- pool_frac * H_total / U_up

  loss_fine <- pool_frac * sl_fine + (1 - pool_frac) * min(1.0, lf_fine * t_foam)
  loss_mid  <- pool_frac * sl_mid  + (1 - pool_frac) * min(1.0, lf_mid  * t_foam)
  loss_crs  <- pool_frac * sl_crs  + (1 - pool_frac) * min(1.0, lf_crs  * t_foam)

  # Retention = 1 - loss
  c(fine = pmax(0, 1 - loss_fine),
    mid  = pmax(0, 1 - loss_mid),
    crs  = pmax(0, 1 - loss_crs))
}

# --- MAIN INTERFACE: foam_wash_column(stream, pars) -> stream ----------------
foam_wash_column <- function(stream, pars = list()) {
  # Default parameters (baseline state from standalone ODE)
  p <- modifyList(list(
    # --- operational ---
    eta_gas      = 0.75,     # fraction of entrained gas washed out
    P_col_Pa     = 3.0e5,    # column operating pressure
    T_col        = 298.15,   # column temperature [K]
    V_tip        = 14.5,     # mixer tip speed [m/s] (fed from upstream mixer)
    D_b_feed_m   = 5.0e-5,   # residual bubble size at column inlet (pre-collapse)
    particle_loss = FALSE,   # include per-class particle retention? (off by default for baseline calibration)

    # --- column geometry ---
    H_total    = 5.0,      # total column height [m]
    H_pool     = 0.45,     # pool depth (below decant overflow)
    U_up       = 8.0e-4,   # upflow velocity [m/s]

    # --- upstream passed properties ---
    mu_ref     = 2.0e-3,   # viscosity at T_ref [Pa s]
    T_ref      = 298.15,   # reference temperature [K]
    E_visc     = 1800,     # Andrade activation temp [K]
    rho_liquid = 1050,     # continuous phase density [kg/m3]
    rho_gas    = 1.8,      # baseline gas density [kg/m3]
    MW_gas     = 0.029,    # molar mass [kg/mol]
    R_gas      = 8.314,    # gas constant [J/mol/K]

    # --- particle classes ---
    d_fine     = 1.5e-5,   # fine size [m]
    d_mid      = 8.0e-5,   # mid size [m]
    d_crs      = 3.0e-4,   # coarse size [m]

    # --- feed state ---
    eps_l_in   = 0.50,     # liquid holdup
    eps_s_in   = 0.020,    # solids holdup
    eps_g_in   = 0.48,     # gas holdup (inlet from mixer)
    phi_smax   = 0.64,     # max packing fraction
    n_RZ       = 4.65,     # Richardson-Zaki exponent

    # --- surface tension + surfactant (step 1) ---
    c_surf     = 5.0,      # bulk surfactant concentration [mol/m3]
    cmc        = 6.0,      # critical micelle concentration [mol/m3]
    Gamma_inf  = 4.0e-6,   # Langmuir plateau [mol/m2]
    K_ads      = 2.0,      # Langmuir adsorption constant [m3/mol]
    K_mic      = 1.0,      # micelle buffering coefficient [-]
    sigma_ref  = 0.045,    # clean-solvent surface tension at T_ref [N/m]
    dsigma_dT  = -1.5e-4,  # temperature coefficient [N/m/K]

    # --- film elasticity (step 2) ---
    E_stab_ref = 0.099,    # reference Gibbs elasticity [N/m]

    # --- film drainage (step 3) ---
    eps_l_dry  = 0.15,     # baseline equilibrium holdup
    tau_drain  = 2000,     # drainage timescale [s]
    hold_exp   = 0.6,      # holdup sensitivity to film stability
    mu_surf    = 1.0e-6,   # surface viscosity [N s/m]
    mu_surf_ref = 1.0e-6,  # reference surface viscosity [N s/m]
    c_pb       = 0.20,     # Plateau-border radius coefficient
    q_film     = 0.5,      # holdup sensitivity to film thickness
    lambda_D   = 1.0e-8,   # Debye length [m]
    Pi_charge  = 1.0e4,    # electrostatic disjoining-pressure scale [Pa]
    h_eq_ref   = 4.195e-8, # reference equilibrium film thickness [m]
    sigma_ref_film = 0.02122, # baseline film surface tension [N/m]

    # --- mixer birth law (Hinze) ---
    V_tip_ref  = 14.5,     # reference tip speed [m/s]
    n_hinze    = 1.2,      # Hinze exponent [-]
    d_b_in     = 2.0e-3,   # born bubble diameter [m]
    d_b_fine_ref = 0.3e-3, # reference fine-mode bubble size [m]
    frac_gas_coarse_ref = 0.40, # inlet gas fraction in coarse mode

    # --- bubble population (coalescence-burst) ---
    K_coal     = 1.3e-4,   # coalescence coefficient [1/s]
    K_break    = 1.2e-4,   # breakage coefficient [1/s]
    d_b_burst  = 3.0e-3,   # critical burst size [m]
    K_burst    = 1.4e-3,   # burst-rate coefficient [1/s]
    a_fs       = 1.0,      # burst size sensitivity to film stability

    # --- particle loss rates ---
    k_det_fine = 1.2e-5,   # fine detachment [1/s]
    k_det_mid  = 6.0e-5,   # mid detachment [1/s]
    k_det_crs  = 3.0e-4,   # coarse detachment [1/s]
    k_burstloss_fine = 0.05, # fine burst-loss coefficient
    k_burstloss_mid = 0.30,  # mid burst-loss coefficient
    k_burstloss_crs = 2.0    # coarse burst-loss coefficient
  ), pars)

  # Derive thermodynamic + surfactant state
  p <- derive_foam_wash_state(p)

  # Gas removal ("washing the foam out")
  stream$alpha_g <- max(0, stream$alpha_g * (1 - p$eta_gas))

  # Pressurized -> Boyle compression of residual bubbles + set bubble field
  P_up <- if (is.finite(stream$P_Pa) && stream$P_Pa > 0) stream$P_Pa else 1.013e5
  stream$D_b_m <- p$D_b_feed_m * (P_up / p$P_col_Pa)^(1/3)
  stream$P_Pa  <- p$P_col_Pa

  # Per-class particle retention (wash-driven dropout) — optional
  # Default off to match handoff baseline; enable for realistic wash column physics
  if (p$particle_loss == TRUE) {
    ret <- retention_fractions(stream, p)
    stream$C_solid <- stream$C_solid * as.numeric(ret["mid"])
  }

  stream
}
