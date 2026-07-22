# =============================================================================
# UP2 MODULE: effervescent-airblast spray drying + particle formation
# =============================================================================
# Physics from morris_sensitivity_analysis.R (spray_dry_model: pressurized
# hold -> two-phase conditioning -> effervescent-airblast exit -> secondary
# breakup -> dryer balance -> particle formation), refactored so that the
# FEED is no longer a set of independent screening factors but the stream
# object handed over by the upstream chain (UP1 mixer -> intermediate stages).
#
#   up2_run_dryer(feed, x) -> named numeric vector (19 outputs)
#
# feed : stream list (see unified/interface_stream.R)
# x    : named list/vector of dryer-side factors:
#        ALR, P_atom_air, P_feed, mdot_L, sigma, mu_L, T_dryer_in,
#        mdot_gas_dry, Y_in, t_hold, D_b, Tg_polymer, n_flow,
#        k_perm_mono, k_perm_plast, k_perm_bind, T_bp_solv
#
# What the stream replaces (vs the standalone screen):
#   alpha_g_0      <- feed$alpha_g          (mixer trapped-gas holdup)
#   C_solid_mass   <- feed$C_solid + feed$C_solid_rigid
#   rho_L          <- feed$rho_slurry
#   T_feed         <- feed$T_K              (mixer exit temperature)
#   C_mono/plas/bind/surf, Delta_pH, I_strength <- stream composition
#   phi_emulsion   <- feed$phi_templ_free   (RTF-depleted liquid template)
#   D_template     <- feed$D_template_um * 1e-6
#   a_prim         <- feed$D_particle_um/2  (was a 0.1 um constant)
#   d_ratio        <- D_particle / D_agg    (was a 0.10 constant)
#   rho_s          <- feed$rho_polymer      (was 1400; UP1 uses 1050)
#   MW_surf, A_molec, HLB                   (were constants; now formulation)
#   D_b            <- feed$D_b_m if set by an intermediate stage, else the
#                     screening factor x$D_b (transfer-line shear unknown yet)
#
# New cross-module couplings (rev38 "collapse risk" wired in):
#   * theta_skin seeding: the mixer wet-skin fraction partially survives
#     atomization (SKIN_SURVIVAL) and pre-seeds the dryer skin state.
#   * Core-absorbed template (feed$w_core = RTF * template load): acts as a
#     core plasticizer — enters phi_solvent (Flory-Huggins Softness), the
#     shell permeability exponent (plasticizer coefficient), and the Fox
#     residual-solvent load; and as a confined volatile it raises the
#     micro-explosion severity weight. This implements the rev38 design
#     intent: plasticized cores holding volatile solvent -> inflation /
#     rupture in the dryer.
# =============================================================================

up2_constants <- function() list(
  P_atm    = 1.013e5,   # ambient pressure                        [Pa]
  R_air    = 287,       # gas constant, air                       [J/kg K]
  gamma_a  = 1.4,       # heat capacity ratio, air                [-]
  A_L      = 1.0e-6,    # liquid passage area                     [m2]
  D_h      = 1.0e-3,    # nozzle hydraulic diameter               [m]
  kB       = 1.380649e-23,
  N_Av     = 6.02214e23,
  phi_m    = 0.63,      # Krieger-Dougherty max packing           [-]
  Tg_solv  = 150,       # solvent (monomer/plasticizer) Tg, Fox   [K]
  C_cham   = 0.80,      # mixing chamber / air supply pressure    [-]
  k_rip0   = 4.0e-15,   # Ostwald ripening rate at 1 atm          [m3/s]
  cp_gas   = 1005,      # dryer gas heat capacity                 [J/kg K]
  h_fg     = 2.30e6,    # latent heat of water evaporation        [J/kg]
  gamma_ref = 100,      # rheometer reference shear rate          [1/s]
  t_res    = 8,         # dryer residence time                    [s]
  Tg_water = 138,       # glass transition of water (Fox)         [K]
  sig_y0   = 4.0e6,     # cake yield strength scale (Rumpf-type)  [Pa]
  rho_solv = 750,       # template solvent liquid density         [kg/m3]
  h_fg_solv = 3.5e5,    # template solvent latent heat            [J/kg]
  k_emu0   = 2.0e-20    # emulsion LSW ripening rate              [m3/s]
)

up2_output_names <- c("d_droplet_um", "d10_um", "d90_um", "d99_um", "span",
                      "BI_bimodal", "D_particle_um", "Dp90_um", "theta_skin_z",
                      "Omega_struct_z", "phi_porosity_z", "rho_tapped",
                      "Tg_product_K", "X_moisture", "Perm_shell_rel",
                      "D_pore_um", "solv_retained", "f_burst_solv",
                      "sigma_y_cake_MPa")

up2_run_dryer <- function(feed, x, cst = up2_constants()) {
  x <- as.list(x)

  ## --- Dryer-side factors ---------------------------------------------------
  ALR    <- x[["ALR"]];         P_G    <- x[["P_atom_air"]]
  P_F    <- x[["P_feed"]]
  mdot_L <- x[["mdot_L"]];      sigma  <- x[["sigma"]]
  mu_L   <- x[["mu_L"]]
  T_in   <- x[["T_dryer_in"]];  mdot_g <- x[["mdot_gas_dry"]]
  Y_in   <- x[["Y_in"]]
  t_hold <- x[["t_hold"]]
  Tg_pol <- x[["Tg_polymer"]];  n_fl   <- x[["n_flow"]]
  k_pm   <- x[["k_perm_mono"]]
  k_pp   <- x[["k_perm_plast"]]
  k_pb   <- x[["k_perm_bind"]]
  T_bpS  <- x[["T_bp_solv"]]

  ## --- Feed state from the stream ------------------------------------------
  alpha0 <- feed$alpha_g
  C_sol  <- feed$C_solid + feed$C_solid_rigid
  rho_L  <- feed$rho_slurry
  T_feed <- feed$T_K
  C_mono <- feed$C_monomer;  C_plas <- feed$C_plasticizer
  C_bind <- feed$C_binder;   C_surf <- feed$C_surfactant
  dpH    <- max(feed$Delta_pH, 0.05)
  I_str  <- feed$ionic_strength
  phi_e  <- feed$phi_templ_free
  D_e    <- feed$D_template_um * 1e-6
  w_core <- feed$w_core
  # Primary radius from the UP1 three-regime closure (milled size past the
  # critical tip speed); fall back to the input primary if the field is absent.
  # Prefer D_primary_phys_um (200 nm physical calibration bead) over
  # D_primary_exit_um (ODE-scaled, 1.25 µm base) for packing/Pe/surfactant.
  D_prim_um <- if (!is.null(feed$D_primary_phys_um) && !is.na(feed$D_primary_phys_um))
                 feed$D_primary_phys_um
               else if (!is.null(feed$D_primary_exit_um) && !is.na(feed$D_primary_exit_um))
                 feed$D_primary_exit_um
               else feed$D_particle_um
  a_prim <- (D_prim_um / 2) * 1e-6
  d_ratio <- min(max(D_prim_um / max(feed$D_agg_um, 1e-6), 0.01), 1.0)
  rho_s  <- feed$rho_polymer
  MW_surf <- feed$MW_surfactant / 1000          # g/mol -> kg/mol
  A_molec <- feed$A_molecule * 1e-18            # nm2 -> m2
  HLB     <- feed$HLB
  theta_seed <- SKIN_SURVIVAL * feed$WetSkin
  D_b <- if (is.finite(feed$D_b_m)) feed$D_b_m else x[["D_b"]]

  P_atm <- cst$P_atm; R_air <- cst$R_air; gamma_a <- cst$gamma_a
  A_L <- cst$A_L; D_h <- cst$D_h; kB <- cst$kB; N_Av <- cst$N_Av
  phi_m <- cst$phi_m; Tg_solv <- cst$Tg_solv; C_cham <- cst$C_cham
  k_rip0 <- cst$k_rip0; cp_gas <- cst$cp_gas; h_fg <- cst$h_fg
  gamma_ref <- cst$gamma_ref; t_res <- cst$t_res; Tg_water <- cst$Tg_water
  sig_y0 <- cst$sig_y0; rho_solv <- cst$rho_solv; h_fg_solv <- cst$h_fg_solv
  k_emu0 <- cst$k_emu0

  ## --- Module 0a: Formulation - Flory-Huggins free-volume swelling ---------
  # Core-absorbed template (w_core, from UP1's RTF) adds to the residual
  # solvent pool: plasticized cores are the rev38 collapse-risk pathway.
  phi_solvent <- C_mono + C_plas + w_core
  Softness    <- (1 + 5 * C_bind) * exp(25 * phi_solvent)

  # Shell permeability: core-absorbed template opens free volume like the
  # plasticizer channel (same permeation coefficient k_pp).
  Perm_shell <- exp(k_pm * C_mono + k_pp * (C_plas + w_core)) / (1 + k_pb * C_bind)

  # Free (interstitial) template emulsion rides as droplets to the dryer
  rho_lm <- (1 - phi_e) * rho_L + phi_e * rho_solv  # liquid mixture density
  w_e    <- phi_e * rho_solv / rho_lm               # free-solvent mass fraction

  ## --- Module 0b: Feed rheology (Krieger-Dougherty + power-law thinning) ---
  phi_s      <- C_sol * rho_L / rho_s             # solid phase volume fraction
  phi_disp   <- phi_s * (1 - phi_e) + phi_e       # solids + emulsion droplets
  mu_serum   <- mu_L * (1 + 10 * C_bind)          # binder thickening (K_fluid)
  mu_slurry0 <- mu_serum *
                (1 - min(phi_disp, 0.60) / phi_m)^(-2.5 * phi_m)
  mu_app <- function(gdot)
    pmax(1e-3, mu_slurry0 * (gdot / gamma_ref)^(n_fl - 1))

  ## --- Module 0c: Surfactant molar stoichiometry ---------------------------
  cap_area   <- (C_surf / MW_surf) * rho_L * N_Av * A_molec   # [m2/m3]
  part_area  <- 3 * phi_s * (1 - phi_e) / a_prim              # [m2/m3]
  gas_area   <- 6 * alpha0 / D_b                              # [m2/m3]
  emu_area   <- if (phi_e > 0) 6 * phi_e / D_e else 0         # [m2/m3]
  theta_surf <- min(1, cap_area / (part_area + gas_area + emu_area))
  # under-stabilized foam sheds entrained gas before the nozzle
  alpha_g    <- alpha0 * (0.3 + 0.7 * theta_surf)

  ## --- Module 0d: Pressurized hold - coalescence / Ostwald ripening --------
  k_rip <- k_rip0 * (P_F / P_atm) * (1 - 0.8 * theta_surf)
  D_b_h <- (D_b^3 + k_rip * t_hold)^(1/3)         # coarsened bubble size
  k_emu <- k_emu0 * (1 - 0.9 * theta_surf)
  D_e_h <- (max(D_e, 1e-8)^3 + k_emu * t_hold)^(1/3)  # coarsened template size

  ## --- Module 0e: DLVO electrostatics --------------------------------------
  E_rep     <- dpH * (1 / (1 + 10 * I_str)) * (1 - theta_surf)
  stability <- 1 + E_rep + 5 * theta_surf         # electrostatic + steric

  ## --- Module 1: Feed properties (at liquid-line pressure P_feed) ----------
  rho_G0    <- P_F / (R_air * T_feed)
  rho_eff0  <- (1 - alpha_g) * rho_lm + alpha_g * rho_G0
  mu_eff0   <- mu_slurry0 * (1 + 2.5 * alpha_g)
  sigma_eff <- sigma * (1 - 0.5 * alpha_g)

  ## --- Module 2: Nozzle pressure ladder & two-phase conditioning -----------
  P_cham   <- C_cham * P_G
  P_throat <- max(P_atm, 0.53 * P_cham)
  r_th     <- max(P_F / P_throat, 1)
  r_exp    <- max(P_F / P_atm, 1)
  r_flash  <- max(P_throat / P_atm, 1)
  alpha_th <- min(alpha_g * r_th / (1 - alpha_g + alpha_g * r_th), 0.97)
  alpha_e  <- min(alpha_g * r_exp / (1 - alpha_g + alpha_g * r_exp), 0.97)
  D_b_e    <- D_b_h * r_exp^(1/3)
  rho_Ge   <- P_atm / (R_air * T_feed)
  rho_eff  <- (1 - alpha_e) * rho_lm + alpha_e * rho_Ge

  ## --- Module 3: Nozzle hydraulics -----------------------------------------
  U_L <- mdot_L / (rho_lm * A_L)
  U_G <- 0.9 * sqrt(pmax(2 * gamma_a / (gamma_a - 1) * R_air * T_feed *
                         (1 - (P_atm / P_G)^((gamma_a - 1) / gamma_a)), 0))
  U_rel <- max(U_G - U_L, 5)
  rho_A <- P_atm / (R_air * T_feed)

  ## --- Module 4: Effervescent-airblast hybrid exit plane -------------------
  t_film <- (D_h / 2) * (1 - sqrt(alpha_th))
  L_c    <- max(2 * t_film, 0.02 * D_h)
  gdot_noz  <- U_rel / max(t_film, 1e-5)
  mu_eff_hs <- mu_app(gdot_noz) * (1 + 2.5 * alpha_g)
  SMD_ab <- 0.48 * L_c * (sigma_eff / (rho_A * U_rel^2 * L_c))^0.4 *
              (1 + 1 / ALR)^0.4 +
            0.15 * L_c * sqrt(mu_eff_hs^2 / (sigma_eff * rho_lm * L_c)) *
              (1 + 1 / ALR)
  D_b_th  <- D_b_h * r_th^(1/3)
  g_in    <- 1 / (1 + D_b_th / t_film)
  F_flash <- (1 + alpha_th * (r_flash - 1))^(-1/3)
  SMD     <- SMD_ab * (1 - g_in * (1 - F_flash))
  S_flashn <- 1 / (1 + exp(-(T_feed - T_bpS) / 5))
  F_fs     <- (1 + 40 * phi_e * S_flashn)^(-1/3)
  SMD      <- SMD * F_fs

  ## --- Module 5a: Secondary breakup (Weber-number correction) --------------
  We_g <- rho_A * U_rel^2 * SMD / sigma_eff
  f_sec <- 0.55 + 0.45 / (1 + We_g / 80)
  d_drop <- SMD * f_sec

  ## --- Module 5b: Droplet size distribution (three-mode volume mixture) ----
  E_a  <- ALR * U_rel^2
  Oh_e <- mu_eff_hs / sqrt(rho_lm * sigma_eff * L_c)
  w_c  <- min(0.6, 0.5 / (1 + E_a / 1.5e5) * (1 + 3 * Oh_e))
  w_f  <- min(0.35, 0.30 / (1 + 6e5 / E_a) +
                    0.15 * g_in * alpha_th +
                    0.50 * phi_e * S_flashn)
  if ((w_c + w_f) > 0.7) { sc <- 0.7 / (w_c + w_f); w_c <- w_c * sc; w_f <- w_f * sc }
  w_m  <- 1 - w_c - w_f

  m_main <- 1.2 * d_drop
  modes_w <- c(w_m, w_c, w_f)
  modes_m <- c(m_main, 3.2 * m_main, 0.25 * m_main)
  modes_s <- c(1.9, 1.5, 1.6)

  gr  <- seq(log(0.02 * m_main), log(20 * m_main), length.out = 400)
  cdf <- Reduce(`+`, Map(function(w, m, s)
           w * pnorm((gr - log(m)) / log(s)), modes_w, modes_m, modes_s))
  q_at <- function(p) exp(approx(cdf, gr, xout = p, ties = "ordered")$y)
  d10 <- q_at(0.10); d50 <- q_at(0.50); d90 <- q_at(0.90); d99 <- q_at(0.99)
  span <- (d90 - d10) / d50
  BI   <- 2 * min(w_f, 1 - w_f)

  ## --- Module 6a: Dryer energy / moisture balance (co-current) -------------
  mdot_w <- mdot_L * (1 - w_e) * (1 - C_sol)
  mdot_e <- mdot_L * w_e
  T_out  <- T_in - (mdot_w * h_fg + mdot_e * h_fg_solv) / (mdot_g * cp_gas)
  T_out  <- max(T_out, T_feed + 2)
  Y_out  <- Y_in + mdot_w / mdot_g
  p_v    <- Y_out / (Y_out + 0.622) * P_atm
  p_sat  <- 610.94 * exp(17.625 * (T_out - 273.15) / (T_out - 273.15 + 243.04))
  RH_out <- min(p_v / p_sat, 0.99)

  ## --- Module 6b: Drying kinetics -------------------------------------------
  kappa <- 5e-8 * pmax(T_out - T_feed, 1) / 100 * (1 - RH_out)
  D_diff <- kB * T_feed / (6 * pi * mu_L * a_prim)
  Pe <- kappa / (8 * D_diff)

  S_skin <- Pe * C_sol * (1 + 1 / stability)
  S_crit <- 500 * (1 + 0.05 * Softness)
  theta_skin <- S_skin / (S_skin + S_crit)
  # Cross-module coupling: the mixer wet-skin fraction that survives
  # atomization pre-seeds the shell (series combination of skin states)
  theta_skin <- 1 - (1 - theta_seed) * (1 - theta_skin)

  D_f <- 1.8 + 0.7 * min((stability - 1) / 5, 1)
  phi_struct <- 0.30 * theta_skin * (1 - d_ratio^(3 - D_f))

  ## --- Module 6c: Per-mode drying & residual moisture ----------------------
  tau_dry <- (modes_m^2 / kappa) * (1 + 20 * theta_skin / Perm_shell)
  X_j     <- exp(-t_res / tau_dry)
  X_moist <- sum(modes_w * X_j)

  ## --- Module 7a: Product glass transition (Fox: solvent + moisture) -------
  T_particle <- 0.85 * T_out + 0.15 * T_feed
  S_bs   <- 1 / (1 + exp(-(T_particle - T_bpS) / 8))
  R_solv <- 1 - S_bs

  # Residual solvent load now includes the core-absorbed template: it only
  # escapes above its boiling point (same R_solv gate as the free emulsion)
  w_res  <- (0.5 * C_mono + C_plas + (w_e + w_core) * R_solv) /
            (0.5 * C_mono + C_plas + (w_e + w_core) * R_solv + C_sol)
  w_wat  <- min(0.12, X_moist * (1 - C_sol) /
                      (X_moist * (1 - C_sol) + C_sol))
  Tg_dry <- 1 / ((1 - w_res) / Tg_pol + w_res / Tg_solv)
  Tg_eff <- 1 / ((1 - w_wat) / Tg_dry + w_wat / Tg_water)

  S_stick <- 1 / (1 + exp(-(T_particle - (Tg_eff + 20)) / 10))

  ## --- Module 7b: Cake (consolidated shell) mechanics -----------------------
  phi_cake <- 0.45 + 0.19 * min((stability - 1) / 5, 1)
  sigma_y  <- sig_y0 * (phi_cake / 0.64)^4 * (1 + 2 / stability) /
              sqrt(Softness) * (1 - 0.85 * S_stick)
  P_cap    <- 2 * sigma_eff / (0.3 * a_prim)
  Pi_col   <- P_cap / max(sigma_y, 1e3)
  f_col    <- Pi_col^2 / (1 + Pi_col^2)

  ## --- Module 7c: Per-mode porosity, sphericity, particle distribution -----
  f_trap_s <- theta_skin / (theta_skin + Perm_shell)
  dP_vap   <- P_atm * pmax(exp(0.025 * (T_particle - T_bpS)) - 1, 0)
  Pi_b     <- dP_vap * f_trap_s * S_bs / max(sigma_y, 1e3)
  f_burst  <- Pi_b^2 / (1 + Pi_b^2)
  B_infl   <- 1.5 * f_trap_s * (1 - f_burst)
  phi_templ <- min(0.6, phi_e * S_bs * (1 + B_infl) * (1 - 0.5 * f_burst))
  D_pore    <- D_e_h * (1 + B_infl)^(1/3)

  a_trap_j <- alpha_e / (1 + D_b_e / modes_m) * (0.5 + 0.5 * theta_surf)
  phi_str_eff <- phi_struct * (1 - 0.7 * f_col)
  f_trap  <- theta_skin / (theta_skin + Perm_shell)
  S_boil  <- 1 / (1 + exp(-(T_particle - 373.15) / 8))
  phi_vac <- 0.30 * f_trap * S_boil
  phi_int <- 1 - (1 - phi_str_eff) * (1 - phi_vac) *
                 (1 - phi_templ)
  phi_j <- a_trap_j + (1 - a_trap_j) * phi_int
  phi_j <- phi_j * (1 - 0.6 * S_stick * (1 - 0.5 * theta_skin))
  phi_porosity <- sum(modes_w * phi_j)

  resist <- (sigma_eff / 0.05) * (mu_eff0 / 0.01)^0.3
  Omega_struct <- 1 - 0.55 * theta_skin * (0.3 + 0.7 * f_col) /
                    ((1 + resist) * (1 + 0.1 * Softness))
  Omega_struct <- Omega_struct + (1 - Omega_struct) * 0.4 * S_stick
  # Burst denting: core-absorbed template flashing inside plasticized cores
  # counts toward the volatile inventory alongside the free emulsion
  vol_load <- phi_e + w_core * rho_lm / rho_solv
  Omega_struct <- Omega_struct * (1 - 0.30 * phi_vac) *
                  (1 - 0.35 * f_burst * min(3 * vol_load, 1))

  Dp_j <- modes_m * ((1 - a_trap_j) * (1 - phi_e) * rho_L * C_sol /
                     (rho_s * (1 - phi_j)))^(1/3)
  grp  <- seq(log(0.05 * min(Dp_j)), log(10 * max(Dp_j)), length.out = 400)
  cdfp <- Reduce(`+`, Map(function(w, m, s)
            w * pnorm((grp - log(m)) / log(s)), modes_w, Dp_j, modes_s))
  qp_at <- function(p) exp(approx(cdfp, grp, xout = p, ties = "ordered")$y)
  Dp50 <- qp_at(0.50); Dp90 <- qp_at(0.90)

  ## --- Module 7d: Tapped density -------------------------------------------
  rho_env <- sum(modes_w * ((1 - phi_j) * rho_s + phi_j * 1.2))
  f_pack  <- 0.64 * (0.55 + 0.45 * Omega_struct) *
             (Dp50 / (Dp50 + 8e-6)) *
             (1 + 0.5 * modes_w[3] * modes_w[1]) *
             (1 - 0.2 * S_stick) *
             (1 - 0.15 * min(Softness / 25, 1)) *
             (1 - 2 * min(X_moist, 0.2))
  rho_tapped <- rho_env * max(f_pack, 0.05)

  c(d_droplet_um   = d50 * 1e6,
    d10_um         = d10 * 1e6,
    d90_um         = d90 * 1e6,
    d99_um         = d99 * 1e6,
    span           = span,
    BI_bimodal     = BI,
    D_particle_um  = Dp50 * 1e6,
    Dp90_um        = Dp90 * 1e6,
    theta_skin_z   = theta_skin,
    Omega_struct_z = Omega_struct,
    phi_porosity_z = phi_porosity,
    rho_tapped     = rho_tapped,
    Tg_product_K   = Tg_eff,
    X_moisture     = X_moist,
    Perm_shell_rel = Perm_shell,
    D_pore_um      = D_pore * 1e6,
    solv_retained  = R_solv,
    f_burst_solv   = f_burst,
    sigma_y_cake_MPa = sigma_y / 1e6)
}
