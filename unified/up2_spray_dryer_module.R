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
  R_GAS    = 8.314,     # universal gas constant                  [J/mol K]
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
                      "f_cr_solv", "f_esc_solv",
                      "sigma_y_cake_MPa",
                      "alpha_g_nucl", "D_b_nucl_um")

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
  k_emu0 <- cst$k_emu0; R_GAS <- cst$R_GAS

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
  # Override: use UP1 mixer exit slurry viscosity directly when available,
  # bypassing the serum-based KD build-up (avoids double-counting floc effects).
  if (!is.null(feed$mu_exit_PaS) && !is.na(feed$mu_exit_PaS)) {
    mu_slurry0 <- feed$mu_exit_PaS
  }
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
  # Runs on the feed-state theta_surf (pre-nucleation): DLVO stability is a
  # property of the slurry entering the dryer, before nozzle pressure drop.
  E_rep     <- dpH * (1 / (1 + 10 * I_str)) * (1 - theta_surf)
  stability <- 1 + E_rep + 5 * theta_surf         # electrostatic + steric

  ## --- Module 0f: Dissolved-gas nucleation at the atomizer -----------------
  # Dissolved gas stays dissolved in the pressurised transfer line (no gas/water
  # interface → zero Ostwald ripening, unlike pre-sparged free bubbles). It
  # nucleates heterogeneously at the polymer-water interface when feed pressure
  # drops at the nozzle. UP1 aggregates (D_agg ~10 µm) are already formed at
  # this point, so bubbles nucleate INTO aggregate pore space — already in
  # position for pore retention. Local skin softening only (surface layer):
  # gaseous alkyl monomer (bp < RT) partitions to the particle surface but
  # cannot accumulate in glassy bulk cores → does NOT raise w_core or Softness.
  #
  # NOTE on C_gas_diss units: the UP1 ODE tracks dissolved gas as a Henry-
  # normalised model state (C_eq = P_local/(5*T_scale) ~ 0.2 at 1 atm), not a
  # strict kg/kg mass fraction. K_alpha_gas converts from ODE model units to gas
  # volume fraction; calibrate experimentally for your gas/pressure system
  # (default 1.0 preserves the original phenomenological calibration).
  r_exp_0f <- max(P_F / P_atm, 1)
  C_gas_diss <- if (!is.null(feed$C_gas_diss) && is.finite(feed$C_gas_diss))
                  pmax(0, feed$C_gas_diss) else 0.0

  alpha_g_nucl <- 0.0
  D_b_nucl     <- D_b_h   # default: no nucleation event, carry existing size

  if (C_gas_diss > 1e-4) {
    # Henry-scaled release: K_alpha_gas maps ODE concentration units to gas
    # volume fraction; (1 - 1/r_exp) is the depressurisation release fraction.
    K_alpha_gas  <- if (!is.null(x[["K_alpha_gas"]])) x[["K_alpha_gas"]] else 1.0
    alpha_g_nucl <- min(C_gas_diss * K_alpha_gas * (1 - 1 / r_exp_0f), 0.90)

    # Critical-nucleus diameter from Laplace balance at heterogeneous cavity.
    dP_flash      <- max(P_F - P_atm, 1e3)                       # [Pa]
    D_b_nucl_crit <- min(max(4 * sigma / dP_flash, 1e-7), 20e-6) # 0.1–20 µm

    # Post-nucleation diffusive growth: redistribute alpha_g_nucl across
    # N_bub nucleation sites per m³ of slurry (default 1e15 m⁻³ gives 1–5 µm
    # final bubble diameters at typical alpha_g_nucl 0.01–0.10).
    # The user may calibrate N_bub_m3 from observed pore-size distributions.
    N_bub_m3  <- if (!is.null(x[["N_bub_m3"]])) x[["N_bub_m3"]] else 1e15
    D_b_grown <- (6 * max(alpha_g_nucl, 1e-12) / (pi * N_bub_m3))^(1 / 3)
    D_b_nucl  <- max(D_b_nucl_crit, D_b_grown)  # grown >= critical nucleus

    # Volume-weighted merge: nucleated bubbles join any pre-existing free gas.
    a_pre  <- alpha_g
    a_tot  <- min(a_pre + alpha_g_nucl, 0.90)
    D_b_h  <- if (a_tot > 1e-8)
                ((a_pre * D_b_h^3 + alpha_g_nucl * D_b_nucl^3) / a_tot)^(1/3)
              else D_b_h
    alpha_g <- a_tot

    # Surfactant debit: nucleation adds new gas/liquid interface; recompute
    # coverage with updated total gas area and shed under-stabilised foam.
    gas_area_post <- 6 * alpha_g / max(D_b_h, 1e-7)
    theta_surf    <- min(1, cap_area / (part_area + gas_area_post + emu_area))
    alpha_g       <- alpha_g * (0.3 + 0.7 * theta_surf)

    # Local skin softening: dissolved gas partitions to particle surface layer
    # (~10–50 nm) only — tracked via theta_seed, not Softness or w_core.
    k_part_gas    <- if (!is.null(x[["k_part_gas"]])) x[["k_part_gas"]] else 0.05
    phi_skin_gas  <- min(C_gas_diss * k_part_gas, 0.10)
    skin_gas_frac <- phi_skin_gas / (phi_skin_gas + 0.03)
    theta_seed    <- 1 - (1 - theta_seed) * (1 - 0.4 * skin_gas_frac)
  }

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

  # Péclet concentration route: softer (more plasticized) surface has lower
  # viscous resistance to consolidation -> Softness LOWERS the threshold.
  # Opposite sign from the old "shell stiffness inhibits" rationale, which
  # conflated the falling-rate shell with the forming surface layer.
  S_crit_pe     <- 500 / (1 + 0.05 * Softness)
  theta_skin_pe <- S_skin / (S_skin + S_crit_pe)

  # Surface fusion route: the constant-rate wet-bulb surface pins at ~100°C
  # while free water is present. Plasticizer/monomer lower the polymer Tg (Fox
  # equation); when Tg_plas drops below 100°C the surface is rubbery throughout
  # the constant-rate period -> continuous coalescence. Strongest differentiator
  # for high-Tg polymers where plasticizer determines whether fusion fires at all.
  # T_particle (the falling-rate temperature) is retained for burst/boil/stick.
  T_surface_cr  <- 373.15   # constant-rate wet-bulb surface temperature [K]
  inv_Tg_plas   <- (1 - phi_solvent) / max(Tg_pol, 200) + phi_solvent / Tg_solv
  Tg_plas       <- 1 / inv_Tg_plas
  theta_skin_fus <- 1 / (1 + exp(-(T_surface_cr - Tg_plas) / 15))

  # Parallel paths: series non-skin probabilities (either route independently fires)
  # Weight 0.5 on fusion route avoids trivial saturation at nominal conditions
  theta_skin <- 1 - (1 - theta_skin_pe) * (1 - 0.5 * theta_skin_fus)

  # Cross-module coupling: mixer wet-skin seeds the dryer skin state
  theta_skin <- 1 - (1 - theta_seed) * (1 - theta_skin)

  D_f <- 1.8 + 0.7 * min((stability - 1) / 5, 1)
  phi_struct <- 0.30 * theta_skin * (1 - d_ratio^(3 - D_f))

  ## --- Module 6c: Per-mode drying & residual moisture ----------------------
  tau_dry <- (modes_m^2 / kappa) * (1 + 20 * theta_skin / Perm_shell)
  X_j     <- exp(-t_res / tau_dry)
  X_moist <- sum(modes_w * X_j)
  # Operational moisture reconciliation: the dryer air is energy-sized to leave
  # <= w_moist_target (0.5%) moisture in the powder, and the product is not
  # hygroscopic, so the final moisture cannot exceed the target regardless of the
  # per-mode drying kinetics. Cap X_moist (fraction of FEED water retained) at the
  # value that yields the target POWDER moisture w_wat at this solids loading:
  #   w_wat = X_moist(1-Csol)/(X_moist(1-Csol)+Csol)  ->  X_cap = t*Csol/((1-Csol)(1-t)).
  # Default 0.005; a future study can sweep w_moist_target up to ~0.01 (1%).
  w_moist_target <- if (!is.null(x[["w_moist_target"]])) x[["w_moist_target"]] else 0.005
  X_moist_cap <- w_moist_target * C_sol / (max(1 - C_sol, 1e-6) * (1 - w_moist_target))
  X_moist <- min(X_moist, X_moist_cap)

  ## --- Module 6c (Raoult): constant-rate partial-pressure evaporation ------
  # During the constant-rate drying period (near dryer inlet) template droplets
  # evaporate by partial-pressure driving force even when T_droplet < T_bp.
  # Droplet surface temperature during constant-rate is approximated as:
  #   T_wb_cr = T_out - 0.8*(T_in - T_out)   [co-current: hot/dry near inlet]
  # Vapour pressure via Clausius-Clapeyron with Trouton's rule (dH_vap = 88*T_bp):
  #   p_cr = P_atm * exp( (dH_vap/R) * (1/T_bp - 1/T_wb_cr) )
  # Escape fraction during constant-rate window (exponential approach to saturation):
  #   f_cr = 1 - exp( -C_cr * (p_cr/P_atm) * Perm_shell )
  # C_cr = 5.0 sets the scale so that BA (p_cr/P_atm~0.075 at T_wb~321K) gives
  # ~35-40% escape and BB (~0.020 at T_wb~321K) gives ~12%, consistent with the
  # 3.6x vapour-pressure ratio between BA and BB at wet-bulb conditions.
  # Constant-rate escape is driven by the SURFACE partial pressure of template
  # solvent minus the partial pressure already carried by the drying gas:
  #     dp = a1 * p_sat_solv(T_wb) - p_solv_gas
  # a1 is the Flory-Huggins solvent activity, a1 = phi*exp((1-phi)+chi*(1-phi)^2):
  #   * free droplets (phase-separated poor solvent, phi ~ 1) -> a1 ~ 1: escape
  #     as pure solvent and template the pores;
  #   * core-absorbed solvent (dissolved at low phi) -> a1 < 1: "held" by the
  #     polymer, the more strongly the LOWER chi (good solvent). Evaluated as the
  #     core dries (phi falls from phi_c0 toward zero) so the chi effect grows
  #     through the constant-rate period (composition-evolving activity).
  T_wb_cr <- T_out - 0.8 * pmax(T_in - T_out, 0)
  dHvap_J <- 88.0 * T_bpS                          # Trouton estimate [J/mol]
  p_sat_s <- pmin(P_atm, P_atm * exp(dHvap_J / R_GAS * (1/T_bpS - 1/T_wb_cr)))

  # Gas-phase solvent partial pressure: evaporated template accumulates in the
  # drying gas (solvent humidity Y_solv), but the vapour space is already partly
  # occupied by steam (water vapour p_v). The solvent partial pressure cannot
  # exceed the remaining headroom OR its own saturation -> a near-water-saturated
  # gas throttles template escape. This is the "how much can actually evaporate"
  # limit set by the drying-gas partial pressures.
  MW_air  <- 0.02896                               # dry air molar mass [kg/mol]
  MW_solv <- if (!is.null(x[["MW_solv"]])) x[["MW_solv"]] else 0.100   # ~100 g/mol
  Y_solv  <- (mdot_L * (w_e + w_core)) / max(mdot_g, 1e-9)  # kg solvent / kg dry gas (cap)
  p_solv0 <- Y_solv / (Y_solv + MW_air / MW_solv) * P_atm
  p_head  <- max(0, P_atm - p_v)                   # headroom left after steam (p_v)
  p_solv_gas <- min(p_solv0, p_head, p_sat_s)

  # chi is temperature-dependent (enthalpic part ~ 1/T): evaluate the input
  # chi_template (quoted at T_CHI_REF) at the constant-rate wet-bulb surface,
  # which is hotter than the mixer -> a poor solvent becomes somewhat less poor
  # as it dries. Entropic part chi_S ~ 0.34 is ~T-independent.
  chi_ref <- if (!is.null(feed$chi_template) && is.finite(feed$chi_template))
               feed$chi_template else 0.5
  CHI_S <- 0.34; T_CHI_REF <- 298.15
  chi_t <- CHI_S + (chi_ref - CHI_S) * (T_CHI_REF / max(T_wb_cr, 1.0))

  ## --- Module 7a: Product glass transition (Fox: solvent + moisture) -------
  T_particle <- 0.85 * T_out + 0.15 * T_feed
  S_bs   <- 1 / (1 + exp(-(T_particle - T_bpS) / 8))   # falling-rate boiling gate

  # Free-droplet template: evaporates as ~pure solvent (a1 ~ 1)
  dp_free    <- max(0, p_sat_s - p_solv_gas)
  f_cr_free  <- if (phi_e > 1e-4)
                  max(0, min(1, 1 - exp(-5.0 * (dp_free / P_atm) * Perm_shell))) else 0.0
  f_esc_free <- f_cr_free + (1 - f_cr_free) * S_bs

  # Core-absorbed template: FH-activity-suppressed, phi-evolving (activity
  # integrated as the core solvent fraction falls from phi_c0 toward dry)
  phi_c0     <- min(0.9, max(w_core, 1e-4))
  phis       <- seq(phi_c0, phi_c0 / 8, length.out = 8)
  a1s        <- phis * exp((1 - phis) + chi_t * (1 - phis)^2)
  dp_core    <- max(0, mean(pmax(0, a1s * p_sat_s - p_solv_gas)))
  f_cr_core  <- if (w_core > 1e-5)
                  max(0, min(1, 1 - exp(-5.0 * (dp_core / P_atm) * Perm_shell))) else 0.0
  f_esc_core <- f_cr_core + (1 - f_cr_core) * S_bs

  # Mass-weighted reported escape over the two solvent pools (free + core)
  m_free <- w_e; m_core <- w_core; m_tot <- m_free + m_core
  f_cr   <- if (m_tot > 1e-9) (m_free * f_cr_free  + m_core * f_cr_core ) / m_tot else f_cr_free
  f_esc  <- if (m_tot > 1e-9) (m_free * f_esc_free + m_core * f_esc_core) / m_tot else f_esc_free
  R_solv <- 1 - f_esc

  # Residual solvent load: each pool retained at its OWN escape fraction (free
  # template escapes as pure solvent -> pores; core template is held by FH
  # activity -> residual / collapse pathway)
  R_solv_free <- 1 - f_esc_free
  R_solv_core <- 1 - f_esc_core
  w_res  <- (0.5 * C_mono + C_plas + w_e * R_solv_free + w_core * R_solv_core) /
            (0.5 * C_mono + C_plas + w_e * R_solv_free + w_core * R_solv_core + C_sol)
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
  # Burst only from falling-rate trapped vapour (S_bs gate); constant-rate
  # evaporation (f_cr) exits freely before the shell seals.
  Pi_b     <- dP_vap * f_trap_s * S_bs / max(sigma_y, 1e3)
  f_burst  <- Pi_b^2 / (1 + Pi_b^2)
  # Shell inflation scaled by (1 - f_cr_free): free template that escaped in the
  # constant-rate period never reaches the sealed-shell stage and does not inflate.
  B_infl   <- 1.5 * f_trap_s * (1 - f_burst) * (1 - f_cr_free)
  # phi_templ: escaped FREE template (f_esc_free) creates pores; post-escape
  # inflation and burst logic as before. Core-absorbed solvent is tracked
  # separately (w_res / collapse), not as a pore former.
  phi_templ <- min(0.6, phi_e * f_esc_free * (1 + B_infl) * (1 - 0.5 * f_burst))
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

  # --- UP1 aggregate size-template regime (MUTED by default) ----------------
  # Two dry-PSD control regimes are seen in the calibration data (visc.xlsx):
  #   * atomizer-control  - dispersed feed (no UP1 aggregate): the droplet-shell
  #     above sets the product size (e.g. cond2, dry d50 19.85 um).
  #   * UP1-control       - a UP1 aggregate templates the particle, so the dry
  #     size tracks the aggregate: dry d50 ~ TEMPLATE_DENSIFY * D_agg (~1.2x wet).
  # This overlay lets the aggregate TEMPLATE the size, preserved as a capability
  # for a future experiment but muted for now via the `size_template` knob [0,1]:
  #   0 -> droplet-shell only (current calibration, bit-for-bit unchanged)
  #   1 -> full template in the aggregated regime
  # The atomizer-control regime is unaffected at ANY setting: d_ratio -> 1 for a
  # dispersed feed, so the template weight w_tmpl -> 0 there automatically.
  size_template <- if (!is.null(x[["size_template"]])) x[["size_template"]] else 0
  if (size_template > 0 && !is.null(feed$D_agg_um) &&
      is.finite(feed$D_agg_um) && feed$D_agg_um > 0) {
    TEMPLATE_DENSIFY <- 1.20                              # dry/wet, UP1-control set
    # Aggregation fraction: how far D_agg sits above the un-aggregated baseline
    # D_primary_exit_um (what D_agg floors to when the feed never aggregates).
    # -> 0 in the atomizer-control regime (dispersed feed, D_agg == baseline)
    #    so the droplet-shell is kept untouched there at any knob setting.
    D_pexit  <- if (!is.null(feed$D_primary_exit_um) && is.finite(feed$D_primary_exit_um))
                  feed$D_primary_exit_um else feed$D_agg_um
    agg_frac <- max(0, 1 - min(D_pexit / feed$D_agg_um, 1))
    w_tmpl    <- agg_frac * min(max(size_template, 0), 1)
    Dp50_sh   <- exp(sum(modes_w * log(Dp_j)))           # shell distribution centre
    Dp_target <- TEMPLATE_DENSIFY * feed$D_agg_um * 1e-6
    Dp_j      <- Dp_j * (Dp_target / Dp50_sh)^w_tmpl     # log-blend toward template
  }

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

  ## --- Particle-morphology recalibration (MUTED by default) -----------------
  # Anchored to the measured powder bulk density (~0.30 g/cc) with the true
  # polymer skeletal density 1.70 g/cc. The current closure returns near-solid
  # particles (phi ~ 0.03) and reaches ~0.22 g/cc tapped only via an implausibly
  # low packing fraction -- right bulk number, wrong mechanism. This overlay
  # replaces porosity, sphericity and tapped density with a physically anchored
  # pair (bulk = rho_skel * (1 - phi) * packing):
  #   * compact UP1-control granule: intra-particle void phi ~ 0.65 (loosely
  #     packed primary aggregate), packing 0.50 -> 1.70*(1-0.65)*0.50 = 0.30 g/cc.
  #   * atomizer-control (dispersed) feed: central-void hollow bump (SEM cond2:
  #     crumpled collapsed shells) -> total void ~0.85, tapped ~0.13 g/cc, and
  #     lower sphericity (irregular vs rounded granule).
  # Regime weight from the same aggregation fraction the size-template uses
  # (D_primary_exit vs D_agg): dispersed feed -> w_disp -> 1, aggregated -> 0.
  # Governed by `morphology_recal` in [0,1] (0 = current closure, unchanged).
  # NOTE: uses skeletal 1.70 g/cc for the DRY powder; the UP1 slurry still uses
  # rho_polymer = 1050. If 1.70 is the true polymer density, UP1 phi_s (hence the
  # viscosity/aggregation-onset calibration) should be revisited -- see DATA_REVIEW.
  morphology_recal <- if (!is.null(x[["morphology_recal"]])) x[["morphology_recal"]] else 0
  if (morphology_recal > 0 && !is.null(feed$D_agg_um) &&
      is.finite(feed$D_agg_um) && feed$D_agg_um > 0) {
    RHO_SKEL <- 1700; F_PACK <- 0.50                 # kg/m3 ; inter-particle packing
    PHI_INTRA_BASE <- 0.65; HOLLOW_MAX <- 0.575      # compact void ; dispersed central-void bump
    OMEGA_COMPACT <- 0.90; OMEGA_HOLLOW_DROP <- 0.35 # rounded granule ; collapse penalty
    AGG_REF <- 0.70
    D_pexit_m  <- if (!is.null(feed$D_primary_exit_um) && is.finite(feed$D_primary_exit_um))
                    feed$D_primary_exit_um else feed$D_agg_um
    agg_frac_m <- max(0, 1 - min(D_pexit_m / feed$D_agg_um, 1))
    w_disp     <- max(0, 1 - min(agg_frac_m / AGG_REF, 1))     # 1 dispersed, 0 aggregated
    phi_hollow <- w_disp * HOLLOW_MAX
    phi_recal  <- 1 - (1 - PHI_INTRA_BASE) * (1 - phi_hollow)  # total void fraction
    omega_recal<- OMEGA_COMPACT - w_disp * OMEGA_HOLLOW_DROP
    tap_recal  <- RHO_SKEL * (1 - phi_recal) * F_PACK
    k <- min(max(morphology_recal, 0), 1)
    phi_porosity <- (1 - k) * phi_porosity + k * phi_recal
    Omega_struct <- (1 - k) * Omega_struct + k * omega_recal
    rho_tapped   <- (1 - k) * rho_tapped   + k * tap_recal
  }

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
    f_cr_solv      = f_cr,
    f_esc_solv     = f_esc,
    sigma_y_cake_MPa = sigma_y / 1e6,
    alpha_g_nucl   = alpha_g_nucl,
    D_b_nucl_um    = D_b_nucl * 1e6)
}
