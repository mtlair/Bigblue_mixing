#!/usr/bin/env Rscript
# =============================================================================
# Morris (Elementary Effects) sensitivity analysis of the reduced-order
# two-fluid (internal-mix) atomization + drying model described in
# "deepresearchreport.md", extended with closures from the Latex Coagulation
# Engine Integrated Master Specification (v30.0.0 / R-Engine v47):
#
#   * Flory-Huggins free-volume swelling: residual solvent (C_monomer +
#     C_plasticizer) softens the matrix, Softness = (1+5 C_binder) exp(25 phi_solvent)
#   * Product glass transition: Fox equation on residual solvent -> Tg_eff,
#     driving stickiness / pore collapse / caking above Tg
#   * Surfactant molar stoichiometry: theta_surf from C_surfactant, MW and
#     A_molecule vs colloid + bubble surface area
#   * DLVO electrostatics: E_repulsion = dpH * Debye screening * (1-theta_surf),
#     gating aggregation (skin onset) and fractal openness (D_f -> porosity)
#   * Krieger-Dougherty crowding viscosity from feed solids (phi_s)
#
# Spray-specific extensions beyond the v47 reactor spec:
#   * Pressurized-hold bubble coarsening: Ostwald ripening + coalescence of
#     the entrained foam during the pressurized residence time t_hold before
#     the nozzle (LSW r^3 ~ t kinetics, Henry's-law accelerated by pressure,
#     damped by surfactant film coverage)
#   * Effervescent-airblast hybrid exit plane: nothing atomizes in the feed
#     line (confined flow); the entrained gas acts AT the bi-fluid exit,
#     where its throat voidage thins the annular liquid film (reducing the
#     characteristic length fed to the airblast correlation) and bubbles
#     small enough to ride inside the films flash-shatter the fragments on
#     the final letdown (multiplicative, gated by bubble size vs film
#     thickness). The nozzle carries its own pressure ladder: liquid line
#     P_feed -> mixing chamber P_chamber = 0.8 P_system -> near-choked
#     throat -> ambient
#   * Dryer energy/moisture balance: dryer gas flow rate mdot_gas_dry and
#     inlet absolute humidity Y_in give co-current outlet temperature and
#     relative humidity, which drive the drying Peclet number and the
#     sticky-point state (replaces the former independent RH_gas knob)
#
# NOTE: the v47 impeller dissipation closure (v_tip) was removed - it belongs
# to an upstream unit operation; its history is carried by the feed bubble
# size D_b and the hold-time coarsening model.
#
# Module chain: Feed Properties -> Pressurized Hold -> Two-Phase
# Conditioning -> Effervescent Stage -> Bi-Fluid Airblast Stage ->
# Secondary Breakup -> Dryer Balance -> Particle Formation.
#
# Outputs screened (nomenclature-sheet symbols where they exist):
#   1. D_particle     - final (dry) particle size                    [um]
#   2. theta_skin_z   - skin network fraction (skin formation)       [-]
#   3. Omega_struct_z - structural memory / sphericity state         [-]
#   4. phi_porosity_z - total porosity (void fraction state)         [-]
#   5. rho_tapped     - powder tapped density (bulk analogue of
#                       rho_colloid_out / SG_out)                    [kg/m3]
#   6. Tg_product     - effective product glass transition T_g,eff   [K]
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
# Symbols follow the master nomenclature sheet where available:
#   C_solid_mass, C_monomer, C_plasticizer, C_binder, C_surfactant : feed
#     composition (wt/wt)
#   alpha_g_0  : total gas holdup entrained in the feed foam        [-]
#   P_system   : atomizing air / system supply pressure             [Pa]
#   T_system   : dryer gas temperature                              [K]
#   Delta_pH   : system delta pH vs isoelectric point               [pH units]
#   I_strength : ionic strength (Debye screening driver)            [M]
# Others use the spray report's notation (sigma, mu_L, ALR, D_b, ...).
# t_hold, mdot_gas_dry and Y_in are spray-line additions (pressurized
# residence before the nozzle; dryer gas flow and inlet humidity).
factors <- data.frame(
  name = c("ALR",           # air-liquid mass ratio  m_G/m_L           [-]
           "P_system",      # atomizing air supply pressure            [Pa]
           "P_feed",        # liquid feed line (hold) pressure         [Pa]
           "mdot_L",        # liquid feed mass flow                    [kg/s]
           "sigma",         # liquid surface tension                   [N/m]
           "mu_L",          # serum (continuous phase) viscosity       [Pa s]
           "rho_L",         # liquid (slurry) density                  [kg/m3]
           "alpha_g_0",     # feed foam quality / gas entrainment      [-]
           "D_b",           # feed bubble diameter (before shear)      [m]
           "C_solid_mass",  # solid mass fraction in feed              [-]
           "T_system",      # dryer gas INLET temperature              [K]
           "mdot_gas_dry",  # dryer gas mass flow rate                 [kg/s]
           "Y_in",          # dryer gas inlet absolute humidity        [kg/kg]
           "T_feed",        # feed / atomizing air temperature         [K]
           "t_hold",        # pressurized hold time before nozzle      [s]
           "C_monomer",     # residual monomer concentration           [-]
           "C_plasticizer", # plasticizer concentration                [-]
           "C_binder",      # binder concentration                     [-]
           "C_surfactant",  # formulated surfactant concentration      [-]
           "Delta_pH",      # delta pH vs isoelectric point            [-]
           "I_strength",    # ionic strength                           [M]
           "Tg_polymer"),   # dry-polymer glass transition             [K]
  min  = c( 1.0, 2.0e5, 1.5e5, 0.002, 0.030, 0.0012, 1000, 0.05, 2.0e-5, 0.05,
            330, 0.10, 0.001, 280,   5, 0.000, 0.000, 0.000, 1.0e-4, 0.2,
            1.0e-3, 280),
  max  = c(10.0, 7.0e5, 1.0e6, 0.020, 0.070, 0.0560, 1300, 0.60, 2.0e-4, 0.40,
            470, 1.00, 0.020, 330, 600, 0.020, 0.050, 0.050, 2.0e-2, 4.0,
            5.0e-1, 380),
  # sample wide-ranging positive factors log-uniformly
  log  = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
           FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE,
           TRUE, FALSE),
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
# 2. Built-in reduced-order model (one row of physical inputs -> 6 outputs)
# -----------------------------------------------------------------------------
# Fixed constants / nozzle geometry (report worked example + v47 defaults)
P_atm    <- 1.013e5   # ambient pressure                        [Pa]
R_air    <- 287       # gas constant, air                       [J/kg K]
gamma_a  <- 1.4       # heat capacity ratio, air                [-]
A_L      <- 1.0e-6    # liquid passage area                     [m2]
A_G      <- 5.0e-6    # gas passage area                        [m2]
D_h      <- 1.0e-3    # nozzle hydraulic diameter               [m]
rho_s    <- 1400      # dry solid (polymer) density             [kg/m3]
a_prim   <- 1.0e-7    # primary colloid particle radius         [m]
kB       <- 1.380649e-23
N_Av     <- 6.02214e23
MW_surf  <- 0.400     # surfactant molecular weight             [kg/mol]
A_molec  <- 0.5e-18   # surfactant molecule area capacity       [m2]
HLB      <- 12        # surfactant HLB (foam stability driver)  [-]
phi_m    <- 0.63      # Krieger-Dougherty max packing           [-]
Tg_solv  <- 150       # solvent (monomer/plasticizer) Tg, Fox   [K]
d_ratio  <- 0.10      # primary particle / aggregate size ratio [-]
C_cham   <- 0.80      # mixing chamber / air supply pressure    [-]
k_rip0   <- 4.0e-15   # Ostwald ripening rate at 1 atm          [m3/s]
cp_gas   <- 1005      # dryer gas heat capacity                 [J/kg K]
h_fg     <- 2.30e6    # latent heat of water evaporation        [J/kg]

spray_dry_model <- function(x) {
  ALR    <- x[["ALR"]];        P_G    <- x[["P_system"]]
  P_F    <- x[["P_feed"]]
  mdot_L <- x[["mdot_L"]];     sigma  <- x[["sigma"]]
  mu_L   <- x[["mu_L"]];       rho_L  <- x[["rho_L"]]
  alpha0 <- x[["alpha_g_0"]];  D_b    <- x[["D_b"]]
  C_sol  <- x[["C_solid_mass"]]
  T_in   <- x[["T_system"]];   mdot_g <- x[["mdot_gas_dry"]]
  Y_in   <- x[["Y_in"]];       T_feed <- x[["T_feed"]]
  t_hold <- x[["t_hold"]]
  C_mono <- x[["C_monomer"]];  C_plas <- x[["C_plasticizer"]]
  C_bind <- x[["C_binder"]];   C_surf <- x[["C_surfactant"]]
  dpH    <- x[["Delta_pH"]];   I_str  <- x[["I_strength"]]
  Tg_pol <- x[["Tg_polymer"]]

  ## --- Module 0a: Formulation - Flory-Huggins free-volume swelling ---------
  # (v47 Sect. 4) residual solvent acts as theta-solvent proxy, softening the
  # matrix exponentially and independently of temperature
  phi_solvent <- C_mono + C_plas
  Softness    <- (1 + 5 * C_bind) * exp(25 * phi_solvent)

  ## --- Module 0b: Feed rheology (Krieger-Dougherty crowding) ---------------
  phi_s     <- C_sol * rho_L / rho_s              # solid phase volume fraction
  mu_serum  <- mu_L * (1 + 10 * C_bind)           # binder thickening (K_fluid)
  mu_slurry <- mu_serum *
               (1 - min(phi_s, 0.60) / phi_m)^(-2.5 * phi_m)

  ## --- Module 0c: Surfactant molar stoichiometry (v47 Sect. 5) -------------
  # monolayer capacity vs total interfacial demand: primary colloid surface
  # (3*phi_s/a per m3, the dominant sink) plus bubble surface (6*alpha/D_b)
  cap_area   <- (C_surf / MW_surf) * rho_L * N_Av * A_molec   # [m2/m3]
  part_area  <- 3 * phi_s / a_prim                            # [m2/m3]
  gas_area   <- 6 * alpha0 / D_b                              # [m2/m3]
  theta_surf <- min(1, cap_area / (part_area + gas_area))
  Foam_Stab  <- 1 + 5 * theta_surf * (HLB / 10)
  # under-stabilized foam sheds entrained gas before the nozzle
  alpha_g    <- alpha0 * (0.3 + 0.7 * theta_surf)

  ## --- Module 0d: Pressurized hold - coalescence / Ostwald ripening --------
  # LSW kinetics D^3 ~ t; Henry's-law gas solubility makes ripening faster
  # under the liquid-line pressure, surfactant coverage retards both
  k_rip <- k_rip0 * (P_F / P_atm) * (1 - 0.8 * theta_surf)
  D_b_h <- (D_b^3 + k_rip * t_hold)^(1/3)         # coarsened bubble size

  ## --- Module 0e: DLVO electrostatics (v47 Sect. 6) ------------------------
  E_rep     <- dpH * (1 / (1 + 10 * I_str)) * (1 - theta_surf)
  stability <- 1 + E_rep + 5 * theta_surf         # electrostatic + steric

  ## --- Module 1: Feed properties (at liquid-line pressure P_feed) ----------
  rho_G0    <- P_F / (R_air * T_feed)             # entrained gas density
  rho_eff0  <- (1 - alpha_g) * rho_L + alpha_g * rho_G0   # rho_eff,z (HEM)
  mu_eff    <- mu_slurry * (1 + 2.5 * alpha_g)    # mu_eff,z (dilute-bubble)
  sigma_eff <- sigma * (1 - 0.5 * alpha_g)        # foam-reduced surface tension

  ## --- Module 2: Nozzle pressure ladder & two-phase conditioning -----------
  # liquid line P_feed -> mixing chamber -> near-choked throat -> ambient;
  # alpha_g_0 is measured at line pressure, so all expansions reference P_feed
  P_cham   <- C_cham * P_G                        # internal chamber pressure
  P_throat <- max(P_atm, 0.53 * P_cham)           # near-choked throat pressure
  r_th     <- max(P_F / P_throat, 1)              # feed -> throat expansion
  r_exp    <- max(P_F / P_atm, 1)                 # feed -> ambient expansion
  r_flash  <- max(P_throat / P_atm, 1)            # final letdown (shattering)
  alpha_th <- min(alpha_g * r_th / (1 - alpha_g + alpha_g * r_th), 0.97)
  alpha_e  <- min(alpha_g * r_exp / (1 - alpha_g + alpha_g * r_exp), 0.97)
  D_b_e    <- D_b_h * r_exp^(1/3)                 # bubble growth D_b ~ P^-1/3
  rho_Ge   <- P_atm / (R_air * T_feed)
  rho_eff  <- (1 - alpha_e) * rho_L + alpha_e * rho_Ge

  ## --- Module 3: Nozzle hydraulics ------------------------------------------
  U_L <- mdot_L / (rho_L * A_L)                   # liquid exit velocity
  # gas: fully-expanded isentropic velocity (choked upstream if P_G/P_atm>1.89)
  U_G <- 0.9 * sqrt(pmax(2 * gamma_a / (gamma_a - 1) * R_air * T_feed *
                         (1 - (P_atm / P_G)^((gamma_a - 1) / gamma_a)), 0))
  U_rel <- max(U_G - U_L, 5)
  rho_A <- P_atm / (R_air * T_feed)               # ambient air density

  ## --- Module 4: Effervescent-airblast hybrid exit plane -------------------
  # Confined flow cannot atomize: the entrained gas does its work AT the
  # bi-fluid exit orifice, simultaneously with the airblast shear.
  # (a) throat voidage of the entrained gas thins the annular liquid film,
  #     shrinking the characteristic length the airblast correlation sees
  t_film <- (D_h / 2) * (1 - sqrt(alpha_th))      # annular film thickness
  L_c    <- max(2 * t_film, 0.02 * D_h)           # characteristic length
  SMD_ab <- 0.48 * L_c * (sigma_eff / (rho_A * U_rel^2 * L_c))^0.4 *
              (1 + 1 / ALR)^0.4 +
            0.15 * L_c * sqrt(mu_eff^2 / (sigma_eff * rho_L * L_c)) *
              (1 + 1 / ALR)
  # (b) bubbles small enough to ride inside the films burst on the final
  #     letdown and subdivide the fragments (multiplicative refinement);
  #     coarse (post-hold) bubbles vent to the gas core instead
  D_b_th  <- D_b_h * r_th^(1/3)                   # bubble size at the throat
  g_in    <- 1 / (1 + D_b_th / t_film)            # fraction riding inside film
  F_flash <- (1 + alpha_th * (r_flash - 1))^(-1/3)
  SMD     <- SMD_ab * (1 - g_in * (1 - F_flash))

  ## --- Module 5: Secondary breakup (Weber-number correction) ---------------
  We_g <- rho_A * U_rel^2 * SMD / sigma_eff       # gas Weber number of drop
  f_sec <- 0.55 + 0.45 / (1 + We_g / 80)          # ~20-45 % diameter reduction
  d_drop <- SMD * f_sec                           # droplet Dv50 entering dryer

  ## --- Module 6a: Dryer energy / moisture balance (co-current) -------------
  mdot_w <- mdot_L * (1 - C_sol)                  # water evaporation load
  T_out  <- T_in - mdot_w * h_fg / (mdot_g * cp_gas)
  T_out  <- max(T_out, T_feed + 2)                # saturated / underpowered dryer
  Y_out  <- Y_in + mdot_w / mdot_g                # outlet absolute humidity
  p_v    <- Y_out / (Y_out + 0.622) * P_atm       # vapour partial pressure
  p_sat  <- 610.94 * exp(17.625 * (T_out - 273.15) / (T_out - 273.15 + 243.04))
  RH_out <- min(p_v / p_sat, 0.99)                # outlet relative humidity

  ## --- Module 6b: Drying kinetics -------------------------------------------
  # d2-law evaporation coefficient, scaled by outlet-state driving force
  kappa <- 5e-8 * pmax(T_out - T_feed, 1) / 100 * (1 - RH_out)  # [m2/s]
  # Stokes-Einstein diffusivity of primary colloid particles
  D_diff <- kB * T_feed / (6 * pi * mu_L * a_prim)              # [m2/s]
  Pe <- kappa / (8 * D_diff)                      # drying Peclet number

  # theta_skin,z : skin network fraction. Colloidal instability (low DLVO
  # stability) accelerates shell aggregation; Flory-Huggins softness delays
  # rigid lock-in (soft particles film-form instead of jamming)
  S_skin <- Pe * C_sol * (1 + 1 / stability)
  S_crit <- 500 * (1 + 0.05 * Softness)
  theta_skin <- S_skin / (S_skin + S_crit)

  # Fractal openness (v47 Sect. 9): unstable colloid -> DLCA, open flocs
  # (low D_f, high structural porosity); stable -> compact RLCA packing
  D_f <- 1.8 + 0.7 * min((stability - 1) / 5, 1)
  phi_struct <- 0.30 * theta_skin * (1 - d_ratio^(3 - D_f))

  # alpha_trap,z : gas retained inside droplets - bubbles comparable to or
  # larger than the droplet are shed during breakup; surfactant-stabilized
  # films (theta_surf) resist rupture and hold gas in
  alpha_trap <- alpha_e / (1 + D_b_e / d_drop) * (0.5 + 0.5 * theta_surf)

  # phi_porosity,z : total porosity (void fraction state)
  phi_porosity <- alpha_trap + (1 - alpha_trap) * phi_struct

  ## --- Module 7: Product glass transition & thermal state ------------------
  # Fox equation on residual solvent retained in the dry particle
  # (monomer partially evaporates; plasticizer stays)
  w_res  <- (0.5 * C_mono + C_plas) / (0.5 * C_mono + C_plas + C_sol)
  Tg_eff <- 1 / ((1 - w_res) / Tg_pol + w_res / Tg_solv)

  # stickiness: particle temperature vs Tg_eff + 20 K sticky-point offset;
  # in co-current drying particles approach the outlet gas temperature
  T_particle <- 0.85 * T_out + 0.15 * T_feed
  S_stick <- 1 / (1 + exp(-(T_particle - (Tg_eff + 20)) / 10))

  # above Tg the matrix flows: pores collapse (unless skin-locked)
  phi_porosity <- phi_porosity * (1 - 0.6 * S_stick * (1 - 0.5 * theta_skin))

  # Omega_struct,z : sphericity; high sigma & mu resist skin buckling, a
  # plasticized (soft) matrix re-rounds, and stickiness anneals dents
  resist <- (sigma_eff / 0.05) * (mu_eff / 0.01)^0.3
  Omega_struct <- 1 - 0.55 * theta_skin / ((1 + resist) * (1 + 0.1 * Softness))
  Omega_struct <- Omega_struct + (1 - Omega_struct) * 0.4 * S_stick

  # Final particle size: solids mass balance on one droplet
  D_particle <- d_drop * ((1 - alpha_trap) * rho_L * C_sol /
                          (rho_s * (1 - phi_porosity)))^(1/3)

  # Envelope density and tapped density (v47 Sect. 10 density closure);
  # sticky / plasticized powders cake and pack worse, fines cohere
  rho_env  <- (1 - phi_porosity) * rho_s + phi_porosity * 1.2
  f_pack   <- 0.64 * (0.55 + 0.45 * Omega_struct) *
              (D_particle / (D_particle + 8e-6)) *
              (1 - 0.2 * S_stick) *
              (1 - 0.15 * min(Softness / 25, 1))
  rho_tapped <- rho_env * f_pack

  c(d_droplet_um   = d_drop * 1e6,
    D_particle_um  = D_particle * 1e6,
    theta_skin_z   = theta_skin,
    Omega_struct_z = Omega_struct,
    phi_porosity_z = phi_porosity,
    rho_tapped     = rho_tapped,
    Tg_product_K   = Tg_eff)
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

titles <- c(d_droplet_um   = "Spray droplet size  Dv50 [um]",
            D_particle_um  = "Final particle size  D_particle [um]",
            theta_skin_z   = "Skin formation  theta_skin,z [-]",
            Omega_struct_z = "Sphericity  Omega_struct,z [-]",
            phi_porosity_z = "Porosity  phi_porosity,z [-]",
            rho_tapped     = "Tapped density  rho_tapped [kg/m3]",
            Tg_product_K   = "Product glass transition  Tg_eff [K]")

plot_morris_panel <- function(st, title, n_label = 10) {
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25),
       ylim = c(0, max(st$sigma) * 1.30),
       pch = 21, bg = "steelblue", cex = 1.3,
       xlab = expression(mu * "*  (mean |elementary effect|)"),
       ylab = expression(sigma * "  (sd of elementary effects)"),
       main = title, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey55")           # sigma = mu* (non-linear)
  abline(0, 0.1, lty = 3, col = "grey70")         # sigma = 0.1 mu* (~linear)
  # label only the strongest factors to keep panels readable
  top <- order(-st$mu.star)[seq_len(min(n_label, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], labels = st$factor[top], pos = 3,
       cex = 0.70, offset = 0.35, xpd = NA)
}

png(file.path("output", "morris_sensitivity_plots.png"),
    width = 3000, height = 1500, res = 160)
op <- par(mfrow = c(2, 4), mar = c(4.5, 4.5, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (o in outputs) plot_morris_panel(stats_list[[o]], titles[[o]])
mtext(sprintf(paste("Morris screening: %d trajectories, %d levels, %d model",
                    "runs | dashed: sigma = mu* (non-linear/interacting),",
                    "dotted: sigma = 0.1 mu* (near-linear) |",
                    "top 10 factors labelled per panel"),
              r_traj, levels, nrow(Y)),
      side = 1, outer = TRUE, cex = 0.75, line = 1)
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
