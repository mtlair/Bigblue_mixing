library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN
#   decant pool (hindered settling) + plug foam + bubble population + gas state
#   + surfactant-driven film elasticity/drainage + T/P thermodynamic state
# =====================================================================
# Preformed, particle-loaded foam (no sparger) enters the base. Two zones split
# at the fixed decant depth H_pool. Properties (mu, rho_foam, d_p, d_b_in,
# d_sep) are passed down from the upstream unit.
#
# NEW in this version:
#  (1) HINDERED SETTLING (Richardson-Zaki) + KRIEGER-DOUGHERTY crowding: the
#      pool Stokes velocity is slowed by holdup, U = U_stokes * (1 - phi)^n_RZ,
#      AND by the local suspension viscosity. As liquid drains, solids
#      concentrate in the Plateau borders (local fraction eps_s/(eps_s+eps_l)),
#      so the effective viscosity mu_eff = mu(T) * (1 - phi_s/phi_smax)^(-2.5 phi_smax)
#      rises and further hinders settling. High local border viscosity is a real
#      retention mechanism, not just holdup.
#  (2) BUBBLE POPULATION: mean bubble size d_b(z) evolves by COALESCENCE
#      (grows, faster in dry/unstable foam) minus BREAKAGE (shear, restores
#      toward the inlet size). Coalescence is now damped by the SURFACTANT-set
#      film elasticity (see steps 1-3 below), not a hand-tuned constant.
#  (3) GAS STATE (no gas generated -- no sparger, no vent). Coalescence and
#      bubble bursting DO happen in the column, but the gas is NOT released
#      here: it stays in the column and is only let out in the DOWNSTREAM
#      open-atmosphere solids-concentration stage ("next up"). So the total
#      gas is CONSERVED and carried out the top -- what changes with height is
#      its STATE, not its amount. Bursting/coalescence convert dispersed foam
#      gas (J_g_foam) into large retained SLUG gas (J_g_slug); the slug gas
#      cannot escape until the next unit. Bursting still coarsens bubbles and
#      still dumps coarse particles (collapse loss tied to the burst rate).
#      Total gas out the top J_g_foam + J_g_slug == gas in the bottom.
#
# SURFACTANT / FILM CHAIN (de-lumps the old scalar film_stability):
#  step 1  Langmuir adsorption of c_surf -> surface excess Gamma; Szyszkowski
#          surface tension sigma(Gamma,T); Gibbs-Marangoni elasticity E_gibbs,
#          rolled off above the CMC by micelle buffering.
#  step 2  film_stability := E_gibbs / E_stab_ref  (normalized ~1 at baseline),
#          so surfactant type/dose -> elasticity -> coalescence & regime.
#  step 3  DRAINAGE is driven by the physical timescale tau_drain(mu(T),
#          surface mobility) toward an equilibrium holdup set by film stability,
#          replacing the hand-set tau_drain / eps_l_dry (now baseline anchors
#          scaled by dimensionless factors that are 1.0 at the baseline state).
#  T/P thermodynamic state: mu_cont(T) (Andrade), sigma(T), rho_gas(P,T) ideal
#          gas. Pressure sets gas density (and, upstream, d_b_in).
#
# States up z: Js_fine, Js_mid, Js_crs [m/s]; C_imp [%]; eps_l [-];
#              d_b [m]; J_g_foam, J_g_slug [m/s]; t_res [s].
# PLACEHOLDER kinetics calibrated to: residence ~1.5-2 h, top solids 3-7%.
# =====================================================================

GRAV <- 9.81

params <- c(
  H_total = 5.0, H_pool = 0.45, U_up = 8.0e-4,

  # --- UPSTREAM-PASSED PROPERTIES ---
  mu_cont = 2.0e-3, rho_liquid = 1050, rho_foam = 575, rho_p = 2500, rho_gas = 1.8,
  d_fine = 1.5e-5, d_mid = 8.0e-5, d_crs = 3.0e-4,
  d_b_in = 2.0e-3, d_sep = 4.5e-5, n_RZ = 4.65, phi_smax = 0.64,

  # --- FEED ---
  eps_l_in = 0.50, eps_s_in = 0.020, eps_g_in = 0.48,
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020, C_imp_in = 100,

  # --- THERMODYNAMIC STATE (T/P) ---
  T_col = 298.15,         # column temperature [K]
  P_col = 1.5e5,          # column absolute pressure [Pa] (~1.5 bar; open-atmosphere is downstream)
  R_gas = 8.314,          # gas constant [J/mol/K]
  MW_gas = 0.029,         # gas molar mass [kg/mol]
  mu_ref = 2.0e-3,        # continuous-phase viscosity at T_ref [Pa s]
  T_ref = 298.15,         # reference temperature [K]
  E_visc = 1800,          # Andrade viscosity activation temp (B in exp(B/T)) [K]
  sigma_ref = 0.045,      # clean-solvent surface tension at T_ref [N/m]
  dsigma_dT = -1.5e-4,    # surface tension temperature coefficient [N/m/K]

  # --- SURFACTANT (step 1) ---
  c_surf = 5.0,           # bulk surfactant concentration [mol/m^3]
  Gamma_inf = 4.0e-6,     # Langmuir plateau surface excess [mol/m^2]
  K_ads = 2.0,            # Langmuir adsorption constant [m^3/mol]
  cmc = 6.0,              # critical micelle concentration [mol/m^3]
  K_mic = 1.0,            # micelle-buffering coeff above CMC [-]
  mu_surf = 1.0e-6,       # surface (dilational) viscosity [N s/m]
  mu_surf_ref = 1.0e-6,   # reference surface viscosity (mobility=1 at baseline) [N s/m]

  # --- FILM ELASTICITY -> STABILITY (step 2) ---
  E_stab_ref = 0.099,     # Gibbs elasticity that maps to film_stability = 1 [N/m]

  # --- FILM DRAINAGE (physical; step 3) ---
  eps_l_dry = 0.15,       # baseline equilibrium holdup anchor [-]
  tau_drain = 2000,       # baseline drainage timescale anchor [s]
  hold_exp = 0.6,         # sensitivity of equilibrium holdup to film stability [-]
  eps_g_pack = 0.74, blend_wz = 0.03,

  # --- FILM THICKNESS / PLATEAU-BORDER GEOMETRY (wired into drainage/holdup) ---
  A_ham = 1.0e-20,        # Hamaker constant [J]
  lambda_D = 1.0e-8,      # Debye length [m]
  Pi_charge = 1.0e4,      # electrostatic disjoining-pressure scale (x coverage) [Pa]
  c_pb = 0.20,            # Plateau-border radius coeff: r_pb = c_pb*d_b*sqrt(eps_l) [-]
  p_perm = 2.0,           # border-drainage permeability exponent (drainage ~ r_pb^2) [-]
  r_pb_ref = 2.69e-4,     # reference border radius (tau normalizer; baseline foam-zone mean) [m]
  h_eq_ref = 4.195e-8,    # reference equilibrium film thickness (holdup normalizer; baseline) [m]
  q_film = 0.5,           # equilibrium-holdup sensitivity to film thickness [-]

  # --- BUBBLE POPULATION (coalescence - breakage) ---
  K_coal = 1.3e-4,        # coalescence coeff [1/s] (grows d_b)
  K_break = 1.2e-4,       # breakage coeff [1/s] (restores toward d_b_in)
  # --- BURST TRIGGER (film rupture -> foam gas becomes slug) ---
  d_b_burst = 3.0e-3,     # baseline critical burst size (at baseline film state) [m]
  K_burst = 1.4e-3,       # burst-rate coeff [1/s] (foam->slug conversion, not a vent)
  a_fs = 1.0,             # critical-size sensitivity to film stability (elasticity) [-]
  a_sig = 1.0,            # critical-size sensitivity to surface tension [-]
  sigma_ref_film = 0.02122, # baseline film surface tension (normalizer) [N/m]
  n_visc = 1.0,           # burst sensitivity to viscosity (drainage speed) [-]
  mu_drain_ref = 2.0e-3,  # fixed reference viscosity for burst normalization [Pa s]
  k_armor = 2.0,          # solids-armoring stabilization (local content) [-]
  k_bridge = 0.5,         # coarse-particle film-bridging destabilization (local size) [-]
  K_bsink = 0.15,         # burst-driven mean-d_b reduction (removes coarsest bubbles) [-]
  d_b_cap = 2.0e-2,       # d_b cap used in border geometry (h_eq robustness) [m]
  h_eq_max = 2.0e-7,      # cap on reported equilibrium film thickness [m]

  # --- WASH ---
  wash_ratio = 0.5, k_drainage = 1.8e-3, channel_cv = 0.35, d_b_ref = 2.0e-3,

  # --- PARTICLE LOSS ---
  k_det_fine = 1.2e-5, k_det_mid = 6.0e-5, k_det_crs = 3.0e-4,   # detachment [1/s]
  k_burstloss_fine = 0.05, k_burstloss_mid = 0.30, k_burstloss_crs = 2.0, # burst dumps [-]

  # --- REGIME (film_stability is DERIVED from surfactant; see derive_state_props) ---
  film_stability = 1.0, load_sens = 1.0, plug_crit = 1.0
)

stokes  <- function(dRho, d, mu) dRho * GRAV * d^2 / (18 * mu)
rho_agg <- function(dp, db, rp, rg) (rp*dp^3 + rg*db^3) / (dp^3 + db^3)

# --- THERMODYNAMIC + SURFACTANT DERIVATION (steps 1-3, T/P) ----------------
# Turns raw inputs (T, P, surfactant, surface viscosity) into the effective
# closures the ODE uses. All dimensionless factors are 1.0 at the baseline
# state, so the baseline reproduces the previously calibrated column exactly.
derive_state_props <- function(p) {
  Tk <- p[["T_col"]]

  # T/P thermodynamic state
  mu_T      <- p[["mu_ref"]] * exp(p[["E_visc"]] * (1/Tk - 1/p[["T_ref"]]))   # Andrade
  rho_gas_T <- p[["P_col"]] * p[["MW_gas"]] / (p[["R_gas"]] * Tk)             # ideal gas
  sigma0_T  <- p[["sigma_ref"]] + p[["dsigma_dT"]] * (Tk - p[["T_ref"]])      # clean-solvent sigma(T)

  # step 1: Langmuir adsorption (monomer activity capped at the CMC)
  c_eff  <- min(p[["c_surf"]], p[["cmc"]])
  theta  <- p[["K_ads"]] * c_eff / (1 + p[["K_ads"]] * c_eff)
  theta  <- min(theta, 0.995)                              # keep sigma/E finite
  Gamma  <- p[["Gamma_inf"]] * theta
  sigma  <- sigma0_T + p[["R_gas"]] * Tk * p[["Gamma_inf"]] * log(1 - theta)  # Szyszkowski
  # Gibbs-Marangoni elasticity with micelle buffering above the CMC
  micelle_buffer <- 1 / (1 + p[["K_mic"]] * max(0, p[["c_surf"]] - p[["cmc"]]) / p[["cmc"]])
  E_gibbs <- p[["R_gas"]] * Tk * p[["Gamma_inf"]] * theta / (1 - theta) * micelle_buffer

  # step 2: elasticity -> film_stability (normalized ~1 at baseline)
  film_stab <- E_gibbs / p[["E_stab_ref"]]

  # step 3: drainage. Surface mobility (rigid/high surface viscosity -> slower).
  mobility      <- 0.5 * (1 + p[["mu_surf"]] / p[["mu_surf_ref"]])
  tau_drain_eff <- p[["tau_drain"]] * (mu_T / p[["mu_ref"]]) * mobility        # thicker/rigid -> slower

  # FILM THICKNESS -> EQUILIBRIUM HOLDUP: the foam stops draining where the
  # disjoining pressure the films can support balances the Plateau-border
  # capillary suction. Evaluate that balance at reference border geometry to
  # get the equilibrium film thickness h_eq_dry, and let a thicker equilibrium
  # film (stronger disjoining / lower sigma) hold a wetter foam. Normalized to
  # 1.0 at the baseline state, so the baseline equilibrium holdup is unchanged.
  r_pb_dry  <- p[["c_pb"]] * p[["d_b_in"]] * sqrt(p[["eps_l_dry"]])            # ref border radius
  P_c_dry   <- sigma / r_pb_dry                                               # capillary suction [Pa]
  Pi_e      <- p[["Pi_charge"]] * theta                                       # electrostatic scale [Pa]
  h_eq_dry  <- p[["lambda_D"]] * log(max(Pi_e / P_c_dry, 1 + 1e-9))           # eq film thickness [m]
  film_wet  <- h_eq_dry / p[["h_eq_ref"]]                                     # ~1 at baseline
  eps_l_dry_eff <- min(p[["eps_l_in"]],
                       p[["eps_l_dry"]] * film_stab^p[["hold_exp"]] * film_wet^p[["q_film"]])

  # overrides consumed by the ODE
  p[["mu_cont"]]        <- mu_T
  p[["rho_gas"]]        <- rho_gas_T
  p[["film_stability"]] <- film_stab
  p[["tau_drain"]]      <- tau_drain_eff
  p[["eps_l_dry"]]      <- eps_l_dry_eff

  # derived reporting quantities
  c(p, sigma_film = sigma, E_gibbs = E_gibbs, Gamma_surf = Gamma, theta_cov = theta,
       mu_cont_T = mu_T, mobility = mobility, h_eq_dry = h_eq_dry, film_wet = film_wet)
}

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    load_rel <- (Js_fine + Js_mid + Js_crs) / (Js_fine_in + Js_mid_in + Js_crs_in)
    eps_s    <- eps_s_in * load_rel
    eps_g    <- 1 - eps_l - eps_s
    U <- U_up
    dt_res <- 1 / U
    foam_frac <- 1 / (1 + exp(-(z - H_pool) / blend_wz))

    # (1) KRIEGER-DOUGHERTY local viscosity: solids concentrated in the drained
    # border liquid raise the effective viscosity that resists settling.
    phi_cond <- eps_s / (eps_s + eps_l)                              # solids fraction in border liquid
    mu_eff   <- mu_cont * (1 - min(phi_cond, 0.999*phi_smax)/phi_smax)^(-2.5*phi_smax)

    # FILM THICKNESS / PLATEAU-BORDER GEOMETRY (local). d_b is bounded (d_b_cap)
    # in the geometry so a runaway mean diameter cannot make the border curvature
    # (and h_eq) blow up -- the burst d_b-sink below keeps d_b physical anyway.
    d_b_eff <- min(d_b, d_b_cap)
    r_pb   <- c_pb * d_b_eff * sqrt(max(eps_l, 1e-3))                # border radius [m]
    P_cap  <- sigma_film / r_pb                                      # capillary suction [Pa]
    Pi_bar <- Pi_charge * theta_cov                                  # disjoining barrier [Pa]
    h_eq   <- min(lambda_D * log(max(Pi_bar / P_cap, 1 + 1e-9)), h_eq_max)  # eq film thickness [m]

    # (2) BUBBLE POPULATION: coalescence (dry/unstable -> faster) - breakage
    coal_rate <- K_coal * (eps_g / eps_g_pack) / film_stability     # [1/s]

    # (3) BURST TRIGGER = FILM RUPTURE. The critical burst size is set by film
    # physics: weak films (low Gibbs elasticity) or high surface tension thin the
    # films and let smaller bubbles rupture, so d_b_crit shrinks -> burst sooner.
    # Rate is modulated by drainage speed (viscosity) and particle effects:
    #   - viscosity  : higher mu -> slower film drainage -> LESS burst
    #   - solids content : particles armor the film -> LESS burst
    #   - coarse-particle size : large particles bridge/rupture the film -> MORE burst
    d_b_crit <- d_b_burst * film_stability^a_fs * (sigma_ref_film / sigma_film)^a_sig
    M_visc   <- (mu_drain_ref / mu_cont)^n_visc                      # low mu -> faster drainage -> more burst
    M_armor  <- 1 / (1 + k_armor * phi_cond)                        # solids content armors film
    M_bridge <- 1 + k_bridge * (d_crs / max(d_b, d_b_in)) * (eps_s / eps_s_in)  # coarse bridging
    burst_rate <- K_burst * max(0, d_b - d_b_crit) / d_b_in * M_visc * M_armor * M_bridge

    # bursting removes the coarsest bubbles -> pulls the mean d_b back down, so
    # coarsening self-limits (bubbles burst before growing without bound).
    dd_b <- (coal_rate * d_b - K_break * (d_b - d_b_in) - K_bsink * burst_rate * d_b) / U

    burst_flux <- burst_rate * J_g_foam / U
    dJ_g_foam <- -burst_flux                                         # dispersed foam gas
    dJ_g_slug <-  burst_flux                                         # retained slug gas (out the top)

    # DRAINAGE: liquid drains through the Plateau-border network, whose
    # permeability scales with r_pb^2 -> wider borders (bigger bubbles / wetter
    # foam) drain FASTER. tau_local = tau_drain * (r_pb_ref/r_pb)^p_perm.
    tau_local <- tau_drain * (r_pb_ref / r_pb)^p_perm
    deps_l <- -(eps_l - eps_l_dry) / tau_local / U

    # (1) POOL decant with HINDERED settling (Richardson-Zaki) + KD viscosity
    phi <- 1 - eps_g                                                 # condensed holdup
    U_settle_sep <- stokes(rho_liquid - rho_foam, d_sep, mu_eff) * (1 - phi)^n_RZ
    dCimp_pool <- -(U_settle_sep / H_pool) * C_imp

    # POOL particle settling (aggregate; uses current d_b and crowded mu_eff)
    settle_loss <- function(dp) {
      d_agg <- (dp^3 + d_b^3)^(1/3)
      Uset  <- stokes(rho_agg(dp, d_b, rho_p, rho_gas) - rho_liquid, d_agg, mu_eff) * (1 - phi)^n_RZ
      max(0, Uset) / H_pool
    }
    sl_fine <- settle_loss(d_fine); sl_mid <- settle_loss(d_mid); sl_crs <- settle_loss(d_crs)

    # FOAM particle loss: detachment + BURST-driven collapse (dumps coarse)
    lf_fine <- k_det_fine + k_burstloss_fine * burst_rate
    lf_mid  <- k_det_mid  + k_burstloss_mid  * burst_rate
    lf_crs  <- k_det_crs  + k_burstloss_crs  * burst_rate

    dJs_fine <- -Js_fine * ((1-foam_frac)*sl_fine + foam_frac*lf_fine) / U
    dJs_mid  <- -Js_mid  * ((1-foam_frac)*sl_mid  + foam_frac*lf_mid ) / U
    dJs_crs  <- -Js_crs  * ((1-foam_frac)*sl_crs  + foam_frac*lf_crs ) / U

    # FOAM washing (channeling scales with current bubble size)
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    dCimp_foam <- -k_drainage * Wash_Efficiency * (wash_ratio/(wash_ratio+1)) * C_imp
    dCimp <- ((1-foam_frac)*dCimp_pool + foam_frac*dCimp_foam) / U

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l, dd_b, dJ_g_foam, dJ_g_slug, dt_res),
         h_eq = h_eq, r_pb = r_pb, tau_local = tau_local,
         d_b_crit = d_b_crit, burst_rate = burst_rate)
  })
}

# --- helpers ---------------------------------------------------------------
regime_of <- function(z, load_rel, p) {
  if (z < p[["H_pool"]]) return("pool")
  idx <- p[["film_stability"]] / (p[["load_sens"]] * load_rel)
  if (idx >= p[["plug_crit"]]) "plug" else "snowglobe"
}
solve_column <- function(p) {
  s0 <- c(Js_fine=p[["Js_fine_in"]], Js_mid=p[["Js_mid_in"]], Js_crs=p[["Js_crs_in"]],
          C_imp=p[["C_imp_in"]], eps_l=p[["eps_l_in"]], d_b=p[["d_b_in"]],
          J_g_foam=p[["eps_g_in"]]*p[["U_up"]], J_g_slug=0, t_res=0)
  z <- seq(0, p[["H_total"]], by = 0.02)
  out <- as.data.frame(ode(s0, z, column_model_psd, p)); names(out)[1] <- "z"
  lr <- (out$Js_fine+out$Js_mid+out$Js_crs)/(p[["Js_fine_in"]]+p[["Js_mid_in"]]+p[["Js_crs_in"]])
  out$J_g_total <- out$J_g_foam + out$J_g_slug     # conserved: carried out the top
  out$eps_s <- p[["eps_s_in"]]*lr
  out$eps_g <- 1 - out$eps_l - out$eps_s
  out$solids_pct <- 100*out$eps_s/(out$eps_s+out$eps_l)
  out$regime <- vapply(seq_len(nrow(out)), function(i) regime_of(out$z[i], lr[i], p), character(1))
  out
}

# --- solve + report --------------------------------------------------------
dpar    <- derive_state_props(params)     # T/P + surfactant -> effective closures
out_psd <- solve_column(dpar)
t_pool <- approx(out_psd$z, out_psd$t_res, dpar[["H_pool"]])$y
t_tot  <- out_psd$t_res[nrow(out_psd)]
Jg_in  <- dpar[["eps_g_in"]]*dpar[["U_up"]]; Jg_out <- out_psd$J_g_total[nrow(out_psd)]
top <- out_psd[nrow(out_psd),]; ini <- out_psd[1,]

cat(sprintf("Thermo: T %.1f K, P %.2f bar | mu(T) %.2f mPa.s, rho_gas %.2f kg/m3, sigma %.1f mN/m\n",
            dpar[["T_col"]], dpar[["P_col"]]/1e5, dpar[["mu_cont"]]*1e3, dpar[["rho_gas"]], dpar[["sigma_film"]]*1e3))
cat(sprintf("Surfactant: c %.1f mol/m3 (CMC %.1f) -> coverage %.2f, Gibbs E %.0f mN/m -> film_stability %.2f\n",
            dpar[["c_surf"]], dpar[["cmc"]], dpar[["theta_cov"]], dpar[["E_gibbs"]]*1e3, dpar[["film_stability"]]))

# Film thickness / Plateau-border geometry (wired into drainage + holdup)
cat(sprintf("Film: eq. thickness h %.0f nm (sets holdup, wet-factor %.2f) | border r_pb %.0f->%.0f um up z\n",
            dpar[["h_eq_dry"]]*1e9, dpar[["film_wet"]], ini$r_pb*1e6, top$r_pb*1e6))
cat(sprintf("      drainage tau %.0f->%.0f s up z (wider borders drain faster) | eq. holdup %.3f\n",
            ini$tau_local, top$tau_local, dpar[["eps_l_dry"]]))

cat(sprintf("Bubble size d_b: inlet %.2f mm -> top %.2f mm (coalescence, burst-limited)\n",
            dpar[["d_b_in"]]*1e3, top$d_b*1e3))
d_b_crit0 <- dpar[["d_b_burst"]] * dpar[["film_stability"]]^dpar[["a_fs"]] *
             (dpar[["sigma_ref_film"]]/dpar[["sigma_film"]])^dpar[["a_sig"]]
cat(sprintf("Burst trigger: film-rupture critical size d_b_crit %.2f mm (film_stab %.2f, sigma %.1f mN/m)\n",
            d_b_crit0*1e3, dpar[["film_stability"]], dpar[["sigma_film"]]*1e3))
cat(sprintf("GAS balance (conserved, no vent): in %.2e -> out %.2e m/s => carried out top %.0f%%\n",
            Jg_in, Jg_out, 100*Jg_out/Jg_in))
cat(sprintf("  gas state at top: foam bubbles %.0f%% + retained slug %.0f%% (slug released next unit up)\n",
            100*top$J_g_foam/Jg_in, 100*top$J_g_slug/Jg_in))
cat(sprintf("Residence: pool %.2f h + foam %.2f h = %.2f h\n",
            t_pool/3600, (t_tot-t_pool)/3600, t_tot/3600))
cat(sprintf("Solids: inlet %.1f%% -> top %.1f%% | plug %.0f%%\n",
            ini$solids_pct, top$solids_pct, 100*mean(out_psd$regime=="plug")))
cat(sprintf("Retention: fine=%.0f%% mid=%.0f%% coarse=%.0f%% | impurity %.0f->%.1f%%\n",
            100*top$Js_fine/ini$Js_fine, 100*top$Js_mid/ini$Js_mid,
            100*top$Js_crs/ini$Js_crs, ini$C_imp, top$C_imp))

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot(out_psd$z, out_psd$Js_mid, type="l", col="green", lwd=3, ylim=c(0,0.05),
     ylab="Solid flux (m/s)", xlab="Height (m)", main="PSD retention")
lines(out_psd$z, out_psd$Js_fine, col="blue", lwd=2, lty=2)
lines(out_psd$z, out_psd$Js_crs,  col="red",  lwd=2, lty=3)
abline(v=dpar[["H_pool"]], col="gray", lty=2)
legend("right", legend=c("Mid","Fine","Coarse"), col=c("green","blue","red"), lty=c(1,2,3), lwd=2)

plot(out_psd$z, out_psd$d_b*1e3, type="l", col="darkorange", lwd=3,
     ylim=c(0, max(out_psd$d_b*1e3, out_psd$d_b_crit*1e3)),
     ylab="Bubble d_b (mm)", xlab="Height (m)", main="Bubble growth & burst trigger")
lines(out_psd$z, out_psd$d_b_crit*1e3, col="red", lwd=2, lty=3)
legend("topleft", legend=c("d_b","d_b_crit (film rupture)"),
       col=c("darkorange","red"), lty=c(1,3), lwd=2)

plot(out_psd$z, out_psd$J_g_total/Jg_in*100, type="l", col="black", lwd=3, ylim=c(0,100),
     ylab="Gas (% of inlet)", xlab="Height (m)", main="Gas state (conserved)")
lines(out_psd$z, out_psd$J_g_foam/Jg_in*100, col="purple",    lwd=2, lty=2)
lines(out_psd$z, out_psd$J_g_slug/Jg_in*100, col="firebrick", lwd=2, lty=3)
abline(v=dpar[["H_pool"]], col="gray", lty=2)
legend("right", legend=c("Total (out top)","Foam bubbles","Slug (retained)"),
       col=c("black","purple","firebrick"), lty=c(1,2,3), lwd=2)

# Film thickness (wired into drainage/holdup) + liquid holdup
plot(out_psd$z, out_psd$h_eq*1e9, type="l", col="steelblue", lwd=3,
     ylab="Film h_eq (nm)", xlab="Height (m)", main="Film thickness & holdup")
par(new=TRUE)
plot(out_psd$z, out_psd$eps_l, type="l", col="darkgreen", lwd=2, lty=2,
     axes=FALSE, xlab="", ylab="", ylim=c(0, max(out_psd$eps_l)))
axis(4); mtext("Liquid holdup eps_l (-)", side=4, line=-1.3, cex=0.7)
abline(v=dpar[["H_pool"]], col="gray", lty=2)
legend("right", legend=c("Film h_eq (L)","Holdup eps_l (R)"),
       col=c("steelblue","darkgreen"), lty=c(1,2), lwd=2)
