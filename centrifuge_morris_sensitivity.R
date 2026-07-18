# =============================================================================
# Morris (Elementary Effects) sensitivity analysis of a unified decanter-
# centrifuge model for a three-phase polymer separation.
#
# FACTORS ARE ORGANISED INTO SWEEP GROUPS (see the `factors` table below):
#   * equipment  - hard-to-change machine geometry (bowl diameter, cylinder
#                  length, beach angle, feed-zone length). MUTED by default
#                  (held at nominal) but available to sweep.
#   * process    - operator setpoints: rpm / g-force, scroll differential,
#                  weir (pool radius), flow, feed composition, degassing,
#                  temperature, feed pressure, dissolved-gas loading.
#   * surface    - surface chemistry: surfactant HLB, surfactant dose, foam
#                  thixotropic recovery.
#   * additive   - formulation additives: plasticizer, residual monomer,
#                  binder, and the Flory-Huggins chi compatibility.
# Set `active_groups` (Section 2) to choose which groups are swept; factors in
# inactive groups are pinned at their nominal value. This lets you screen one
# group in isolation or all process-side factors together.
#
# THREE-PHASE separation: light (recovered) + heavy (recovered) + centrate
# (lost) => Yield_Loss_Centrate = 1 - (Heavy_Solid_Yield + Light_Solid_Yield).
#
# ZONE-SEPARATED RESIDENCE: feed/acceleration (length_feed), cylindrical
# clarification (L_cyl -> t_clarify, where settling is captured), beach
# dewatering (t_gap). Feed-zone acceleration intensifies aggregate rupture.
# Outlet classification: per-size capture -> grade-efficiency cut size d50c
# and cake solids median d50.
#
# Carrier is water (temperature-dependent, Vogel); the discharged wet cake /
# mousse is a Herschel-Bulkley yield-stress material. Thermodynamic layer:
# Henry's-law dissolved-gas flash template + Clausius-Clapeyron monomer
# volatility (both need the explicit T / P state).
#
# SURFACE CHEMISTRY IS COUPLED: the surfactant monolayer coverage (theta_surf,
# from dose vs interfacial demand) now drives foam stability, which (a) retards
# the mechanical pop degassing, (b) slows centrifugal FOAM DRAINAGE, (c) boosts
# the foam-cushion shear protection, and (d) keeps the light (foam-borne) solid
# phase from shedding to the centrate. Centrate_Foam_Ratio remains the residual
# downstream foam-potential output.
#
# FOAM DRAINAGE: liquid drains from the foam films under the centrifugal body
# force (Plateau-border drainage), retarded by the surfactant-rigidified
# interfaces AND the bulk film viscosity (polymer / binder / plasticizer). A
# high-g bowl breaks the foam (extra degassing) unless a stable foam resists it.
#
# FOAM STABILITY is a composite of four mechanisms, all driven by feed /
# upstream colloidal properties: surfactant coverage, (1) Pickering armoring by
# fine solids (contact_angle_deg), (2) DLVO disjoining pressure (Delta_pH vs the
# isoelectric point, screened by I_strength), (3) Gibbs-Marangoni film
# elasticity, and (4) the bubble size D_b (which coarsens by coalescence /
# Ostwald ripening over the pressurized hold). (5) Electrolyte (I_strength) also
# salts out the dissolved gas (Sechenov) and raises the clean surface tension.
#
# UPSTREAM/DOWNSTREAM: D_b, I_strength, Delta_pH and contact_angle_deg are feed
# properties from the prior step; the model in turn hands downstream the coarsened
# Bubble_Size_um, the Serum_Surface_Tension_N_m and the Foam_Stability state.
#
# Method: Morris OAT screening on a [0,1] hypercube over the ACTIVE factors
# (comparable elementary effects), mapped to physical ranges, all outputs from
# one design. CSV + mu* vs sigma panel grid written to output/.
#
# Uses the `sensitivity` package when installed; otherwise falls back to a
# built-in base-R Morris OAT design (no hard dependency). A smoke test evaluates
# the model at all-nominal first.  Run: Rscript centrifuge_morris_sensitivity.R
# =============================================================================

# install.packages("sensitivity")   # optional; a base-R Morris fallback is built in

# ---- Temperature-dependent physical property helpers ------------------------
mu_water_T <- function(T) 2.414e-5 * 10^(247.8 / (T - 140))          # water viscosity [Pa s]
k_H_T      <- function(T) 6.5e-6 * exp(1700 * (1/T - 1/298))          # Henry const [mol/(m3 Pa)]

# =========================================================================
# 1. THE UNIFIED COMPREHENSIVE CENTRIFUGE ENGINE
# =========================================================================
unified_centrifuge_model <- function(run) {

  # ---------------------------------------------------------
  # A. EXTRACT PARAMETERS (all groups; muted ones arrive at nominal)
  # ---------------------------------------------------------
  # process
  rpm              <- run[["rpm"]]
  delta_rpm        <- run[["delta_rpm"]]
  r_pool           <- run[["r_pool"]]            # weir / pond radius       [m]
  flow_rate        <- run[["flow_rate_lpm"]]
  feed_solid_frac  <- run[["feed_solid_frac"]]
  feed_gas_frac    <- run[["feed_gas_frac"]]
  pop_frac         <- run[["pop_frac"]]
  T_process        <- run[["T_process"]]         # process temperature      [K]
  P_feed           <- run[["P_feed"]]            # feed-line pressure       [Pa]
  gas_sat_frac     <- run[["gas_sat_frac"]]      # dissolved-gas saturation [-]
  D_b              <- run[["D_b"]]               # feed bubble diameter     [m]
  # surface chemistry / colloidal
  HLB_value        <- run[["HLB_value"]]
  surf_dose_kg_m3  <- run[["surf_dose_kg_m3"]]   # surfactant dose          [kg/m3]
  t_recover_sec    <- run[["t_recover_sec"]]     # thixotropic foam recovery[s]
  I_strength       <- run[["I_strength"]]        # ionic strength (Debye)   [mol/L]
  Delta_pH         <- run[["Delta_pH"]]          # pH offset from IEP       [-]
  contact_angle_deg<- run[["contact_angle_deg"]] # particle wettability (Pickering)
  # additives / material
  plasticizer_frac <- run[["plasticizer_frac"]]
  C_monomer        <- run[["C_monomer"]]
  C_binder         <- run[["C_binder"]]          # binder concentration     [wt/wt]
  chi_parameter    <- run[["chi_parameter"]]
  S_base           <- run[["S_base"]]            # compacted-sludge solids floor (vol frac)
  S_ceiling        <- run[["S_ceiling"]]         # max dewatered solids (vol frac)
  # equipment (muted by default)
  L_cyl            <- run[["L_cyl"]]             # cylinder length          [m]
  r_bowl           <- run[["r_bowl"]]            # bowl radius (diameter)   [m]
  beach_angle_deg  <- run[["beach_angle_deg"]]
  length_feed      <- run[["length_feed"]]       # feed-zone length         [m]

  # ---------------------------------------------------------
  # B. FIXED CONSTANTS & DERIVED GEOMETRY
  # ---------------------------------------------------------
  r_discharge <- 0.110        # solids discharge radius (fixed) [m]; deeper-pond bowl
  scroll_pitch <- 0.050       # scroll pitch (fixed)            [m]
  beach_angle <- beach_angle_deg * pi / 180
  # keep the pond physical even if equipment geometry is swept
  r_pool <- min(max(r_pool, r_discharge + 0.005), r_bowl - 0.005)

  rho_poly <- 1800; rho_liq <- 1000     # dry solids density (no air) [kg/m3]; water
  R_gas <- 8.314; R_air <- 287
  P_atm <- 101325; omega <- rpm * (2 * pi / 60)

  visc_water <- mu_water_T(T_process)

  # Core/template density coupled to plasticizer (theta-solvent core)
  rho_plast_theta <- 820
  plast_lo <- 0.05; plast_hi <- 0.25
  plast_norm <- max(0, min(1, (plasticizer_frac - plast_lo) / (plast_hi - plast_lo)))
  template_density <- rho_liq - (rho_liq - rho_plast_theta) * plast_norm

  # ---------------------------------------------------------
  # C. MASS BALANCE, SURFACE CHEMISTRY, FOAM DRAINAGE & DEGASSING
  # ---------------------------------------------------------
  Q_in_m3s <- flow_rate / 60000
  feed_liq_frac <- max(0, 1.0 - (feed_solid_frac + feed_gas_frac))

  Q_solid <- Q_in_m3s * feed_solid_frac
  Q_liq   <- Q_in_m3s * feed_liq_frac
  Q_gas_atm <- Q_in_m3s * feed_gas_frac
  rho_gas <- P_atm / (R_air * T_process)

  # Zone-separated residence. Axial transport is liquid-dominated (entrained
  # gas rises to the air core), so the pond velocity uses the liquid+solid flow.
  A_pond <- pi * (r_bowl^2 - r_pool^2)
  Q_ls   <- Q_solid + Q_liq
  U_ax   <- Q_ls / A_pond
  t_feed    <- length_feed / U_ax
  t_clarify <- L_cyl / U_ax
  t_pond    <- t_feed + t_clarify
  dist <- r_bowl - r_pool

  # Hydrostatic (centrifugal) pressure field
  P_hydro <- 0.5 * rho_liq * omega^2 * (r_bowl^2 - r_pool^2)
  P_local <- P_atm + P_hydro

  # --- Surfactant monolayer coverage (feed interfacial demand) ----------
  # Low-HLB surfactant prefers the (hydrophobic) polymer surface; bubble area
  # uses the feed bubble size D_b (carried from the upstream step).
  MW_surf <- 500; a0_angstrom2 <- 50; CMC_kg_m3 <- 1.0
  kg_per_m2 <- (MW_surf / 1000) / ((a0_angstrom2 * 1e-20) * 6.022e23)
  affinity_poly <- 1 / (1 + exp((HLB_value - 10) / 2))       # adsorption onto polymer
  SA_poly_feed  <- (6 / 50e-6) * Q_solid                     # solids interfacial area
  SA_gas_feed   <- (6 / D_b) * Q_gas_atm                     # bubble interfacial area
  surf_demand   <- (SA_gas_feed + SA_poly_feed * affinity_poly) * kg_per_m2
  surf_supply   <- Q_liq * surf_dose_kg_m3
  theta_surf    <- min(1, surf_supply / max(surf_demand, 1e-12))   # fractional coverage
  foam_HLB      <- 1 / (1 + exp(-(HLB_value - 10) / 2))            # hydrophilic -> foams

  # --- Composite foam stability: four mechanisms ------------------------
  # (3) Gibbs-Marangoni film elasticity: peaks at intermediate coverage.
  elasticity <- 4 * theta_surf * (1 - theta_surf)
  # (2) DLVO disjoining pressure: surface charge (pH offset from the isoelectric
  #     point) stabilizes the thin films; ionic strength screens it (Debye).
  #     Both Delta_pH and I_strength are feed / upstream colloidal properties.
  I_ref <- 0.01; dpH_ref <- 2.0
  debye_screen <- 1 / (1 + sqrt(I_strength / I_ref))
  dlvo_raw <- (Delta_pH / dpH_ref) * debye_screen
  dlvo <- dlvo_raw / (1 + dlvo_raw)                           # bounded 0..1
  # (1) Pickering: fine solids adsorb at the gas-liquid interface and armor the
  #     films, most effective near a 90 deg contact angle.
  pickering_wet <- exp(-((contact_angle_deg - 90) / 40)^2)
  solid_avail   <- SA_poly_feed / max(SA_gas_feed, 1e-12)
  pickering <- min(1, 0.5 * solid_avail * pickering_wet)
  # composite, gated by a foam-favorable (hydrophilic) HLB
  stab_raw <- 0.35 * theta_surf + 0.20 * elasticity + 0.20 * dlvo + 0.25 * pickering
  foam_stability <- foam_HLB * stab_raw                        # 0..1

  # --- (5) Electrolyte + surfactant surface tension ---------------------
  # Salt raises the clean surface tension slightly; surfactant lowers it toward
  # the CMC plateau. sigma_eff is the surface-tension state handed downstream.
  sigma_cmc   <- 0.030
  sigma_clean <- 0.072 + 0.02 * I_strength / (I_strength + I_ref)
  sigma_eff   <- sigma_clean - (sigma_clean - sigma_cmc) * theta_surf

  # Film viscosity resisting drainage (polymer makeup + stabilizers + T)
  mu_film <- visc_water * (1 + 12 * C_binder + 3 * plasticizer_frac)

  # --- (4) Bubble size: coalescence / Ostwald ripening over the hold -----
  # Feed bubbles coarsen during the pressurized residence, faster under pressure
  # (Henry) and retarded by a well-covered film. Feeds drainage and foam yield.
  k_rip0 <- 4e-15
  k_rip  <- k_rip0 * (P_local / P_atm) * (1 - 0.8 * theta_surf)
  D_b_eff <- (D_b^3 + k_rip * t_pond)^(1/3)
  R_bub   <- D_b_eff / 2

  # --- Foam drainage (centrifugal Plateau-border drainage) --------------
  # Body-force drainage retarded by interfacial rigidity (surfactant) and bulk
  # film viscosity; drained foam collapses and vents its gas to the air core.
  g_eff    <- r_pool * omega^2
  eps_foam <- 0.20
  r_pb     <- R_bub * sqrt(eps_foam)
  C_drain  <- 50 * (1 + 20 * theta_surf)
  v_drain  <- rho_liq * g_eff * r_pb^2 / (C_drain * mu_film)
  t_drain  <- dist / max(v_drain, 1e-9)
  foam_drained_frac <- 1 - exp(-t_pond / t_drain)

  # --- Drift-flux slip: buoyant rise of FREE bubbles to the air core -----
  # Bubbles are far lighter than the serum, so in the centrifugal field they
  # slip radially inward at a Stokes rise velocity relative to the liquid and
  # escape to the air core (the light-phase analogue of solids settling). This
  # is a genuine gas-liquid slip term, separate from foam drainage.
  u_slip     <- D_b_eff^2 * (rho_liq - rho_gas) * g_eff / (18 * mu_film)
  f_gas_rise <- 1 - exp(-u_slip * t_pond / dist)          # free-bubble escape

  # --- Degassing: pop, then THREE parallel fates for the entrained gas ---
  # After mechanical popping the gas partitions three ways:
  #   * particle-attached (flotation): bubbles attach to the hydrophobic polymer
  #     particles and RIDE WITH THE SOLIDS - they do not rise freely, so this
  #     gas is RETAINED in the product. This is the residual air that keeps the
  #     exit slurry buoyant (still floats in water, denser than a gravity foam).
  #   * foam-bound: leaves by drainage collapse.
  #   * free bubbles: leave by the drift-flux buoyant rise above.
  pop_frac_eff <- pop_frac * (1 - 0.6 * foam_stability)
  Q_gas_after_pop <- Q_gas_atm * (1 - pop_frac_eff)
  # flotation attachment: hydrophobic (low-HLB) particles at a favourable
  # contact angle capture the most gas
  p_attach   <- min(1, pickering_wet * affinity_poly)
  Q_gas_attach <- Q_gas_after_pop * p_attach              # retained with solids
  Q_gas_disp   <- Q_gas_after_pop * (1 - p_attach)
  Q_gas_foam   <- Q_gas_disp * foam_stability
  Q_gas_free   <- Q_gas_disp * (1 - foam_stability)
  Q_gas_surviving_atm <- Q_gas_attach +
                         Q_gas_foam * (1 - foam_drained_frac) +
                         Q_gas_free * (1 - f_gas_rise)

  # Compression (Boyle)
  Q_gas_local <- Q_gas_surviving_atm * (P_atm / P_local)

  # Local suspension properties
  Q_local_m3s <- Q_solid + Q_liq + Q_gas_local
  rho_feed_local <- ((Q_solid * rho_poly) + (Q_liq * rho_liq) + (Q_gas_surviving_atm * rho_gas)) / Q_local_m3s

  phi_dispersed <- (Q_solid + Q_gas_local) / Q_local_m3s
  safe_phi <- min(phi_dispersed, 0.64)
  mu_apparent <- visc_water * (1 - (safe_phi / 0.65))^(-2.5 * 0.65)

  # Foam cushion: thixotropic recovery over the residence, boosted by the
  # surfactant-stabilized film network.
  foam_dampening_full <- max(1.0, mu_apparent / visc_water)
  recovery_frac <- 1 - exp(-t_pond / t_recover_sec)
  foam_dampening_factor <- 1.0 + (foam_dampening_full - 1.0) * recovery_frac * (1 + 1.5 * foam_stability)

  k_feed_accel <- 80
  n_rev_feed   <- omega * t_feed / (2 * pi)
  feed_shear_factor <- 1 + k_feed_accel / max(n_rev_feed, 1e-6)

  # ---------------------------------------------------------
  # C.5 THERMODYNAMIC FLASHING (T, P dependent)
  # ---------------------------------------------------------
  # (5) electrolyte salting-out lowers gas solubility (Sechenov)
  kH <- k_H_T(T_process) * 10^(-0.2 * I_strength)
  C_gas_loaded <- gas_sat_frac * kH * P_feed
  C_gas_final  <- kH * P_atm
  flash_mol_per_m3 <- max(0, C_gas_loaded - C_gas_final)
  V_flash_per_liq  <- flash_mol_per_m3 * R_gas * T_process / P_atm

  Hvap_mono <- 32000; T_bp_mono <- 325
  p_mono <- P_atm * exp((Hvap_mono / R_gas) * (1 / T_bp_mono - 1 / T_process))
  mono_flash_frac <- 1 / (1 + exp(-(p_mono - P_atm) / (0.2 * P_atm)))
  C_mono_retained <- C_monomer * (1 - mono_flash_frac)

  # ---------------------------------------------------------
  # D. PARTICLE TRACKING: SHEAR RUPTURE, MIGRATION & CLASSIFICATION
  # ---------------------------------------------------------
  # Log-spaced primary size grid down to sub-micron so the fine tail that
  # escapes to the centrate (and sets the classification cut size) is resolved.
  size_bins <- exp(seq(log(0.5), log(200), length.out = 24))
  core_bins <- seq(0.1, 0.6, length.out = 15)
  n_size <- length(size_bins)

  # log-normal feed PSD (median ~25 um) with a real fine tail
  size_pdf   <- dnorm(log(size_bins), log(25), 1.0)
  weight_sum <- sum(size_pdf %o% dnorm(core_bins, 0.3, 0.1))

  heavy_solid_yield <- 0
  light_solid_yield <- 0
  intact_template_yield <- 0
  agg_survival <- 0
  size_capture   <- numeric(n_size)
  size_feed_mass <- numeric(n_size)

  phi_s_local <- Q_solid / Q_local_m3s
  hinder <- (max(0, 1 - phi_s_local))^4.65
  # dist (settling distance) already defined in Section C

  for (i in seq_len(n_size)) {
    d_um <- size_bins[i]
    for (core_frac in core_bins) {

      weight <- size_pdf[i] * dnorm(core_frac, 0.3, 0.1)
      if(is.na(weight) || weight <= 0) next

      d_primary_m <- d_um * 1e-6
      agg_multiplier <- 1 + (10 * plasticizer_frac)
      d_intact_m <- d_primary_m * agg_multiplier

      # Shear rupture: plasticizer AND binder bind the aggregate; feed-zone
      # acceleration intensifies the disrupting shear.
      cohesive_strength <- 1000 + (50000 * plasticizer_frac) + (30000 * C_binder)
      shear_stress <- (abs(rho_poly - template_density) * r_pool * omega^2 * d_intact_m) *
                      feed_shear_factor / foam_dampening_factor
      survival_prob <- 1 / (1 + exp(-(cohesive_strength - shear_stress) / 500))
      agg_survival <- agg_survival + (survival_prob * weight)

      rho_intact <- (core_frac * template_density) + ((1 - core_frac) * rho_poly)
      v_c_intact <- (d_intact_m^2 * (rho_intact - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent) * hinder
      d_ruptured <- d_primary_m * (agg_multiplier * 0.5)
      v_c_rupt   <- (d_ruptured^2 * (rho_poly - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent) * hinder

      rec_intact_heavy <- ifelse(v_c_intact > 0, 1 - exp(-abs(v_c_intact) * t_clarify / dist), 0)
      rec_intact_light <- ifelse(v_c_intact < 0, 1 - exp(-abs(v_c_intact) * t_clarify / dist), 0)
      rec_rupt_heavy   <- ifelse(v_c_rupt > 0,   1 - exp(-abs(v_c_rupt) * t_clarify / dist), 0)

      heavy_rec <- (rec_intact_heavy * survival_prob) + (rec_rupt_heavy * (1 - survival_prob))

      intact_template_yield <- intact_template_yield + ((rec_intact_heavy + rec_intact_light) * survival_prob * weight)
      heavy_solid_yield <- heavy_solid_yield + (heavy_rec * weight)
      light_solid_yield <- light_solid_yield + (rec_intact_light * survival_prob * weight)

      size_capture[i]   <- size_capture[i] + (heavy_rec * weight)
      size_feed_mass[i] <- size_feed_mass[i] + weight
    }
  }

  intact_template_yield <- intact_template_yield / weight_sum
  heavy_solid_yield     <- heavy_solid_yield / weight_sum
  light_solid_yield     <- light_solid_yield / weight_sum
  agg_survival          <- agg_survival / weight_sum

  # Light (foam-borne) solids only report to the light phase while the foam
  # survives; drained / collapsed foam sheds them to the centrate (a loss).
  foam_vehicle <- 0.3 + 0.7 * (1 - foam_drained_frac)
  light_solid_yield <- light_solid_yield * foam_vehicle

  # Outlet classification: grade efficiency -> cut size + cake median
  grade_eff <- size_capture / pmax(size_feed_mass, 1e-12)
  feed_mean_size <- sum(size_bins * size_feed_mass) / sum(size_feed_mass)
  d50c <- tryCatch(approx(grade_eff, size_bins, xout = 0.5, ties = "ordered", rule = 2)$y,
                   error = function(e) NA_real_)
  if (is.na(d50c)) d50c <- feed_mean_size
  if (sum(size_capture) > 1e-12) {
    cake_cdf <- cumsum(size_capture) / sum(size_capture)
    cake_d50 <- tryCatch(approx(cake_cdf, size_bins, xout = 0.5, ties = "ordered", rule = 2)$y,
                         error = function(e) feed_mean_size)
  } else {
    cake_d50 <- feed_mean_size
  }

  # ---------------------------------------------------------
  # E. WET-CAKE / MOUSSE RHEOLOGY & RETENTION
  # ---------------------------------------------------------
  v_convey <- (delta_rpm / 60) * scroll_pitch
  # Beach dewatering residence. Per plant observation, a SMALLER beach angle
  # gives a SHORTER effective dry-drainage zone -> WETTER cake (and a steeper
  # beach drains more -> drier). This is opposite to the naive "shallow cone =
  # long beach" view; decanter beach-angle conventions/designs differ, so the
  # dry-beach length is anchored to a 10 deg reference and scales with the angle.
  beach_ref <- 10 * pi / 180
  beach_axial <- (r_pool - r_discharge) / tan(beach_ref) * (tan(beach_angle) / tan(beach_ref))
  t_gap <- beach_axial / v_convey

  # --- Dewatering: cake solids content (gas-free basis) -----------------
  # Solids rise toward the achievable ceiling (~50 %, a decanter dewatering
  # limit) with beach residence t_gap; tacky/plasticized/binder cakes and a
  # high serum surface tension hold moisture, while surfactant (low sigma) aids
  # dewatering. Cake solids cannot exceed the ceiling.
  # S_base (compaction floor) and S_ceiling (max dewatered) are material factors
  S_ceiling <- max(S_ceiling, S_base + 0.02)          # keep the ceiling above the floor
  dewater <- 1 - exp(-t_gap / 8.0)                     # residence approach to ceiling
  aid <- 1 - 0.40 * (plasticizer_frac / 0.25) - 0.40 * (C_binder / 0.05) -
             0.20 * ((sigma_eff / sigma_clean) - 0.33) / 0.67
  aid <- max(0.2, min(1, aid))
  solids_gasfree <- S_base + (S_ceiling - S_base) * dewater * aid   # <= 0.50
  liq_gasfree    <- 1 - solids_gasfree

  # --- Cake gas: sparged holdup + dissolved-gas flash (carried in liquid) --
  gas_frac_local   <- Q_gas_local / Q_local_m3s
  cake_gas_sparged <- (0.02 + 0.20 * exp(-t_gap / 3.0)) * (0.5 + 2.5 * gas_frac_local)
  cake_gas_flash   <- V_flash_per_liq * liq_gasfree + 0.5 * C_monomer * mono_flash_frac
  cake_gas_frac    <- min(0.5, cake_gas_sparged + cake_gas_flash)

  # Compose the three phases (solids stay at / below the ~50 % ceiling)
  cake_solid_frac    <- solids_gasfree * (1 - cake_gas_frac)
  cake_moisture_frac <- liq_gasfree * (1 - cake_gas_frac)

  # Herschel-Bulkley wet-cake rheology (yield-stress paste / mousse); binder
  # bridging stiffens the frictional network.
  phi_pack <- 0.64; p_exp <- 3.0
  tau_y0   <- 5.0e4          # frictional yield scale [Pa]  (R_bub, sigma_eff from Section C)
  n_hb <- 0.5; K_hb <- 200; gamma_disc <- 10

  phi_ratio  <- min(cake_solid_frac, 0.98 * phi_pack) / phi_pack
  tau_y_fric <- tau_y0 * phi_ratio^p_exp * (1 + 5.0 * C_binder)
  tau_y_foam <- (sigma_eff / R_bub) * max(cake_gas_frac - 0.10, 0)   # Princen foam yield
  soften     <- exp(-2.0 * (plasticizer_frac + C_mono_retained)) * exp(-0.02 * (T_process - 293))
  cake_yield_stress <- (tau_y_fric + tau_y_foam) * soften

  paste_viscosity_pa_s <- cake_yield_stress / gamma_disc + K_hb * gamma_disc^(n_hb - 1)

  # Sprayability: apparent (Herschel-Bulkley) viscosity at an atomization shear
  # rate vs a practical sprayable ceiling. As dewatering pushes solids toward
  # 50 %, the yield stress climbs and the cake becomes too viscous to atomize
  # without dilution / heating downstream (Spray_Visc_Ratio > 1).
  gamma_spray    <- 1e4              # atomization shear rate      [1/s]
  mu_spray_limit <- 0.5              # practical atomization ceiling [Pa s]
  mu_spray <- cake_yield_stress / gamma_spray + K_hb * gamma_spray^(n_hb - 1)
  Spray_Visc_Ratio <- mu_spray / mu_spray_limit

  # Plasticizer Retention (Flory-Huggins)
  Q_liq_safe  <- max(Q_liq, 1e-12)
  Q_liq_heavy <- (Q_solid * heavy_solid_yield / cake_solid_frac) * cake_moisture_frac
  K_vol_plast <- exp(1 - 2 * chi_parameter)

  Cap_Poly <- K_vol_plast * Q_solid
  Cap_Liq  <- 1.0 * Q_liq
  Cap_Tot  <- max(Cap_Poly + Cap_Liq, 1e-12)
  Plast_frac_solid  <- Cap_Poly / Cap_Tot
  Plast_frac_liquid <- Cap_Liq / Cap_Tot
  Total_Plast_Retained <- (Plast_frac_solid * heavy_solid_yield) + (Plast_frac_liquid * (Q_liq_heavy / Q_liq_safe))

  # ---------------------------------------------------------
  # F. SURFACTANT PARTITIONING -> CENTRATE FOAM POTENTIAL (surface group)
  # ---------------------------------------------------------
  # Reuse the Section-C monolayer coverage: surfactant not adsorbed on the
  # polymer or bubble interfaces stays dissolved in the centrate, whose
  # concentration vs CMC gives the residual (downstream) foam potential.
  # The coverage also drives foam stability / drainage / degassing in Section C.
  Surf_on_Gas  <- SA_gas_feed  * kg_per_m2
  Surf_on_Poly <- SA_poly_feed * kg_per_m2 * affinity_poly
  Surf_Dissolved <- max(0, surf_supply - Surf_on_Gas - Surf_on_Poly)

  Q_liq_centrate <- Q_liq - Q_liq_heavy - (Q_gas_surviving_atm * 0.05)
  Centrate_Conc_kg_m3 <- Surf_Dissolved / max(Q_liq_centrate, 1e-6)
  Foam_Potential_Ratio <- Centrate_Conc_kg_m3 / CMC_kg_m3

  # ---------------------------------------------------------
  # F.5 GAS CARRIED FORWARD IN THE SLURRY (handoff to next unit)
  # ---------------------------------------------------------
  # Entrained gas retained in the RECOVERED product (heavy + light), mostly the
  # particle-attached share, as a volume fraction - the feed foam quality / gas
  # holdup the downstream unit inherits. Plus the dissolved gas (Henry) still in
  # solution at discharge, a latent flash source for the next unit.
  Q_solid_rec  <- Q_solid * (heavy_solid_yield + light_solid_yield)
  liq_to_solid <- cake_moisture_frac / cake_solid_frac         # cake L/S ratio
  Q_liq_rec    <- Q_solid_rec * liq_to_solid
  alpha_g_out  <- Q_gas_surviving_atm /
                  max(Q_gas_surviving_atm + Q_solid_rec + Q_liq_rec, 1e-12)
  dissolved_gas_mol_m3 <- min(C_gas_loaded, C_gas_final)

  # Exit density of the recovered product (gas-free matrix aerated by the
  # retained gas). Below 1000 kg/m3 -> still floats in water; a gravity-decanter
  # foam is ~600 kg/m3, so a denser-but-still-floating exit is expected.
  rho_pf <- (Q_solid_rec * rho_poly + Q_liq_rec * rho_liq) /
            max(Q_solid_rec + Q_liq_rec, 1e-12)
  exit_density_kg_m3 <- (1 - alpha_g_out) * rho_pf + alpha_g_out * rho_gas
  # product solids on a MASS basis (spray model's C_solid_mass), and the
  # continuous-phase serum viscosity handed downstream (spray model's mu_L)
  solids_mass_frac <- (Q_solid_rec * rho_poly) /
                      max(Q_solid_rec * rho_poly + Q_liq_rec * rho_liq, 1e-12)

  # ---------------------------------------------------------
  # G. RETURN ALL VARIABLES AS A LIST
  # ---------------------------------------------------------
  return(list(
    Intact_Template_Yield = intact_template_yield,
    Aggregate_Survival    = agg_survival,
    Aggregate_Destruction = 1 - agg_survival,
    Heavy_Solid_Yield     = heavy_solid_yield,
    Yield_Loss_Centrate   = 1 - (heavy_solid_yield + light_solid_yield),
    Cut_Size_d50c_um      = d50c,
    Cake_d50_um           = cake_d50,
    Wet_Cake_Moisture     = cake_moisture_frac,
    Beach_Residence_s     = t_gap,
    Gas_Template_Voidage  = cake_gas_frac,
    Entrained_Gas_Holdup  = alpha_g_out,
    Exit_Density_kg_m3    = exit_density_kg_m3,
    Gasfree_Density_kg_m3 = rho_pf,
    Product_Solids_MassFrac = solids_mass_frac,
    Serum_Viscosity_Pa_s  = mu_film,
    Bubble_Rise_Escape    = f_gas_rise,
    Dissolved_Gas_mol_m3  = dissolved_gas_mol_m3,
    Cake_Yield_Stress_Pa  = cake_yield_stress,
    Paste_Viscosity_Pa_s  = paste_viscosity_pa_s,
    Spray_Visc_Ratio      = Spray_Visc_Ratio,
    Plasticizer_Retained  = Total_Plast_Retained,
    Monomer_Retained      = 1 - mono_flash_frac,
    Foam_Stability        = foam_stability,
    Foam_Drainage_Frac    = foam_drained_frac,
    Bubble_Size_um        = D_b_eff * 1e6,
    Serum_Surface_Tension_N_m = sigma_eff,
    Centrate_Foam_Ratio   = Foam_Potential_Ratio
  ))
}

# =========================================================================
# 1b. CENTRIFUGE -> SPRAY-DRYER HANDOFF
# =========================================================================
# The discharged cake is too concentrated to atomize (see Spray_Visc_Ratio), so
# a reslurry / dilution step sits between the units. This maps every spray-model
# input the centrifuge STREAM determines, after diluting to a target sprayable
# solids mass fraction. Spray-unit operating settings (ALR, air/feed pressures,
# dryer gas flow/temperature/humidity, hold time, Tg, permeabilities, emulsion
# template) are NOT set here - they belong to the spray dryer.
#
# Dilution assumptions (flagged where approximate): solids mass is conserved;
# everything dissolved/entrained in the liquid is diluted by the added water;
# gas-free slurry density is the exact two-density mixture at the target solids;
# bubble size is unchanged; serum surface tension and viscosity move toward the
# clean-water limits as the surfactant/polymer are diluted (first-order).
centrifuge_to_spray <- function(run, target_solid_mass = 0.30) {
  o <- unified_centrifuge_model(run)
  rho_poly <- 1800; rho_liq <- 1000; sigma_clean0 <- 0.072

  Cs_cake <- o$Product_Solids_MassFrac
  target  <- min(target_solid_mass, Cs_cake)               # only dilute, never concentrate
  # liquid dilution factor: added-water ratio to drop solids from cake -> target
  L_over_S_cake   <- (1 - Cs_cake) / max(Cs_cake, 1e-9)
  L_over_S_target <- (1 - target)  / max(target, 1e-9)
  dil <- max(1, L_over_S_target / max(L_over_S_cake, 1e-9))  # >=1 = water added

  rho_L <- 1 / (target / rho_poly + (1 - target) / rho_liq)  # gas-free slurry density

  c(
    rho_L         = rho_L,                                  # slurry density   [kg/m3]
    C_solid_mass  = target,                                 # solids mass frac [-]
    alpha_g_0     = o$Entrained_Gas_Holdup / dil,           # gas holdup (diluted)
    sigma         = o$Serum_Surface_Tension_N_m +
                    (sigma_clean0 - o$Serum_Surface_Tension_N_m) * (1 - 1/dil),  # -> clean on dilution
    D_b           = o$Bubble_Size_um * 1e-6,                # bubble diameter  [m]
    mu_L          = 1e-3 + (o$Serum_Viscosity_Pa_s - 1e-3) / dil,  # serum viscosity [Pa s]
    C_solid_massfrac_cake = Cs_cake,                        # (pre-dilution, for reference)
    C_monomer     = run[["C_monomer"]]     * o$Monomer_Retained / dil,
    C_plasticizer = run[["plasticizer_frac"]] * o$Plasticizer_Retained / dil,
    C_binder      = run[["C_binder"]] / dil,
    I_strength    = run[["I_strength"]] / dil,
    Delta_pH      = run[["Delta_pH"]],
    dilution_x    = dil
  )
}

# =========================================================================
# 2. FACTOR TABLE, SWEEP GROUPS, AND [0,1] -> PHYSICAL MAPPING
# =========================================================================
fac <- function(name, lo, hi, nominal, group, log = FALSE)
  data.frame(name = name, lo = lo, hi = hi, nominal = nominal,
             group = group, log = log, stringsAsFactors = FALSE)

factors <- rbind(
  # --- PROCESS (operator setpoints) ------------------------------------
  fac("rpm",             1500,  4500,  3000,  "process"),   # bowl speed / g-force
  fac("delta_rpm",       1.0,   20.0,  8.0,   "process"),   # scroll differential
  fac("r_pool",          0.118, 0.175, 0.150, "process"),   # weir / pond radius (deeper pond low end)
  fac("flow_rate_lpm",   20,    150,   60,    "process"),
  fac("feed_solid_frac", 0.01,  0.25,  0.15,  "process"),   # down to 1% (heavy feed dilution)
  fac("feed_gas_frac",   0.10,  0.50,  0.25,  "process"),
  fac("pop_frac",        0.20,  0.90,  0.50,  "process"),   # mechanical degassing
  fac("T_process",       288,   343,   303,   "process"),
  fac("P_feed",          1.0e5, 1.0e6, 3.0e5, "process"),
  fac("gas_sat_frac",    0.0,   1.0,   0.5,   "process"),
  fac("D_b",             2.0e-5,2.0e-4,5.0e-5,"process", log = TRUE),  # feed bubble diameter
  # --- SURFACE CHEMISTRY / COLLOIDAL ----------------------------------
  fac("HLB_value",       4.0,   18.0,  12.0,  "surface"),
  fac("surf_dose_kg_m3", 0.5,   4.0,   2.0,   "surface", log = TRUE),
  fac("t_recover_sec",   0.5,   5.0,   2.0,   "surface", log = TRUE),  # foam thixotropic recovery
  fac("I_strength",      1.0e-3,5.0e-1,5.0e-2,"surface", log = TRUE),  # ionic strength (DLVO / salting)
  fac("Delta_pH",        0.2,   4.0,   2.0,   "surface"),   # pH offset from IEP (DLVO)
  fac("contact_angle_deg",30,   150,   90,    "surface"),   # particle wettability (Pickering)
  # --- ADDITIVES ------------------------------------------------------
  fac("plasticizer_frac",0.05,  0.25,  0.12,  "additive"),
  fac("C_monomer",       0.000, 0.020, 0.005, "additive"),
  fac("C_binder",        0.000, 0.050, 0.010, "additive"),
  fac("chi_parameter",   0.1,   0.9,   0.5,   "additive"),
  fac("S_base",          0.10,  0.35,  0.25,  "additive"),   # compacted-sludge floor (vol solids)
  fac("S_ceiling",       0.45,  0.60,  0.50,  "additive"),   # max dewatered (vol solids)

  # --- EQUIPMENT (MUTED by default: hard to change) -------------------
  fac("L_cyl",           0.30,  1.00,  0.60,  "equipment"),
  fac("beach_angle_deg", 3.0,   15.0,  10.0,  "equipment"),
  fac("r_bowl",          0.15,  0.25,  0.20,  "equipment"),
  fac("length_feed",     0.05,  0.20,  0.10,  "equipment")
)

# ---> CHOOSE WHICH GROUPS TO SWEEP HERE <---
# e.g. c("surface") to screen surface chemistry alone, or add "equipment".
active_groups <- c("process", "surface", "additive")
# individual factors to sweep even if their group is muted (e.g. one piece of
# equipment geometry without un-muting all of it)
active_extra  <- c("beach_angle_deg")

factors$active <- (factors$group %in% active_groups) | (factors$name %in% active_extra)
act <- factors[factors$active, ]
params_active <- act$name
k_act <- nrow(act)
nominal_vec <- setNames(factors$nominal, factors$name)

# Build a full physical parameter row from a [0,1] design row over the active
# factors; muted factors stay at their nominal value. Factors flagged log=TRUE
# are mapped geometrically (log-uniform) so a [0,1] step is a fixed number of
# decades - appropriate for wide-range multiplicative factors (e.g. ionic
# strength, bubble size), where the response varies per-decade not per-unit.
build_row <- function(x01) {
  full <- nominal_vec
  x <- x01[act$name]
  vals <- ifelse(act$log,
                 exp(log(act$lo) + x * (log(act$hi) - log(act$lo))),
                 act$lo + x * (act$hi - act$lo))
  full[act$name] <- vals
  full
}

# =========================================================================
# 3. SMOKE TEST + MORRIS SCREENING (sensitivity pkg, else base-R fallback)
# =========================================================================
# --- Smoke test: evaluate the model at all-nominal and print every output ----
smoke <- unlist(unified_centrifuge_model(nominal_vec))
cat("== Smoke test: all factors at nominal ==\n")
print(round(smoke, 5))
if (!all(is.finite(smoke)))
  stop("Smoke test FAILED - non-finite outputs: ",
       paste(names(smoke)[!is.finite(smoke)], collapse = ", "))
cat(sprintf("Smoke test PASSED: %d finite outputs\n\n", length(smoke)))

# --- Base-R Morris OAT design + elementary effects (no package needed) --------
r_traj <- 25L; levels <- 6L; grid_jump <- 3L

build_oat_design <- function(k, r, levels, grid_jump) {
  delta <- grid_jump / (levels - 1)
  grid  <- seq(0, 1 - delta, length.out = max(2L, levels - grid_jump))
  X <- matrix(0, r * (k + 1), k); info <- vector("list", r); row <- 1L
  for (t in seq_len(r)) {
    base <- sample(grid, k, replace = TRUE); ord <- sample.int(k)
    dirs <- numeric(k); x <- base; X[row, ] <- x
    for (s in seq_len(k)) {
      j <- ord[s]
      d <- if (x[j] + delta <= 1) delta else -delta
      x[j] <- x[j] + d; dirs[s] <- d; X[row + s, ] <- x
    }
    info[[t]] <- list(order = ord, dirs = dirs); row <- row + k + 1L
  }
  list(X = X, info = info)
}

elementary_effects <- function(design, Y, k, r) {
  n_out <- ncol(Y)
  ee <- lapply(seq_len(n_out), function(i) matrix(NA_real_, r, k))
  for (t in seq_len(r)) {
    off <- (t - 1) * (k + 1); ord <- design$info[[t]]$order; dirs <- design$info[[t]]$dirs
    for (s in seq_len(k)) {
      j <- ord[s]; dy <- Y[off + s + 1, ] - Y[off + s, ]
      for (i in seq_len(n_out)) ee[[i]][t, j] <- dy[i] / dirs[s]
    }
  }
  ee
}

set.seed(42)
use_pkg <- requireNamespace("sensitivity", quietly = TRUE)

if (use_pkg) {
  message("Using sensitivity::morris()")
  mor  <- sensitivity::morris(model = NULL, factors = params_active, r = r_traj,
                              design = list(type = "oat", levels = levels,
                                            grid.jump = grid_jump),
                              binf = 0, bsup = 1)
  X01  <- mor$X; colnames(X01) <- params_active
  Y <- t(apply(X01, 1, function(x01) unlist(unified_centrifuge_model(build_row(x01)))))
  ee_list <- lapply(seq_len(ncol(Y)), function(i) sensitivity::tell(mor, Y[, i])$ee)
} else {
  message("Package 'sensitivity' not found - using built-in Morris OAT design")
  design <- build_oat_design(k_act, r_traj, levels, grid_jump)
  X01 <- design$X; colnames(X01) <- params_active
  Y <- t(apply(X01, 1, function(x01) unlist(unified_centrifuge_model(build_row(x01)))))
  ee_list <- elementary_effects(design, Y, k_act, r_traj)
}
outputs <- colnames(Y)

dir.create("output", showWarnings = FALSE)

grp_of <- setNames(factors$group, factors$name)
stats_list <- lapply(seq_along(outputs), function(i) {
  ee <- ee_list[[i]]
  data.frame(output  = outputs[i],
             group   = grp_of[params_active],
             factor  = params_active,
             mu      = colMeans(ee),
             mu.star = colMeans(abs(ee)),
             sigma   = apply(ee, 2, sd),
             row.names = NULL)
})
names(stats_list) <- outputs

all_stats <- do.call(rbind, stats_list)
write.csv(all_stats, file.path("output", "centrifuge_morris_indices.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Morris panel grid (mu* vs sigma), one panel per output, coloured by group
# -----------------------------------------------------------------------------
grp_col <- c(process = "steelblue", surface = "darkorange",
             additive = "forestgreen", equipment = "firebrick")

plot_panel <- function(st, title, n_label = 6) {
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25 + 1e-12),
       ylim = c(0, max(st$sigma)   * 1.30 + 1e-12),
       pch = 21, bg = grp_col[st$group], cex = 1.3,
       xlab = expression(mu * "*  (mean |EE|)"),
       ylab = expression(sigma * "  (sd EE)"),
       main = title, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey60")
  top <- order(-st$mu.star)[seq_len(min(n_label, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], labels = st$factor[top],
       pos = 3, cex = 0.68, offset = 0.35, xpd = NA)
}

n_out  <- length(outputs)
ncol_p <- 4L
nrow_p <- ceiling(n_out / ncol_p)
png(file.path("output", "centrifuge_morris_plots.png"),
    width = 850 * ncol_p, height = 750 * nrow_p, res = 150)
op <- par(mfrow = c(nrow_p, ncol_p), mar = c(4.5, 4.5, 2.5, 1),
          oma = c(3.5, 0, 1.5, 0))
for (o in outputs) plot_panel(stats_list[[o]], o)
mtext(sprintf("Centrifuge Morris: groups {%s} active | %d factors, r=25, %d runs",
              paste(active_groups, collapse = ", "), k_act, nrow(Y)),
      side = 1, outer = TRUE, cex = 0.8, line = 0.5)
legend_lab <- names(grp_col)[names(grp_col) %in% act$group]
mtext("", side = 1, outer = TRUE)
par(op)
# group colour legend along the bottom
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("bottom", legend = legend_lab, pt.bg = grp_col[legend_lab],
       pch = 21, horiz = TRUE, bty = "n", cex = 0.9)
dev.off()

# -----------------------------------------------------------------------------
# 5. Console summary
# -----------------------------------------------------------------------------
for (o in outputs) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat("\n==", o, "==\n")
  print(st[, c("group", "factor", "mu.star", "sigma")], row.names = FALSE, digits = 3)
}
cat("\nActive groups:", paste(active_groups, collapse = ", "),
    "|", k_act, "factors\n")
cat("Wrote output/centrifuge_morris_plots.png and output/centrifuge_morris_indices.csv\n")

# -----------------------------------------------------------------------------
# 6. Spray-dryer handoff (nominal case, reslurried to 30% solids)
# -----------------------------------------------------------------------------
spray_ranges <- list(rho_L=c(1000,1300), C_solid_mass=c(0.05,0.40),
                     alpha_g_0=c(0.05,0.60), sigma=c(0.030,0.070),
                     D_b=c(2e-5,2e-4), mu_L=c(0.0012,0.056),
                     C_monomer=c(0,0.02), C_plasticizer=c(0,0.05),
                     C_binder=c(0,0.05), I_strength=c(1e-3,0.5), Delta_pH=c(0.2,4.0))
sf <- centrifuge_to_spray(nominal_vec, target_solid_mass = 0.30)
cat("\n== Centrifuge -> spray-dryer feed (nominal, diluted to 30% solids) ==\n")
for (nm in names(spray_ranges)) {
  rng <- spray_ranges[[nm]]
  flag <- if (sf[[nm]] < rng[1] || sf[[nm]] > rng[2]) "  <-- outside spray range" else ""
  cat(sprintf("  %-14s %11.5g   [%.4g, %.4g]%s\n", nm, sf[[nm]], rng[1], rng[2], flag))
}
cat(sprintf("  (cake was %.1f%% solids; diluted %.1fx to reach the spray feed)\n",
            100 * sf[["C_solid_massfrac_cake"]], sf[["dilution_x"]]))
