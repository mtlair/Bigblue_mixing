# =============================================================================
# Morris (Elementary Effects) sensitivity analysis of a unified decanter-
# centrifuge model for a three-phase polymer separation.
#
# Process context: this is a THREE-PHASE separation.
#   * Light phase  - foam / low-density solids that report to the overflow
#                    side but are RECOVERED (combined downstream with heavy)
#   * Heavy phase  - settled solids forming the cake                  RECOVERED
#   * Centrate     - the middle liquid phase                          LOST
# Hence Yield_Loss_Centrate = 1 - (Heavy_Solid_Yield + Light_Solid_Yield),
# with both yields on a common normalized [0,1] scale.
#
# ZONE-SEPARATED RESIDENCE (three decanter zones):
#   * Feed / acceleration zone (length_feed): the slurry is spun up to bowl
#     speed. High shear, no net clarification. Aggregates accelerated within
#     fewer revolutions feel a sharper spin-up shear (feed_shear_factor), so
#     short feed dwell / high throughput ruptures more of them.
#   * Cylindrical clarification zone (L_cyl, SWEPT): the main settling length.
#     Radial migration is captured here over t_clarify ~ L_cyl / U_axial, so a
#     LONGER CYLINDER gives more residence and higher fine-particle capture.
#   * Conical beach / dewatering zone: cake conveyed up the beach, drained over
#     t_gap (set by the scroll differential delta_rpm and the beach geometry).
#
# OUTLET PARTICLE DISTRIBUTION: the per-size capture is accumulated across the
# size bins to give a grade-efficiency curve, reported as the classification
# cut size d50c (size at 50 % capture) and the cake solids median d50.
#
# Carrier: the pond continuous phase is water (feed carrier plus upstream
# dilution water), so the pond viscosity is water-based but TEMPERATURE
# dependent (Vogel). The dewatered product is a wet cake / wet mousse, a
# yield-stress material - Herschel-Bulkley closure, not a Newtonian viscosity.
#
# Thermodynamic layer (T and P explicit unit-op state):
#   * Dissolved-gas flash template: gas loaded at the feed pressure P_feed
#     (Henry's law, T-dependent) comes out of solution on the pond -> discharge
#     letdown, nucleating template voids in the cake - in parallel with the
#     mechanically sparged / entrained gas held as a pore former.
#   * Residual monomer volatility (Clausius-Clapeyron): sets flashed vs
#     retained monomer; retained monomer plasticizes / softens the cake.
#   * Temperature drives the carrier viscosity and the cake softening.
#
# Aggregate integrity: primary solids aggregate in the first formulation step
# (size grows with plasticizer_frac); Aggregate_Destruction is the
# distribution-weighted fraction ruptured by centrifugal shear. The core /
# template density is coupled to plasticizer_frac as a theta-solvent phase
# (lighter than bulk liquid), setting the density contrast that drives shear.
#
# Method: Morris OAT screening on a [0,1] hypercube (so elementary effects are
# comparable across factors of very different physical scale), mapped to
# physical ranges, evaluated for ALL outputs from one design. CSV + a mu* vs
# sigma panel grid are written to output/.
#
# Requires the `sensitivity` package (install.packages("sensitivity")).
# Run:  Rscript centrifuge_morris_sensitivity.R
# =============================================================================

# install.packages("sensitivity")
library(sensitivity)

# ---- Temperature-dependent physical property helpers ------------------------
# Water dynamic viscosity vs temperature (Vogel correlation), [Pa s]
mu_water_T <- function(T) 2.414e-5 * 10^(247.8 / (T - 140))
# Henry solubility constant vs temperature (N2-like, van't Hoff), [mol/(m3 Pa)]
# Higher at lower T (gas more soluble cold); scale is gas-identity dependent.
k_H_T <- function(T) 6.5e-6 * exp(1700 * (1/T - 1/298))

# =========================================================================
# 1. THE UNIFIED COMPREHENSIVE CENTRIFUGE ENGINE
# =========================================================================
unified_centrifuge_model <- function(run) {

  # ---------------------------------------------------------
  # A. EXTRACT SWEPT PARAMETERS
  # ---------------------------------------------------------
  # [[ ]] returns a plain unnamed scalar; [ ] would keep the factor name and
  # propagate a stale name onto every downstream result.
  rpm               <- run[["rpm"]]
  delta_rpm         <- run[["delta_rpm"]]
  flow_rate         <- run[["flow_rate_lpm"]]
  feed_solid_frac   <- run[["feed_solid_frac"]]
  feed_gas_frac     <- run[["feed_gas_frac"]]
  pop_frac          <- run[["pop_frac"]]          # Mechanical degassing efficiency
  plasticizer_frac  <- run[["plasticizer_frac"]]  # Polymer tackiness/cohesion
  chi_parameter     <- run[["chi_parameter"]]     # Flory-Huggins compatibility
  HLB_value         <- run[["HLB_value"]]         # Surfactant affinity
  t_recover_sec     <- run[["t_recover_sec"]]     # Thixotropic foam recovery time
  T_process         <- run[["T_process"]]         # Process temperature          [K]
  P_feed            <- run[["P_feed"]]            # Feed-line (gas-loading) pressure [Pa]
  gas_sat_frac      <- run[["gas_sat_frac"]]      # Dissolved-gas saturation      [-]
  C_monomer         <- run[["C_monomer"]]         # Residual monomer conc.    [wt/wt]
  L_cyl             <- run[["L_cyl"]]             # Cylindrical clarification length [m]

  # ---------------------------------------------------------
  # B. FIXED EQUIPMENT GEOMETRY & CONSTANTS
  # ---------------------------------------------------------
  r_bowl <- 0.200; r_pool <- 0.150; r_discharge <- 0.120
  length_feed <- 0.100                            # feed / acceleration zone length [m]
  beach_angle <- 10 * pi/180; scroll_pitch <- 0.050

  rho_poly <- 1200; rho_liq <- 1000
  R_gas <- 8.314; R_air <- 287
  P_atm <- 101325; omega <- rpm * (2 * pi / 60)

  # Carrier (water) viscosity at the process temperature.
  visc_water <- mu_water_T(T_process)

  # Core/template density COUPLED to plasticizer content (theta-solvent core).
  rho_plast_theta <- 820           # theta-solvent core density (max plasticizer) [kg/m3]
  plast_lo <- 0.05; plast_hi <- 0.25
  plast_norm <- max(0, min(1, (plasticizer_frac - plast_lo) / (plast_hi - plast_lo)))
  template_density <- rho_liq - (rho_liq - rho_plast_theta) * plast_norm

  # ---------------------------------------------------------
  # C. MASS BALANCE, COMPRESSION, & DEGASSING
  # ---------------------------------------------------------
  Q_in_m3s <- flow_rate / 60000
  feed_liq_frac <- max(0, 1.0 - (feed_solid_frac + feed_gas_frac))

  Q_solid <- Q_in_m3s * feed_solid_frac
  Q_liq   <- Q_in_m3s * feed_liq_frac
  Q_gas_atm <- Q_in_m3s * feed_gas_frac

  rho_gas <- P_atm / (R_air * T_process)          # entrained gas density at T, 1 atm

  Q_gas_vented <- Q_gas_atm * pop_frac
  Q_gas_surviving_atm <- Q_gas_atm - Q_gas_vented

  P_hydro <- 0.5 * rho_liq * omega^2 * (r_bowl^2 - r_pool^2)
  P_local <- P_atm + P_hydro
  Q_gas_local <- Q_gas_surviving_atm * (P_atm / P_local)

  Q_local_m3s <- Q_solid + Q_liq + Q_gas_local
  rho_feed_local <- ((Q_solid * rho_poly) + (Q_liq * rho_liq) + (Q_gas_surviving_atm * rho_gas)) / Q_local_m3s

  # Dynamic Foam Cushion (Krieger-Dougherty, water base at T)
  phi_dispersed <- (Q_solid + Q_gas_local) / Q_local_m3s
  safe_phi <- min(phi_dispersed, 0.64)
  mu_apparent <- visc_water * (1 - (safe_phi / 0.65))^(-2.5 * 0.65)

  # --- Zone-separated residence times -----------------------------------
  A_pond <- pi * (r_bowl^2 - r_pool^2)            # annular pond cross-section [m2]
  U_ax   <- Q_local_m3s / A_pond                  # axial pond velocity        [m/s]
  t_feed    <- length_feed / U_ax                 # feed / acceleration dwell  [s]
  t_clarify <- L_cyl / U_ax                        # cylindrical settling dwell [s]
  t_pond    <- t_feed + t_clarify                  # total pond residence       [s]

  # Thixotropic foam recovery over the full pond residence
  foam_dampening_full <- max(1.0, mu_apparent / visc_water)
  recovery_frac <- 1 - exp(-t_pond / t_recover_sec)
  foam_dampening_factor <- 1.0 + (foam_dampening_full - 1.0) * recovery_frac

  # Feed-zone acceleration: aggregates spun up within fewer revolutions feel a
  # sharper shear. k_feed_accel is a lumped acceleration-severity calibration.
  k_feed_accel <- 80
  n_rev_feed   <- omega * t_feed / (2 * pi)        # revolutions during feed dwell
  feed_shear_factor <- 1 + k_feed_accel / max(n_rev_feed, 1e-6)

  # ---------------------------------------------------------
  # C.5 THERMODYNAMIC FLASHING (T, P dependent)
  # ---------------------------------------------------------
  # Dissolved-gas template released on the pond -> discharge letdown.
  kH <- k_H_T(T_process)                          # mol/(m3 Pa)
  C_gas_loaded <- gas_sat_frac * kH * P_feed       # dissolved gas at feed  [mol/m3]
  C_gas_final  <- kH * P_atm                        # solubility at discharge [mol/m3]
  flash_mol_per_m3 <- max(0, C_gas_loaded - C_gas_final)
  V_flash_per_liq  <- flash_mol_per_m3 * R_gas * T_process / P_atm  # m3 gas / m3 liq

  # Residual monomer volatility (Clausius-Clapeyron).
  Hvap_mono <- 32000; T_bp_mono <- 325            # moderately volatile monomer
  p_mono <- P_atm * exp((Hvap_mono / R_gas) * (1 / T_bp_mono - 1 / T_process))
  mono_flash_frac <- 1 / (1 + exp(-(p_mono - P_atm) / (0.2 * P_atm)))
  C_mono_retained <- C_monomer * (1 - mono_flash_frac)

  # ---------------------------------------------------------
  # D. PARTICLE TRACKING: SHEAR RUPTURE, MIGRATION & CLASSIFICATION
  # ---------------------------------------------------------
  size_bins <- seq(10, 150, length.out = 15)      # primary diameter [um]
  core_bins <- seq(0.1, 0.6, length.out = 15)     # template (core) volume fraction
  n_size <- length(size_bins)

  weight_sum <- sum(dnorm(size_bins, 50, 20) %o% dnorm(core_bins, 0.3, 0.1))

  heavy_solid_yield <- 0
  light_solid_yield <- 0
  intact_template_yield <- 0
  agg_survival <- 0
  size_capture   <- numeric(n_size)               # heavy (cake) mass per size
  size_feed_mass <- numeric(n_size)               # feed mass per size

  # Hindered settling (Richardson-Zaki, Stokes regime)
  phi_s_local <- Q_solid / Q_local_m3s
  hinder <- (max(0, 1 - phi_s_local))^4.65
  dist <- r_bowl - r_pool

  for (i in seq_len(n_size)) {
    d_um <- size_bins[i]
    for (core_frac in core_bins) {

      weight <- dnorm(d_um, 50, 20) * dnorm(core_frac, 0.3, 0.1)
      if(is.na(weight) || weight <= 0) next

      d_primary_m <- d_um * 1e-6
      agg_multiplier <- 1 + (10 * plasticizer_frac)
      d_intact_m <- d_primary_m * agg_multiplier

      # Shear rupture, intensified by the feed-zone acceleration
      cohesive_strength <- 1000 + (50000 * plasticizer_frac)
      shear_stress <- (abs(rho_poly - template_density) * r_pool * omega^2 * d_intact_m) *
                      feed_shear_factor / foam_dampening_factor
      survival_prob <- 1 / (1 + exp(-(cohesive_strength - shear_stress) / 500))
      agg_survival <- agg_survival + (survival_prob * weight)

      # Path A: survives (intact core-shell); Path B: ruptures (bare polymer)
      rho_intact <- (core_frac * template_density) + ((1 - core_frac) * rho_poly)
      v_c_intact <- (d_intact_m^2 * (rho_intact - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent) * hinder
      d_ruptured <- d_primary_m * (agg_multiplier * 0.5)
      v_c_rupt   <- (d_ruptured^2 * (rho_poly - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent) * hinder

      # Recoveries over the CLARIFICATION-zone residence (settling happens in
      # the cylinder, not the feed / acceleration zone)
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

  # Normalize yields on the common denominator (shared [0,1] scale)
  intact_template_yield <- intact_template_yield / weight_sum
  heavy_solid_yield     <- heavy_solid_yield / weight_sum
  light_solid_yield     <- light_solid_yield / weight_sum
  agg_survival          <- agg_survival / weight_sum

  # --- Outlet classification: grade efficiency -> cut size + cake median ---
  grade_eff <- size_capture / pmax(size_feed_mass, 1e-12)   # T(d): fraction to cake
  feed_mean_size <- sum(size_bins * size_feed_mass) / sum(size_feed_mass)
  d50c <- tryCatch(approx(grade_eff, size_bins, xout = 0.5,
                          ties = "ordered", rule = 2)$y,
                   error = function(e) NA_real_)
  if (is.na(d50c)) d50c <- feed_mean_size
  if (sum(size_capture) > 1e-12) {
    cake_cdf <- cumsum(size_capture) / sum(size_capture)
    cake_d50 <- tryCatch(approx(cake_cdf, size_bins, xout = 0.5,
                                ties = "ordered", rule = 2)$y,
                         error = function(e) feed_mean_size)
  } else {
    cake_d50 <- feed_mean_size
  }

  # ---------------------------------------------------------
  # E. WET-CAKE / MOUSSE RHEOLOGY & RETENTION
  # ---------------------------------------------------------
  # Beach / dewatering zone residence
  v_convey <- (delta_rpm / 60) * scroll_pitch
  t_gap <- ((r_pool - r_discharge) / tan(beach_angle)) / v_convey

  cake_moisture_frac <- 0.15 + 0.35 * exp(-t_gap / 5.0) + 0.10 * plasticizer_frac

  # Cake gas = sparged / entrained holdup + dissolved-gas flash template
  gas_frac_local   <- Q_gas_local / Q_local_m3s
  cake_gas_sparged <- (0.02 + 0.20 * exp(-t_gap / 3.0)) * (0.5 + 2.5 * gas_frac_local)
  cake_gas_flash   <- V_flash_per_liq * cake_moisture_frac + 0.5 * C_monomer * mono_flash_frac
  cake_gas_frac    <- min(0.5, cake_gas_sparged + cake_gas_flash)

  cake_solid_frac <- max(1e-3, 1.0 - cake_moisture_frac - cake_gas_frac)

  # Herschel-Bulkley wet-cake rheology (yield-stress paste / mousse)
  phi_pack <- 0.64; p_exp <- 3.0
  tau_y0   <- 5.0e4                # frictional yield scale            [Pa]
  sigma_lg <- 0.030; R_bub <- 50e-6  # foam film tension / cell radius
  n_hb <- 0.5; K_hb <- 200; gamma_disc <- 10   # HB flow index / consistency / discharge shear

  phi_ratio  <- min(cake_solid_frac, 0.98 * phi_pack) / phi_pack
  tau_y_fric <- tau_y0 * phi_ratio^p_exp
  tau_y_foam <- (sigma_lg / R_bub) * max(cake_gas_frac - 0.10, 0)
  soften     <- exp(-2.0 * (plasticizer_frac + C_mono_retained)) * exp(-0.02 * (T_process - 293))
  cake_yield_stress <- (tau_y_fric + tau_y_foam) * soften

  paste_viscosity_pa_s <- cake_yield_stress / gamma_disc + K_hb * gamma_disc^(n_hb - 1)

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
  # F. SURFACTANT PARTITIONING & SECONDARY FOAMING
  # ---------------------------------------------------------
  surf_feed_kg_m3 <- 2.0; CMC_kg_m3 <- 1.0; MW <- 500; a0_angstrom2 <- 50
  kg_per_m2 <- (MW / 1000) / ((a0_angstrom2 * 1e-20) * 6.022e23)

  SA_poly <- (6 / 50e-6) * Q_solid
  SA_gas  <- (6 / 50e-6) * Q_gas_surviving_atm

  Surf_on_Gas <- SA_gas * kg_per_m2
  affinity_factor <- 1 / (1 + exp((HLB_value - 10) / 2))
  Surf_on_Poly <- SA_poly * kg_per_m2 * affinity_factor

  Total_Surf_In <- Q_liq * surf_feed_kg_m3
  Surf_Dissolved <- max(0, Total_Surf_In - Surf_on_Gas - Surf_on_Poly)

  Q_liq_centrate <- Q_liq - Q_liq_heavy - (Q_gas_surviving_atm * 0.05)
  Centrate_Conc_kg_m3 <- Surf_Dissolved / max(Q_liq_centrate, 1e-6)
  Foam_Potential_Ratio <- Centrate_Conc_kg_m3 / CMC_kg_m3

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
    Gas_Template_Voidage  = cake_gas_frac,
    Cake_Yield_Stress_Pa  = cake_yield_stress,
    Paste_Viscosity_Pa_s  = paste_viscosity_pa_s,
    Plasticizer_Retained  = Total_Plast_Retained,
    Monomer_Retained      = 1 - mono_flash_frac,
    Centrate_Foam_Ratio   = Foam_Potential_Ratio
  ))
}

# =========================================================================
# 2. FACTORS, RANGES, AND [0,1] -> PHYSICAL MAPPING
# =========================================================================
# template_density is NOT a factor: it is derived from plasticizer_frac.
params <- c("rpm", "delta_rpm", "flow_rate_lpm", "feed_solid_frac",
            "feed_gas_frac", "pop_frac", "plasticizer_frac",
            "chi_parameter", "HLB_value", "t_recover_sec",
            "T_process", "P_feed", "gas_sat_frac", "C_monomer", "L_cyl")

#          rpm  drpm  flow  fsol  fgas  pop  plast  chi  HLB  trec  T_K    P_feed  gsat  Cmono  Lcyl
lower <- c(1500,  1.0,  20, 0.05, 0.10, 0.20, 0.05, 0.1,  4.0, 0.5,  288, 1.0e5,  0.0, 0.000, 0.30)
upper <- c(4500, 20.0, 150, 0.25, 0.50, 0.90, 0.25, 0.9, 18.0, 5.0,  343, 1.0e6,  1.0, 0.020, 1.00)
k <- length(params)

# Map a unit-hypercube design to physical values (linear per factor)
scale_design <- function(X01) {
  X <- X01
  for (j in seq_len(k)) X[, j] <- lower[j] + X01[, j] * (upper[j] - lower[j])
  colnames(X) <- params
  X
}

# =========================================================================
# 3. EXECUTE MORRIS SCREENING FOR ALL OUTPUTS
# =========================================================================
set.seed(42)

mor <- morris(model = NULL, factors = params, r = 25,
              design = list(type = "oat", levels = 6, grid.jump = 3),
              binf = 0, bsup = 1)

Xphys <- scale_design(mor$X)
Y <- t(apply(Xphys, 1, function(row) unlist(unified_centrifuge_model(row))))
outputs <- colnames(Y)

dir.create("output", showWarnings = FALSE)

stats_list <- lapply(outputs, function(o) {
  m  <- sensitivity::tell(mor, Y[, o])
  ee <- m$ee                                       # r x k elementary effects
  data.frame(output  = o,
             factor  = params,
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
# 4. Morris panel grid (mu* vs sigma), one panel per output
# -----------------------------------------------------------------------------
plot_panel <- function(st, title, n_label = 6) {
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25 + 1e-12),
       ylim = c(0, max(st$sigma)   * 1.30 + 1e-12),
       pch = 21, bg = "steelblue", cex = 1.2,
       xlab = expression(mu * "*  (mean |EE|)"),
       ylab = expression(sigma * "  (sd EE)"),
       main = title, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey60")            # sigma = mu* (non-linear / interacting)
  top <- order(-st$mu.star)[seq_len(min(n_label, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], labels = st$factor[top],
       pos = 3, cex = 0.70, offset = 0.35, xpd = NA)
}

n_out  <- length(outputs)
ncol_p <- 4L
nrow_p <- ceiling(n_out / ncol_p)
png(file.path("output", "centrifuge_morris_plots.png"),
    width = 850 * ncol_p, height = 750 * nrow_p, res = 150)
op <- par(mfrow = c(nrow_p, ncol_p), mar = c(4.5, 4.5, 2.5, 1),
          oma = c(2.5, 0, 1.5, 0))
for (o in outputs) plot_panel(stats_list[[o]], o)
mtext(sprintf("Centrifuge Morris screening: %d factors, r=25, 6 levels, %d runs | dashed: sigma = mu*",
              k, nrow(Y)),
      side = 1, outer = TRUE, cex = 0.8, line = 1)
par(op); dev.off()

# -----------------------------------------------------------------------------
# 5. Console summary
# -----------------------------------------------------------------------------
for (o in outputs) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat("\n==", o, "==\n")
  print(st[, c("factor", "mu.star", "sigma")], row.names = FALSE, digits = 3)
}
cat("\nWrote output/centrifuge_morris_plots.png and output/centrifuge_morris_indices.csv\n")
