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
# Module chain: Feed mass balance / compression -> foam cushion (Krieger-
# Dougherty + thixotropic recovery) -> particle tracking (shear rupture +
# core-shell migration) -> heavy-cake rheology & plasticizer retention ->
# surfactant partitioning & secondary foaming.
#
# Aggregate integrity: the primary solids aggregate in the first formulation
# step (aggregate size grows with plasticizer_frac); Aggregate_Destruction
# reports the distribution-weighted fraction of those aggregates ruptured by
# centrifugal shear. The core/template density is coupled to plasticizer_frac
# as a theta-solvent phase (lighter than bulk liquid), which sets the density
# contrast abs(rho_poly - template_density) driving that shear - so it is
# derived here, not swept as an independent factor.
#
# Requires the `sensitivity` package (install.packages("sensitivity")).
#
# Run:  Rscript centrifuge_morris_sensitivity.R
# =============================================================================

# install.packages("sensitivity")
library(sensitivity)

# =========================================================================
# 1. THE UNIFIED COMPREHENSIVE CENTRIFUGE ENGINE
# =========================================================================
unified_centrifuge_model <- function(run) {

  # ---------------------------------------------------------
  # A. EXTRACT SWEPT PARAMETERS
  # ---------------------------------------------------------
  # NOTE: use [[ ]] so each value comes out as a plain unnamed scalar;
  # single-bracket [ ] keeps the factor name attached and propagates a
  # stale name onto every downstream result.
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
  # template_density is NOT swept independently - it is derived from
  # plasticizer_frac below (theta-solvent coupling).

  # ---------------------------------------------------------
  # B. FIXED EQUIPMENT GEOMETRY & CONSTANTS
  # ---------------------------------------------------------
  r_bowl <- 0.200; r_pool <- 0.150; r_discharge <- 0.120
  length_cyl <- 0.600; length_pond <- 0.500; length_feed <- 0.100
  beach_angle <- 10 * pi/180; scroll_pitch <- 0.050

  rho_poly <- 1200; rho_liq <- 1000; rho_gas <- 1.2
  visc_water <- 0.001; tau_max <- 25.0
  P_atm <- 101325; omega <- rpm * (2 * pi / 60)

  # Core/template density COUPLED to plasticizer content (theta-solvent core).
  # More plasticizer -> a more solvent-rich, lighter core that sits below the
  # 1000 kg/m3 bulk liquid. Linear from bulk-liquid density at the low-
  # plasticizer limit down to a theta-solvent density at the high-plasticizer
  # limit. A lighter core widens abs(rho_poly - template_density), the
  # centrifugal contrast that stresses (and can rupture) the aggregates.
  rho_plast_theta <- 820           # theta-solvent core density (max plasticizer) [kg/m3]
  plast_lo <- 0.05; plast_hi <- 0.25
  plast_norm <- max(0, min(1, (plasticizer_frac - plast_lo) / (plast_hi - plast_lo)))
  template_density <- rho_liq - (rho_liq - rho_plast_theta) * plast_norm

  # ---------------------------------------------------------
  # C. MASS BALANCE, COMPRESSION, & DEGASSING
  # ---------------------------------------------------------
  Q_in_m3s <- flow_rate / 60000
  # guard: never let the rounding of the fractions drive the liquid share
  # negative (safe under the current bounds, cheap insurance against edits)
  feed_liq_frac <- max(0, 1.0 - (feed_solid_frac + feed_gas_frac))

  Q_solid <- Q_in_m3s * feed_solid_frac
  Q_liq   <- Q_in_m3s * feed_liq_frac
  Q_gas_atm <- Q_in_m3s * feed_gas_frac

  # Venting (The Pop) & Hydrostatic Compression
  Q_gas_vented <- Q_gas_atm * pop_frac
  Q_gas_surviving_atm <- Q_gas_atm - Q_gas_vented

  P_hydro <- 0.5 * rho_liq * omega^2 * (r_bowl^2 - r_pool^2)
  P_local <- P_atm + P_hydro
  Q_gas_local <- Q_gas_surviving_atm * (P_atm / P_local)

  # Local flow and properties inside the pond
  Q_local_m3s <- Q_solid + Q_liq + Q_gas_local
  rho_feed_local <- ((Q_solid * rho_poly) + (Q_liq * rho_liq) + (Q_gas_surviving_atm * rho_gas)) / Q_local_m3s

  # Dynamic Foam Cushion (Krieger-Dougherty Apparent Viscosity)
  phi_dispersed <- (Q_solid + Q_gas_local) / Q_local_m3s
  safe_phi <- min(phi_dispersed, 0.64)
  mu_apparent <- visc_water * (1 - (safe_phi / 0.65))^(-2.5 * 0.65)

  t_pond <- (pi * (r_bowl^2 - r_pool^2) * length_pond) / Q_local_m3s

  # Thixotropic foam recovery: the foam cushion is sheared down at the feed
  # zone and rebuilds its structure over the characteristic recovery time
  # t_recover_sec. Only the fraction of structure that recovers within the
  # pond residence time t_pond contributes to the protective dampening, so a
  # short residence / slow-recovering foam offers less shear protection.
  foam_dampening_full <- max(1.0, mu_apparent / visc_water)
  recovery_frac <- 1 - exp(-t_pond / t_recover_sec)
  foam_dampening_factor <- 1.0 + (foam_dampening_full - 1.0) * recovery_frac

  # ---------------------------------------------------------
  # D. PARTICLE TRACKING: SHEAR RUPTURE & MIGRATION
  # ---------------------------------------------------------
  size_bins <- seq(10, 150, length.out = 15)
  core_bins <- seq(0.1, 0.6, length.out = 15) # Template volume fraction

  # common normalization denominator (total distribution weight); used for
  # ALL three yields so they share one [0,1] scale
  weight_sum <- sum(dnorm(size_bins, 50, 20) %o% dnorm(core_bins, 0.3, 0.1))

  heavy_solid_yield <- 0
  light_solid_yield <- 0
  intact_template_yield <- 0
  agg_survival <- 0   # distribution-weighted fraction surviving shear intact

  for (d_um in size_bins) {
    for (core_frac in core_bins) {

      weight <- dnorm(d_um, 50, 20) * dnorm(core_frac, 0.3, 0.1)
      if(is.na(weight) || weight <= 0) next

      # Aggregate Structure & Stickiness
      d_primary_m <- d_um * 1e-6
      agg_multiplier <- 1 + (10 * plasticizer_frac)
      d_intact_m <- d_primary_m * agg_multiplier

      # 1. Shear Rupture Check (Protected by foam)
      cohesive_strength <- 1000 + (50000 * plasticizer_frac)
      shear_stress <- (abs(rho_poly - template_density) * r_pool * omega^2 * d_intact_m) / foam_dampening_factor

      survival_prob <- 1 / (1 + exp(-(cohesive_strength - shear_stress) / 500))

      # How much of the first-step aggregate survives the centrifuge shear
      agg_survival <- agg_survival + (survival_prob * weight)

      # Determine migration physics based on survival
      # Path A: It Survives (Migrates as intact core-shell)
      rho_intact <- (core_frac * template_density) + ((1 - core_frac) * rho_poly)
      v_c_intact <- (d_intact_m^2 * (rho_intact - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent)

      # Path B: It Ruptures (Template lost, core collapses to bare polymer)
      d_ruptured <- d_primary_m * (agg_multiplier * 0.5)
      v_c_rupt <- (d_ruptured^2 * (rho_poly - rho_feed_local) * r_bowl * omega^2) / (18 * mu_apparent)

      # Calculate recoveries
      dist <- r_bowl - r_pool
      rec_intact_heavy <- ifelse(v_c_intact > 0, 1 - exp(-abs(v_c_intact) * t_pond / dist), 0)
      rec_intact_light <- ifelse(v_c_intact < 0, 1 - exp(-abs(v_c_intact) * t_pond / dist), 0)
      rec_rupt_heavy   <- ifelse(v_c_rupt > 0,   1 - exp(-abs(v_c_rupt) * t_pond / dist), 0)

      # Accumulate outputs (raw weighted sums; normalized once after the loop)
      intact_template_yield <- intact_template_yield + ((rec_intact_heavy + rec_intact_light) * survival_prob * weight)
      heavy_solid_yield <- heavy_solid_yield + (rec_intact_heavy * survival_prob * weight) + (rec_rupt_heavy * (1 - survival_prob) * weight)
      light_solid_yield <- light_solid_yield + (rec_intact_light * survival_prob * weight)
    }
  }

  # Normalize ALL yields on the common denominator so heavy, light and intact
  # are directly comparable fractions and the 3-phase centrate balance closes.
  intact_template_yield <- intact_template_yield / weight_sum
  heavy_solid_yield     <- heavy_solid_yield / weight_sum
  light_solid_yield     <- light_solid_yield / weight_sum
  agg_survival          <- agg_survival / weight_sum

  # ---------------------------------------------------------
  # E. HEAVY CAKE RHEOLOGY & PLASTICIZER RETENTION
  # ---------------------------------------------------------
  # 1. Drainage & Paste Viscosity
  v_convey <- (delta_rpm / 60) * scroll_pitch
  t_gap <- ((r_pool - r_discharge) / tan(beach_angle)) / v_convey

  # Cake moisture: base drainage decay in conveying time, plus extra held
  # moisture from a tacky / plasticized cake (plasticizer raises cohesion and
  # capillary retention, so the cake drains less completely).
  cake_moisture_frac <- 0.15 + 0.35 * exp(-t_gap / 5.0) + 0.10 * plasticizer_frac

  # Cake gas: base bridging decay scaled by how much gas actually survived to
  # the pond (compressed local gas fraction) rather than a fixed amount, so
  # the cake gas content is coupled to feed gas, popping and bowl speed.
  gas_frac_local <- Q_gas_local / Q_local_m3s
  cake_gas_frac  <- (0.02 + 0.20 * exp(-t_gap / 3.0)) * (0.5 + 2.5 * gas_frac_local)
  cake_gas_frac  <- min(cake_gas_frac, 0.35)

  cake_solid_frac <- max(1e-3, 1.0 - cake_moisture_frac - cake_gas_frac)

  safe_solid_phi <- min(cake_solid_frac, 0.64)
  mu_friction <- visc_water * (1 - (safe_solid_phi / 0.65))^(-2.5 * 0.65)
  mu_bridging <- 1.0 + (50.0 * cake_gas_frac)
  paste_viscosity_pa_s <- mu_friction * mu_bridging

  # 2. Plasticizer Retention (Flory-Huggins)
  Q_liq_safe <- max(Q_liq, 1e-12)
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
    Wet_Cake_Moisture     = cake_moisture_frac,
    Paste_Viscosity_Pa_s  = paste_viscosity_pa_s,
    Plasticizer_Retained  = Total_Plast_Retained,
    Centrate_Foam_Ratio   = Foam_Potential_Ratio
  ))
}

# =========================================================================
# 2. THE MORRIS WRAPPER FUNCTION
# =========================================================================
# Morris requires a single numeric output. This wrapper lets you choose
# which output from the unified list to analyze.

run_morris_wrapper <- function(X, target_variable = "Intact_Template_Yield") {
  results <- apply(X, 1, function(row) {
    full_output <- unified_centrifuge_model(row)
    return(full_output[[target_variable]])
  })
  return(results)
}

# =========================================================================
# 3. SETUP & EXECUTE MORRIS SENSITIVITY
# =========================================================================
# template_density is not listed here: it is derived inside the model from
# plasticizer_frac (theta-solvent coupling), so it is not an independent knob.
params <- c("rpm", "delta_rpm", "flow_rate_lpm", "feed_solid_frac",
            "feed_gas_frac", "pop_frac", "plasticizer_frac",
            "chi_parameter", "HLB_value", "t_recover_sec")

# Define Boundaries for your Process
lower <- c(1500,  1.0,  20, 0.05, 0.10, 0.20, 0.05, 0.1,  4.0, 0.5)
upper <- c(4500, 20.0, 150, 0.25, 0.50, 0.90, 0.25, 0.9, 18.0, 5.0)

set.seed(42)

# ---> CHANGE THE TARGET VARIABLE HERE <---
# Options: Intact_Template_Yield, Aggregate_Survival, Aggregate_Destruction,
# Heavy_Solid_Yield, Yield_Loss_Centrate, Wet_Cake_Moisture,
# Paste_Viscosity_Pa_s, Plasticizer_Retained, Centrate_Foam_Ratio

# Default screens what destroys the first-step aggregates in the centrifuge.
target <- "Aggregate_Destruction"

morris_out <- morris(model = function(X) run_morris_wrapper(X, target),
                     factors = params, r = 25,
                     design = list(type = "oat", levels = 6, grid.jump = 3),
                     binf = lower, bsup = upper)

print(morris_out)
plot(morris_out, main = paste("Morris Sensitivity:", target))
