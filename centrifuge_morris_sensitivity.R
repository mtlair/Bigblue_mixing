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
# NOTE on surface chemistry (Section F): the surfactant monolayer partition is
# grounded (coverage = MW/(a0 N_A), geometric interfacial areas, dose balance
# vs CMC), but it currently only produces Centrate_Foam_Ratio - it does NOT
# feed back into the foam cushion, degassing or aggregate cohesion. Dose is now
# a swept surface factor; coupling it back is a future refinement.
#
# Method: Morris OAT screening on a [0,1] hypercube over the ACTIVE factors
# (comparable elementary effects), mapped to physical ranges, all outputs from
# one design. CSV + mu* vs sigma panel grid written to output/.
#
# Requires the `sensitivity` package.  Run: Rscript centrifuge_morris_sensitivity.R
# =============================================================================

# install.packages("sensitivity")
library(sensitivity)

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
  # surface chemistry
  HLB_value        <- run[["HLB_value"]]
  surf_dose_kg_m3  <- run[["surf_dose_kg_m3"]]   # surfactant dose          [kg/m3]
  t_recover_sec    <- run[["t_recover_sec"]]     # thixotropic foam recovery[s]
  # additives
  plasticizer_frac <- run[["plasticizer_frac"]]
  C_monomer        <- run[["C_monomer"]]
  C_binder         <- run[["C_binder"]]          # binder concentration     [wt/wt]
  chi_parameter    <- run[["chi_parameter"]]
  # equipment (muted by default)
  L_cyl            <- run[["L_cyl"]]             # cylinder length          [m]
  r_bowl           <- run[["r_bowl"]]            # bowl radius (diameter)   [m]
  beach_angle_deg  <- run[["beach_angle_deg"]]
  length_feed      <- run[["length_feed"]]       # feed-zone length         [m]

  # ---------------------------------------------------------
  # B. FIXED CONSTANTS & DERIVED GEOMETRY
  # ---------------------------------------------------------
  r_discharge <- 0.120        # solids discharge radius (fixed) [m]
  scroll_pitch <- 0.050       # scroll pitch (fixed)            [m]
  beach_angle <- beach_angle_deg * pi / 180
  # keep the pond physical even if equipment geometry is swept
  r_pool <- min(max(r_pool, r_discharge + 0.005), r_bowl - 0.005)

  rho_poly <- 1200; rho_liq <- 1000
  R_gas <- 8.314; R_air <- 287
  P_atm <- 101325; omega <- rpm * (2 * pi / 60)

  visc_water <- mu_water_T(T_process)

  # Core/template density coupled to plasticizer (theta-solvent core)
  rho_plast_theta <- 820
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

  rho_gas <- P_atm / (R_air * T_process)

  Q_gas_vented <- Q_gas_atm * pop_frac
  Q_gas_surviving_atm <- Q_gas_atm - Q_gas_vented

  P_hydro <- 0.5 * rho_liq * omega^2 * (r_bowl^2 - r_pool^2)
  P_local <- P_atm + P_hydro
  Q_gas_local <- Q_gas_surviving_atm * (P_atm / P_local)

  Q_local_m3s <- Q_solid + Q_liq + Q_gas_local
  rho_feed_local <- ((Q_solid * rho_poly) + (Q_liq * rho_liq) + (Q_gas_surviving_atm * rho_gas)) / Q_local_m3s

  phi_dispersed <- (Q_solid + Q_gas_local) / Q_local_m3s
  safe_phi <- min(phi_dispersed, 0.64)
  mu_apparent <- visc_water * (1 - (safe_phi / 0.65))^(-2.5 * 0.65)

  # Zone-separated residence
  A_pond <- pi * (r_bowl^2 - r_pool^2)
  U_ax   <- Q_local_m3s / A_pond
  t_feed    <- length_feed / U_ax
  t_clarify <- L_cyl / U_ax
  t_pond    <- t_feed + t_clarify

  foam_dampening_full <- max(1.0, mu_apparent / visc_water)
  recovery_frac <- 1 - exp(-t_pond / t_recover_sec)
  foam_dampening_factor <- 1.0 + (foam_dampening_full - 1.0) * recovery_frac

  k_feed_accel <- 80
  n_rev_feed   <- omega * t_feed / (2 * pi)
  feed_shear_factor <- 1 + k_feed_accel / max(n_rev_feed, 1e-6)

  # ---------------------------------------------------------
  # C.5 THERMODYNAMIC FLASHING (T, P dependent)
  # ---------------------------------------------------------
  kH <- k_H_T(T_process)
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
  size_bins <- seq(10, 150, length.out = 15)
  core_bins <- seq(0.1, 0.6, length.out = 15)
  n_size <- length(size_bins)

  weight_sum <- sum(dnorm(size_bins, 50, 20) %o% dnorm(core_bins, 0.3, 0.1))

  heavy_solid_yield <- 0
  light_solid_yield <- 0
  intact_template_yield <- 0
  agg_survival <- 0
  size_capture   <- numeric(n_size)
  size_feed_mass <- numeric(n_size)

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
  t_gap <- ((r_pool - r_discharge) / tan(beach_angle)) / v_convey

  # Binder holds interstitial moisture in the cake
  cake_moisture_frac <- 0.15 + 0.35 * exp(-t_gap / 5.0) + 0.10 * plasticizer_frac + 0.30 * C_binder

  gas_frac_local   <- Q_gas_local / Q_local_m3s
  cake_gas_sparged <- (0.02 + 0.20 * exp(-t_gap / 3.0)) * (0.5 + 2.5 * gas_frac_local)
  cake_gas_flash   <- V_flash_per_liq * cake_moisture_frac + 0.5 * C_monomer * mono_flash_frac
  cake_gas_frac    <- min(0.5, cake_gas_sparged + cake_gas_flash)

  cake_solid_frac <- max(1e-3, 1.0 - cake_moisture_frac - cake_gas_frac)

  # Herschel-Bulkley wet-cake rheology (yield-stress paste / mousse); binder
  # bridging stiffens the frictional network.
  phi_pack <- 0.64; p_exp <- 3.0
  tau_y0   <- 5.0e4
  sigma_lg <- 0.030; R_bub <- 50e-6
  n_hb <- 0.5; K_hb <- 200; gamma_disc <- 10

  phi_ratio  <- min(cake_solid_frac, 0.98 * phi_pack) / phi_pack
  tau_y_fric <- tau_y0 * phi_ratio^p_exp * (1 + 5.0 * C_binder)
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
  # F. SURFACTANT PARTITIONING & SECONDARY FOAMING (surface group)
  # ---------------------------------------------------------
  # Monolayer partition: coverage = MW / (a0 N_A); surfactant splits across the
  # polymer and bubble interfaces (HLB-weighted) with the remainder dissolved
  # in the centrate, whose concentration vs CMC gives the foam potential.
  # NOTE: currently a terminal output - not fed back into foam/cohesion.
  CMC_kg_m3 <- 1.0; MW <- 500; a0_angstrom2 <- 50
  kg_per_m2 <- (MW / 1000) / ((a0_angstrom2 * 1e-20) * 6.022e23)

  SA_poly <- (6 / 50e-6) * Q_solid
  SA_gas  <- (6 / 50e-6) * Q_gas_surviving_atm

  Surf_on_Gas <- SA_gas * kg_per_m2
  affinity_factor <- 1 / (1 + exp((HLB_value - 10) / 2))
  Surf_on_Poly <- SA_poly * kg_per_m2 * affinity_factor

  Total_Surf_In <- Q_liq * surf_dose_kg_m3
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
# 2. FACTOR TABLE, SWEEP GROUPS, AND [0,1] -> PHYSICAL MAPPING
# =========================================================================
fac <- function(name, lo, hi, nominal, group)
  data.frame(name = name, lo = lo, hi = hi, nominal = nominal,
             group = group, stringsAsFactors = FALSE)

factors <- rbind(
  # --- PROCESS (operator setpoints) ------------------------------------
  fac("rpm",             1500,  4500,  3000,  "process"),   # bowl speed / g-force
  fac("delta_rpm",       1.0,   20.0,  8.0,   "process"),   # scroll differential
  fac("r_pool",          0.130, 0.170, 0.150, "process"),   # weir / pond radius
  fac("flow_rate_lpm",   20,    150,   60,    "process"),
  fac("feed_solid_frac", 0.05,  0.25,  0.15,  "process"),
  fac("feed_gas_frac",   0.10,  0.50,  0.25,  "process"),
  fac("pop_frac",        0.20,  0.90,  0.50,  "process"),   # mechanical degassing
  fac("T_process",       288,   343,   303,   "process"),
  fac("P_feed",          1.0e5, 1.0e6, 3.0e5, "process"),
  fac("gas_sat_frac",    0.0,   1.0,   0.5,   "process"),
  # --- SURFACE CHEMISTRY ----------------------------------------------
  fac("HLB_value",       4.0,   18.0,  12.0,  "surface"),
  fac("surf_dose_kg_m3", 0.5,   4.0,   2.0,   "surface"),
  fac("t_recover_sec",   0.5,   5.0,   2.0,   "surface"),   # foam thixotropic recovery
  # --- ADDITIVES ------------------------------------------------------
  fac("plasticizer_frac",0.05,  0.25,  0.12,  "additive"),
  fac("C_monomer",       0.000, 0.020, 0.005, "additive"),
  fac("C_binder",        0.000, 0.050, 0.010, "additive"),
  fac("chi_parameter",   0.1,   0.9,   0.5,   "additive"),
  # --- EQUIPMENT (MUTED by default: hard to change) -------------------
  fac("L_cyl",           0.30,  1.00,  0.60,  "equipment"),
  fac("beach_angle_deg", 5.0,   15.0,  10.0,  "equipment"),
  fac("r_bowl",          0.15,  0.25,  0.20,  "equipment"),
  fac("length_feed",     0.05,  0.20,  0.10,  "equipment")
)

# ---> CHOOSE WHICH GROUPS TO SWEEP HERE <---
# e.g. c("surface") to screen surface chemistry alone, or add "equipment".
active_groups <- c("process", "surface", "additive")

factors$active <- factors$group %in% active_groups
act <- factors[factors$active, ]
params_active <- act$name
k_act <- nrow(act)
nominal_vec <- setNames(factors$nominal, factors$name)

# Build a full physical parameter row from a [0,1] design row over the active
# factors; muted factors stay at their nominal value.
build_row <- function(x01) {
  full <- nominal_vec
  full[act$name] <- act$lo + x01[act$name] * (act$hi - act$lo)
  full
}

# =========================================================================
# 3. EXECUTE MORRIS SCREENING FOR ALL OUTPUTS
# =========================================================================
set.seed(42)

mor <- morris(model = NULL, factors = params_active, r = 25,
              design = list(type = "oat", levels = 6, grid.jump = 3),
              binf = 0, bsup = 1)

Y <- t(apply(mor$X, 1, function(x01) unlist(unified_centrifuge_model(build_row(x01)))))
outputs <- colnames(Y)

dir.create("output", showWarnings = FALSE)

grp_of <- setNames(factors$group, factors$name)
stats_list <- lapply(outputs, function(o) {
  m  <- sensitivity::tell(mor, Y[, o])
  ee <- m$ee
  data.frame(output  = o,
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
