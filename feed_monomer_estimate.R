# =============================================================================
# Free gaseous-monomer estimate from the upstream UP1 feed density
# =============================================================================
# up1_rho is measured on the stream ENTERING UP1 (colloid + gaseous monomer +
# liquid, mostly water), UPSTREAM of the foaming air (up1_scfh, added inside UP1).
# Monomer is gaseous at room temperature (MW ~ 70 g/mol); as free gas it lowers
# up1_rho below the colloid(1700) + water baseline. That deficit converts to
# monomer mass by the ideal-gas law:
#
#   rho_theo = 1 / (s/rho_s + (1-s)/rho_w)          # colloid + water, gas-free
#   rho_gas  = P * MW / (R * T)                      # ~2.86 kg/m3 at 1 atm, 25 C
#   x_mono   = (1/rho_meas - 1/rho_theo) / (1/rho_gas - 1/rho_w)   # mass frac
#
# CAVEATS (see DATA_REVIEW.md / review/up3_density_monomer_note.md):
#   * Bounds FREE-gas monomer only. Dissolved monomer occupies liquid-like
#     volume, leaves no PV=nRT signature, and needs TGA/GC.
#   * Hypersensitive: 1 kg/m3 of deficit ~ 2.3 ppm monomer (1 atm, 25 C), so the
#     sign is set by the colloid+water baseline (feed T -> rho_w, dissolved
#     surfactant/solids, exact solid fraction). Pin the baseline with a MEASURED
#     monomer-free feed density on the same rig to turn the bound into a signed
#     number. With current data the estimate is ~ -15 ppm, i.e. < ~30 ppm.
#
# Run:  Rscript feed_monomer_estimate.R   ->  data/feed_monomer_estimate.csv
# =============================================================================

R_GAS  <- 8.314        # J/mol/K
MW_MONO <- 0.070       # kg/mol (~70 g/mol, confirmed)
ATM    <- 101325       # Pa
RHO_S  <- 1700         # kg/m3, bare polymer colloid (= dried-product skeletal)

# Water density [kg/m3] vs temperature [C] (Kell fit, 0-100 C, 1 atm).
rho_water <- function(T_C) {
  (999.83952 + 16.945176*T_C - 7.9870401e-3*T_C^2 - 46.170461e-6*T_C^3 +
     105.56302e-9*T_C^4 - 280.54253e-12*T_C^5) / (1 + 16.879850e-3*T_C)
}

# Free gaseous-monomer mass fraction [ppm of feed] from a feed density.
#   rho_meas   measured up1_rho [kg/m3]
#   s          solid (colloid) mass fraction     (default 0.25)
#   feed_T_C   feed temperature [C]              (default 25; ambient, upstream)
#   feed_P_atm feed pressure [atm]               (default 1)
free_monomer_ppm <- function(rho_meas, s = 0.25, feed_T_C = 25, feed_P_atm = 1,
                             rho_s = RHO_S) {
  T_K     <- feed_T_C + 273.15
  rho_w   <- rho_water(feed_T_C)
  rho_gas <- feed_P_atm * ATM * MW_MONO / (R_GAS * T_K)
  rho_theo <- 1 / (s/rho_s + (1 - s)/rho_w)
  x_mono   <- (1/rho_meas - 1/rho_theo) / (1/rho_gas - 1/rho_w)
  data.frame(colloid_water_kgm3 = rho_theo,
             deficit_kgm3        = rho_theo - rho_meas,
             free_monomer_ppm    = x_mono * 1e6)
}

# --- apply to the 4 clean chain feeds + historical steady feeds --------------
if (sys.nframe() == 0) {
  feed_T_C <- 25; feed_P_atm <- 1   # <-- set to true upstream feed conditions

  chain <- read.csv("data/cond_up1234.csv")
  ch <- cbind(source = "chain", cond = chain$cond,
              up1_rho_kgm3 = round(chain$up1_feed_rho_kgm3, 1),
              feed_T_C = feed_T_C, feed_P_atm = feed_P_atm,
              round(free_monomer_ppm(chain$up1_feed_rho_kgm3, 0.25,
                                     feed_T_C, feed_P_atm), 1))

  out <- ch
  # Historical feed rows (up1_2_3_solids_rho.xlsx) if readxl is available:
  if (requireNamespace("readxl", quietly = TRUE)) {
    h <- readxl::read_excel("data/up1_2_3_solids_rho.xlsx", sheet = "Sheet1")
    rho <- suppressWarnings(as.numeric(h$up1feed_rho))
    keep <- !is.na(rho) & rho > 1080          # drop startup/dilution transients
    hh <- cbind(source = "historical", cond = "",
                up1_rho_kgm3 = round(rho[keep], 1),
                feed_T_C = feed_T_C, feed_P_atm = feed_P_atm,
                round(free_monomer_ppm(rho[keep], 0.25, feed_T_C, feed_P_atm), 1))
    out <- rbind(ch, hh)
  }

  write.csv(out, "data/feed_monomer_estimate.csv", row.names = FALSE)
  cat(sprintf("Wrote data/feed_monomer_estimate.csv (%d rows)\n", nrow(out)))
  cat(sprintf("Sensitivity: 1 kg/m3 deficit ~ %.2f ppm (%.0f C, %.0f atm)\n",
              free_monomer_ppm(1000, 0.25, feed_T_C, feed_P_atm)$free_monomer_ppm -
              free_monomer_ppm(1001, 0.25, feed_T_C, feed_P_atm)$free_monomer_ppm,
              feed_T_C, feed_P_atm))
  cat("chain-condition estimates:\n"); print(ch)
}
