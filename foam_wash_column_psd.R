library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN
#   decant pool (hindered settling) + plug foam + bubble population + gas loss
# =====================================================================
# Preformed, particle-loaded foam (no sparger) enters the base. Two zones split
# at the fixed decant depth H_pool. Properties (mu, rho_foam, d_p, d_b_in,
# d_sep) are passed down from the upstream unit.
#
# NEW in this version:
#  (1) HINDERED SETTLING (Richardson-Zaki): the pool Stokes settling velocity
#      is slowed by crowding, U = U_stokes * (1 - phi)^n_RZ, phi = condensed
#      holdup (1 - eps_g). Dense pools settle much slower -> longer decant.
#  (2) BUBBLE POPULATION: mean bubble size d_b(z) evolves by COALESCENCE
#      (grows, faster in dry/unstable foam) minus BREAKAGE (shear, restores
#      toward the inlet size). This is the key operational lever.
#  (3) GAS BALANCE / GAS-LOSS REGIME: no gas is generated in the column, so the
#      gas flux J_g only DECREASES -- bubbles that coarsen past d_b_burst
#      rupture and vent gas. Track J_g in vs out (gas recovery) and the gas
#      lost. Bursting is also what dumps coarse particles (collapse loss is now
#      tied to the burst rate, not an empirical dryness proxy).
#
# States up z: Js_fine, Js_mid, Js_crs [m/s]; C_imp [%]; eps_l [-];
#              d_b [m]; J_g [m/s]; t_res [s].
# PLACEHOLDER kinetics calibrated to: residence ~1.5-2 h, top solids 3-7%.
# =====================================================================

GRAV <- 9.81

params <- c(
  H_total = 5.0, H_pool = 0.45, U_up = 8.0e-4,

  # --- UPSTREAM-PASSED PROPERTIES ---
  mu_cont = 2.0e-3, rho_liquid = 1050, rho_foam = 575, rho_p = 2500, rho_gas = 1.8,
  d_fine = 1.5e-5, d_mid = 8.0e-5, d_crs = 3.0e-4,
  d_b_in = 2.0e-3, d_sep = 4.5e-5, n_RZ = 4.65,

  # --- FEED ---
  eps_l_in = 0.50, eps_s_in = 0.020, eps_g_in = 0.48,
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020, C_imp_in = 100,

  # --- FILM DRAINAGE (time-based) ---
  eps_l_dry = 0.15, tau_drain = 2000, eps_g_pack = 0.74, blend_wz = 0.03,

  # --- BUBBLE POPULATION (coalescence - breakage) ---
  K_coal = 1.3e-4,        # coalescence coeff [1/s] (grows d_b)
  K_break = 1.2e-4,       # breakage coeff [1/s] (restores toward d_b_in)
  # --- GAS LOSS (bursting of over-coarsened bubbles) ---
  d_b_burst = 3.0e-3,     # bubble size above which bursting accelerates [m]
  K_burst = 1.0e-3,       # burst-rate coeff [1/s]

  # --- WASH ---
  wash_ratio = 0.5, k_drainage = 1.8e-3, channel_cv = 0.35, d_b_ref = 2.0e-3,

  # --- PARTICLE LOSS ---
  k_det_fine = 1.2e-5, k_det_mid = 6.0e-5, k_det_crs = 3.0e-4,   # detachment [1/s]
  k_burstloss_fine = 0.05, k_burstloss_mid = 0.30, k_burstloss_crs = 2.0, # burst dumps [-]

  # --- REGIME ---
  film_stability = 1.0, load_sens = 1.0, plug_crit = 1.0
)

stokes  <- function(dRho, d, mu) dRho * GRAV * d^2 / (18 * mu)
rho_agg <- function(dp, db, rp, rg) (rp*dp^3 + rg*db^3) / (dp^3 + db^3)

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    load_rel <- (Js_fine + Js_mid + Js_crs) / (Js_fine_in + Js_mid_in + Js_crs_in)
    eps_s    <- eps_s_in * load_rel
    eps_g    <- 1 - eps_l - eps_s
    U <- U_up
    dt_res <- 1 / U
    foam_frac <- 1 / (1 + exp(-(z - H_pool) / blend_wz))

    # (2) BUBBLE POPULATION: coalescence (dry/unstable -> faster) - breakage
    coal_rate <- K_coal * (eps_g / eps_g_pack) / film_stability     # [1/s]
    dd_b <- (coal_rate * d_b - K_break * (d_b - d_b_in)) / U

    # (3) GAS LOSS: bursting when bubbles coarsen past d_b_burst (gas-only sink)
    burst_rate <- K_burst * max(0, d_b - d_b_burst) / d_b_in         # [1/s]
    dJ_g <- -burst_rate * J_g / U

    # film drainage (foam dries -> eps_l falls toward eps_l_dry)
    deps_l <- -(eps_l - eps_l_dry) / tau_drain / U

    # (1) POOL decant with HINDERED settling (Richardson-Zaki)
    phi <- 1 - eps_g                                                 # condensed holdup
    U_settle_sep <- stokes(rho_liquid - rho_foam, d_sep, mu_cont) * (1 - phi)^n_RZ
    dCimp_pool <- -(U_settle_sep / H_pool) * C_imp

    # POOL particle settling (aggregate; uses current d_b)
    settle_loss <- function(dp) {
      d_agg <- (dp^3 + d_b^3)^(1/3)
      Uset  <- stokes(rho_agg(dp, d_b, rho_p, rho_gas) - rho_liquid, d_agg, mu_cont) * (1 - phi)^n_RZ
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

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l, dd_b, dJ_g, dt_res))
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
          J_g=p[["eps_g_in"]]*p[["U_up"]], t_res=0)
  z <- seq(0, p[["H_total"]], by = 0.02)
  out <- as.data.frame(ode(s0, z, column_model_psd, p)); names(out)[1] <- "z"
  lr <- (out$Js_fine+out$Js_mid+out$Js_crs)/(p[["Js_fine_in"]]+p[["Js_mid_in"]]+p[["Js_crs_in"]])
  out$eps_s <- p[["eps_s_in"]]*lr
  out$eps_g <- 1 - out$eps_l - out$eps_s
  out$solids_pct <- 100*out$eps_s/(out$eps_s+out$eps_l)
  out$regime <- vapply(seq_len(nrow(out)), function(i) regime_of(out$z[i], lr[i], p), character(1))
  out
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
t_pool <- approx(out_psd$z, out_psd$t_res, params[["H_pool"]])$y
t_tot  <- out_psd$t_res[nrow(out_psd)]
Jg_in  <- params[["eps_g_in"]]*params[["U_up"]]; Jg_out <- out_psd$J_g[nrow(out_psd)]
top <- out_psd[nrow(out_psd),]; ini <- out_psd[1,]

cat(sprintf("Bubble size d_b: inlet %.2f mm -> top %.2f mm (coalescence)\n",
            params[["d_b_in"]]*1e3, top$d_b*1e3))
cat(sprintf("GAS balance: in %.2e -> out %.2e m/s  => recovery %.0f%%, lost %.0f%%\n",
            Jg_in, Jg_out, 100*Jg_out/Jg_in, 100*(1-Jg_out/Jg_in)))
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
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(out_psd$z, out_psd$Js_mid, type="l", col="green", lwd=3, ylim=c(0,0.05),
     ylab="Solid flux (m/s)", xlab="Height (m)", main="PSD retention")
lines(out_psd$z, out_psd$Js_fine, col="blue", lwd=2, lty=2)
lines(out_psd$z, out_psd$Js_crs,  col="red",  lwd=2, lty=3)
abline(v=params[["H_pool"]], col="gray", lty=2)
legend("right", legend=c("Mid","Fine","Coarse"), col=c("green","blue","red"), lty=c(1,2,3), lwd=2)

plot(out_psd$z, out_psd$d_b*1e3, type="l", col="darkorange", lwd=3,
     ylab="Bubble d_b (mm)", xlab="Height (m)", main="Bubble growth & burst")
abline(h=params[["d_b_burst"]]*1e3, col="red", lty=3)
legend("topleft", legend=c("d_b","burst onset"), col=c("darkorange","red"), lty=c(1,3), lwd=2)

plot(out_psd$z, out_psd$J_g/Jg_in*100, type="l", col="purple", lwd=3, ylim=c(0,100),
     ylab="Gas remaining (% of inlet)", xlab="Height (m)", main="Gas-loss regime")
abline(v=params[["H_pool"]], col="gray", lty=2)
