library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN (decant pool + plug-flow foam)
# =====================================================================
# Preformed, particle-loaded foam (no sparger) enters the base as a wet
# dispersion. Two zones, split at the fixed decant depth H_pool:
#
#   BOTTOM  z < H_pool : GRAVITY-DECANT POOL. The pool VELOCITY is a proper
#           density-differential SETTLING law (Stokes), computed from the
#           passed-down dispersed size d_sep, the phase density difference
#           (rho_liquid - rho_foam), and the continuous-phase viscosity mu:
#               U_settle = (rho_liquid - rho_foam) * g * d_sep^2 / (18 mu)
#           This is slow when the density differential is low -> the ~30 min
#           decant timescale. No more free tau_sep / pool_mobility.
#   TOP     z >= H_pool: THROUGHPUT plug-flow foam at the (very low) feed
#           velocity U_up, washed by liquid added at the top.
#
# THROUGHPUT + GATE: residence is throughput-based (H / U_up per zone). The
# decant COMPLETES only if the clarifying settling keeps up with the upflow,
# i.e. completion = min(1, U_settle / U_up). Very low upward velocity -> good
# separation; too-fast feed -> impurity carryover. This replaces "pool
# residence forced = decant time" with the physical competition of the two.
#
# Particle loss: coarse aggregates (particle + attached bubble) settle out of
# the pool only if their net density beats the liquid (Stokes on the aggregate
# with bubble size d_b); with normal bubbles all classes are buoyant, so
# coarse loss is the foam-zone detachment/collapse below. The pool settling
# term is included and self-activates if the bubbles are small enough.
#
# Passed down from the UPSTREAM unit (inputs, not fitted): mu_cont, rho_foam,
# d_fine/d_mid/d_crs, d_b, d_sep.
#
# States up z: Js_fine, Js_mid, Js_crs [m/s]; C_imp [%]; eps_l [-]; t_res [s].
# PLACEHOLDER foam loss/wash constants calibrated to: residence 1.5-2.5 h and
# top solids 3-7%.
# =====================================================================

GRAV <- 9.81

params <- c(
  H_total = 5.0, H_pool = 0.45,          # column & decant depth [m] (~1.5 ft)
  U_up = 8.0e-4,                         # throughput superficial upward velocity [m/s]

  # --- UPSTREAM-PASSED PHYSICAL PROPERTIES (inputs, not fitted) ---
  mu_cont   = 2.0e-3,                    # continuous-phase viscosity [Pa s] (=2 cP)
  rho_liquid = 1050, rho_foam = 575, rho_p = 2500, rho_gas = 1.8,  # [kg/m3]
  d_fine = 1.5e-5, d_mid = 8.0e-5, d_crs = 3.0e-4,  # particle diameters [m]
  d_b    = 2.0e-3,                       # bubble diameter (set upstream by P) [m]
  d_sep  = 4.5e-5,                       # decant-controlling dispersed size [m]

  # --- FEED LOADING ---
  eps_l_in = 0.50, eps_s_in = 0.020,
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020, C_imp_in = 100,

  # --- FILM DRAINAGE (foam dries toward a wet equilibrium) ---
  eps_l_dry = 0.15, L_drain = 1.5, eps_g_pack = 0.74, blend_wz = 0.03,
  d_b_ref = 2.0e-3,                      # bubble-size reference for channeling ratio

  # --- WASH LIQUID ---
  wash_ratio = 0.5,

  # --- FOAM-ZONE PARTICLE LOSS (per time): detachment + collapse ---
  k_det_fine = 1.2e-5, k_det_mid = 6.0e-5, k_det_crs = 3.0e-4,   # [1/s]
  k_col_fine = 5.0e-5, k_col_mid = 2.5e-4, k_col_crs = 1.5e-3,   # [1/s]

  # --- CHANNELING / WASH DISTRIBUTION ---
  k_drainage = 1.8e-3, channel_cv = 0.35,

  # --- REGIME (snow-globe vs plug): material stability vs loading ---
  film_stability = 1.0, load_sens = 1.0, plug_crit = 1.0
)

# Stokes settling velocity [m/s]; positive = settles, negative = rises
stokes <- function(dRho, d, mu) dRho * GRAV * d^2 / (18 * mu)
# net density of a particle+bubble aggregate (volume-weighted spheres)
rho_agg <- function(dp, db, rp, rg) (rp*dp^3 + rg*db^3) / (dp^3 + db^3)

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    load_rel <- (Js_fine + Js_mid + Js_crs) / (Js_fine_in + Js_mid_in + Js_crs_in)
    eps_s    <- eps_s_in * load_rel
    eps_g    <- 1 - eps_l

    U <- U_up                                   # throughput velocity (very low)
    dt_res <- 1 / U

    foam_frac <- 1 / (1 + exp(-(z - H_pool) / blend_wz))   # 0 pool -> 1 foam
    deps_l <- -(eps_l - eps_l_dry) / L_drain               # film drainage

    # ---- POOL: gravity decant of impurity (density-differential settling) ----
    U_settle_sep <- stokes(rho_liquid - rho_foam, d_sep, mu_cont)
    # removal rate per height = (settling velocity / depth) / U ; over the pool
    # residence H_pool/U this integrates to exp(-U_settle_sep/U) -> the
    # completion = min(1, U_settle/U_up) competition.
    dCimp_pool <- -(U_settle_sep / H_pool) * C_imp

    # ---- POOL: particle settling (aggregate; coarse only if it beats liquid)
    settle_loss <- function(dp) {
      d_agg <- (dp^3 + d_b^3)^(1/3)
      Uset  <- stokes(rho_agg(dp, d_b, rho_p, rho_gas) - rho_liquid, d_agg, mu_cont)
      max(0, Uset) / H_pool                      # [1/s]; >0 only if aggregate sinks
    }
    sl_fine <- settle_loss(d_fine); sl_mid <- settle_loss(d_mid); sl_crs <- settle_loss(d_crs)

    # ---- FOAM: detachment + collapse (attached-particle losses) ----
    coalescence <- max(0, eps_g - eps_g_pack)
    lf_fine <- k_det_fine + k_col_fine * coalescence
    lf_mid  <- k_det_mid  + k_col_mid  * coalescence
    lf_crs  <- k_det_crs  + k_col_crs  * coalescence

    # blend pool settling and foam detachment; convert per-time -> per-height
    dJs_fine <- -Js_fine * ((1-foam_frac)*sl_fine + foam_frac*lf_fine) / U
    dJs_mid  <- -Js_mid  * ((1-foam_frac)*sl_mid  + foam_frac*lf_mid ) / U
    dJs_crs  <- -Js_crs  * ((1-foam_frac)*sl_crs  + foam_frac*lf_crs ) / U

    # ---- FOAM: washing, penalized by Plateau-border channeling ----
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    wash_strength   <- wash_ratio / (wash_ratio + 1)
    dCimp_foam <- -k_drainage * Wash_Efficiency * wash_strength * C_imp

    dCimp <- ((1-foam_frac)*dCimp_pool + foam_frac*dCimp_foam) / U

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l, dt_res))
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
          C_imp=p[["C_imp_in"]], eps_l=p[["eps_l_in"]], t_res=0)
  z <- seq(0, p[["H_total"]], by = 0.02)
  out <- as.data.frame(ode(s0, z, column_model_psd, p)); names(out)[1] <- "z"
  out$eps_g <- 1 - out$eps_l
  lr <- (out$Js_fine+out$Js_mid+out$Js_crs)/(p[["Js_fine_in"]]+p[["Js_mid_in"]]+p[["Js_crs_in"]])
  out$eps_s <- p[["eps_s_in"]]*lr
  out$solids_pct <- 100*out$eps_s/(out$eps_s+out$eps_l)
  out$regime <- vapply(seq_len(nrow(out)), function(i) regime_of(out$z[i], lr[i], p), character(1))
  out
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
U_set  <- stokes(params[["rho_liquid"]]-params[["rho_foam"]], params[["d_sep"]], params[["mu_cont"]])
completion <- min(1, U_set / params[["U_up"]])
t_pool <- approx(out_psd$z, out_psd$t_res, params[["H_pool"]])$y
t_tot  <- out_psd$t_res[nrow(out_psd)]
top <- out_psd[nrow(out_psd),]; ini <- out_psd[1,]

cat(sprintf("Upward velocity U_up   = %.2e m/s  |  settling U_settle = %.2e m/s\n",
            params[["U_up"]], U_set))
cat(sprintf("Decant completion      = %.0f%%  (min(1, U_settle/U_up))\n", 100*completion))
cat(sprintf("Decant time  H/U_settle= %.2f h  [Stokes]   vs heuristic 100mu/dRho= %.2f h\n",
            params[["H_pool"]]/U_set/3600,
            100*(params[["mu_cont"]]*1000)/(params[["rho_liquid"]]-params[["rho_foam"]])))
cat(sprintf("Residence: pool %.2f h + foam %.2f h = %.2f h (throughput)\n",
            t_pool/3600, (t_tot-t_pool)/3600, t_tot/3600))
cat(sprintf("Solids content: inlet %.1f%% -> top %.1f%% | plug %.0f%%\n",
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

plot(out_psd$z, out_psd$solids_pct, type="l", col="darkorange", lwd=3, ylim=c(0,50),
     ylab="%", xlab="Height (m)", main="Solids content & impurity")
lines(out_psd$z, out_psd$C_imp/2, col="purple", lwd=2, lty=2)
abline(v=params[["H_pool"]], col="gray", lty=2)
legend("topright", legend=c("Solids %","Impurity %/2"), col=c("darkorange","purple"), lty=c(1,2), lwd=2)

plot(out_psd$z, out_psd$t_res/3600, type="l", col="black", lwd=3,
     ylab="Cumulative residence (h)", xlab="Height (m)", main="Residence (throughput)")
abline(v=params[["H_pool"]], col="gray", lty=2)
plug_z <- out_psd$z[out_psd$regime=="plug"]; if (length(plug_z)) rug(plug_z, col="forestgreen", lwd=2)
legend("topleft", legend=c("t_res","interface","plug"), col=c("black","gray","forestgreen"), lty=c(1,2,1), lwd=2)
