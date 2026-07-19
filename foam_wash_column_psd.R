library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN (decant pool + plug-flow foam)
# =====================================================================
# The foam is PREFORMED and PARTICLE-LOADED upstream (no gas sparger). It
# enters the base as a WET dispersion. The column has two zones set by the
# local gas fraction eps_g = 1 - eps_l:
#
#   BOTTOM  z < H_pool : GRAVITY-DECANT LIQUID POOL. Residence is set by the
#           decanter heuristic  t_decant[h] = 100*mu/(rho_A - rho_B)  (mu in cP,
#           rho in kg/m3) -- typically ~0.5 h, rising toward ~1.5 h as the
#           density differential shrinks. The interface (decant depth H_pool)
#           is treated as FIXED, since it is observed to be stable in plug flow.
#   TOP     z >= H_pool: drained PLUG-FLOW FOAM, washed by liquid added at the
#           top. At slow feed this zone adds ~1 h, so total residence ~1.5 h.
#
# TIME SCALE: total residence is HOURS. Because residence >> fast film
# drainage, low top-of-column solids (~3-7%) is NOT "drainage ran out of
# time" -- the foam drains to a WET equilibrium holdup eps_l_eq (capillary /
# particle-held), which is simply wet. The slow steps are the gravity decant
# (heuristic above) and the slow plug-flow feed.
#
# States marched up z: Js_fine, Js_mid, Js_crs [m/s]; C_imp [%];
#                      eps_l (liquid holdup) [-]; t_res (cumulative time) [s].
#
# Snow-globe vs plug is still a MATERIAL property (film stability vs loading).
# PLACEHOLDER closures throughout -- constants calibrated to two anchors:
#   total residence 1.5-2.5 h  and  top solids content 3-7%.
# =====================================================================

params <- c(
  H_total = 5.0,

  # --- PREFORMED FOAM FEED (enters as a wet dispersion) ---
  J_foam   = 8.3e-4,     # superficial foam velocity in the plug zone [m/s]
  eps_l_in = 0.50,       # inlet liquid holdup (wet dispersion) [-]
  eps_s_in = 0.020,      # inlet solid volume fraction [-] (for solids %)
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020,  # inlet loading [m/s]
  C_imp_in   = 100,      # inlet entrained impurity [%]

  # --- GRAVITY DECANT (bottom zone) : heuristic  t[h] = 100*mu/(rho_A-rho_B) ---
  H_pool = 1.5,          # decant zone depth [m] (stable interface in plug flow)
  mu_cont = 2.0,         # continuous-phase viscosity [cP]
  rho_A   = 1050,        # heavier separating liquid [kg/m3]
  rho_B   = 650,         # lighter separating liquid [kg/m3]

  # --- PREFORMING PRESSURE -> INLET BUBBLE SIZE ---
  P = 150000, P_ref = 150000, d_b_ref = 2.0e-3,

  # --- FILM DRAINAGE (foam dries toward a WET equilibrium in the plug zone) ---
  eps_l_dry_ref = 0.15,  # fully-drained (wet) residual holdup at P_ref [-]
  L_drain       = 1.5,   # film-drainage length in the foam [m]
  eps_g_pack    = 0.74,  # packing gas fraction (collapse / regime threshold)
  blend_wz      = 0.10,  # smooth pool<->foam velocity switch width [m]

  # --- WASH LIQUID (added at top, correlated to the foam feed) ---
  wash_ratio = 0.5,      # wash-water : entrained-feed-liquid ratio

  # --- PSD LOSS: buoyancy-failure detachment + bubble-collapse (per unit TIME)
  k_det_fine = 1.2e-5, k_det_mid = 6.0e-5, k_det_crs = 3.0e-4,   # [1/s]
  k_col_fine = 5.0e-5, k_col_mid = 2.5e-4, k_col_crs = 1.5e-3,   # [1/s]

  # --- CHANNELING / WASH DISTRIBUTION ---
  k_drainage = 1.8e-3,   # base impurity sweep rate [1/s]
  channel_cv = 0.35,     # CV of Plateau-border widths (0 = ideal plug wash)

  # --- REGIME (snow-globe vs plug): material stability vs loading ---
  film_stability = 1.0, load_sens = 1.0, plug_crit = 1.0
)

# NOTE: loss/wash rates are now per-unit-TIME [1/s]; the model converts to
# per-height with the local velocity (d/dz = (1/U) d/dt), so the hours-long
# residence -- not an arbitrary length -- sets how much is lost/washed.
column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    Js_tot_in <- Js_fine_in + Js_mid_in + Js_crs_in
    load_rel  <- (Js_fine + Js_mid + Js_crs) / Js_tot_in
    eps_s     <- eps_s_in * load_rel

    d_b      <- d_b_ref * (P_ref / P)^(1 / 3)
    eps_l_eq <- eps_l_dry_ref * (d_b_ref / d_b)

    eps_g <- 1 - eps_l

    # Zone velocities. Pool residence is set by the decanter heuristic
    #   t_decant[h] = 100 * mu / (rho_A - rho_B);  U_pool = H_pool / t_decant.
    # The interface is a fixed decant depth (stable in plug flow), so the
    # switch is geometric on z, not on holdup.
    t_decant_s <- (100 * mu_cont / max(rho_A - rho_B, 1)) * 3600
    U_pool <- H_pool / t_decant_s
    U_foam <- J_foam / max(eps_g, 1e-3)                          # feed-driven plug
    foam_frac <- 1 / (1 + exp(-(z - H_pool) / blend_wz))         # 0 pool -> 1 foam
    U <- max((1 - foam_frac) * U_pool + foam_frac * U_foam, 1e-7)

    # cumulative residence time
    dt_res <- 1 / U

    # film drainage toward the (wet) equilibrium over a length in the foam
    deps_l <- -(eps_l - eps_l_eq) / L_drain

    # size loss (per time -> per height via 1/U): detach + dryness-gated collapse
    coalescence <- max(0, eps_g - eps_g_pack)
    dJs_fine <- -Js_fine * (k_det_fine + k_col_fine * coalescence) / U
    dJs_mid  <- -Js_mid  * (k_det_mid  + k_col_mid  * coalescence) / U
    dJs_crs  <- -Js_crs  * (k_det_crs  + k_col_crs  * coalescence) / U

    # washing (foam zone only): swept by wash liquid, penalized by channeling
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    wash_strength   <- wash_ratio / (wash_ratio + 1)
    dCimp <- -foam_frac * k_drainage * Wash_Efficiency * wash_strength * C_imp / U

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l, dt_res))
  })
}

# --- helpers ---------------------------------------------------------------
d_b_of_P    <- function(p) p[["d_b_ref"]] * (p[["P_ref"]] / p[["P"]])^(1/3)
t_decant_h  <- function(p) 100 * p[["mu_cont"]] / max(p[["rho_A"]] - p[["rho_B"]], 1)
# regime: pool below the (fixed) decant depth; above it, plug vs snow-globe
# from film stability vs solid loading.
regime_of <- function(z, load_rel, p) {
  if (z < p[["H_pool"]]) return("pool")
  plug_index <- p[["film_stability"]] / (p[["load_sens"]] * load_rel)
  if (plug_index >= p[["plug_crit"]]) "plug" else "snowglobe"
}

solve_column <- function(p) {
  state_init <- c(Js_fine = p[["Js_fine_in"]], Js_mid = p[["Js_mid_in"]],
                  Js_crs = p[["Js_crs_in"]], C_imp = p[["C_imp_in"]],
                  eps_l = p[["eps_l_in"]], t_res = 0)
  z_seq <- seq(0, p[["H_total"]], by = 0.02)
  out <- as.data.frame(ode(y = state_init, times = z_seq,
                           func = column_model_psd, parms = p))
  names(out)[1] <- "z"
  out$eps_g <- 1 - out$eps_l
  Js_tot_in <- p[["Js_fine_in"]] + p[["Js_mid_in"]] + p[["Js_crs_in"]]
  load_rel  <- (out$Js_fine + out$Js_mid + out$Js_crs) / Js_tot_in
  out$eps_s <- p[["eps_s_in"]] * load_rel
  out$solids_pct <- 100 * out$eps_s / (out$eps_s + out$eps_l)
  out$regime <- vapply(seq_len(nrow(out)),
                       function(i) regime_of(out$z[i], load_rel[i], p), character(1))
  out
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
H_int   <- params[["H_pool"]]
t_pool  <- approx(out_psd$z, out_psd$t_res, H_int)$y
t_tot   <- out_psd$t_res[nrow(out_psd)]
top <- out_psd[nrow(out_psd), ]; ini <- out_psd[1, ]

cat(sprintf("Decant time (heuristic)= %.2f h  [100*%.1f/(%.0f-%.0f)]\n",
            t_decant_h(params), params[["mu_cont"]], params[["rho_A"]], params[["rho_B"]]))
cat(sprintf("Pool depth / interface = %.2f m (fixed decant depth)\n", H_int))
cat(sprintf("Residence: pool %.2f h + foam %.2f h = %.2f h total\n",
            t_pool/3600, (t_tot - t_pool)/3600, t_tot/3600))
cat(sprintf("Solids content: inlet %.1f%% -> top %.1f%%\n", ini$solids_pct, top$solids_pct))
cat(sprintf("Plug-flow height frac  = %.0f%%\n", 100*mean(out_psd$regime=="plug")))
cat(sprintf("Retention: fine=%.0f%% mid=%.0f%% coarse=%.0f%% | impurity %.0f->%.1f%%\n",
            100*top$Js_fine/ini$Js_fine, 100*top$Js_mid/ini$Js_mid,
            100*top$Js_crs/ini$Js_crs, ini$C_imp, top$C_imp))

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

plot(out_psd$z, out_psd$Js_mid, type = "l", col = "green", lwd = 3, ylim = c(0, 0.05),
     ylab = "Solid flux (m/s)", xlab = "Height (m)", main = "PSD retention")
lines(out_psd$z, out_psd$Js_fine, col = "blue", lwd = 2, lty = 2)
lines(out_psd$z, out_psd$Js_crs,  col = "red",  lwd = 2, lty = 3)
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)
legend("right", legend = c("Mid","Fine","Coarse"), col = c("green","blue","red"),
       lty = c(1,2,3), lwd = 2)

plot(out_psd$z, out_psd$solids_pct, type = "l", col = "darkorange", lwd = 3,
     ylim = c(0, 50), ylab = "%", xlab = "Height (m)", main = "Solids content & impurity")
lines(out_psd$z, out_psd$C_imp / 2, col = "purple", lwd = 2, lty = 2)
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)
legend("topright", legend = c("Solids %","Impurity %/2"), col = c("darkorange","purple"),
       lty = c(1,2), lwd = 2)

plot(out_psd$z, out_psd$t_res/3600, type = "l", col = "black", lwd = 3,
     ylab = "Cumulative residence (h)", xlab = "Height (m)",
     main = "Residence (decant-limited)")
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)
plug_z <- out_psd$z[out_psd$regime == "plug"]
if (length(plug_z)) rug(plug_z, col = "forestgreen", lwd = 2)
legend("topleft", legend = c("t_res","interface","plug region"),
       col = c("black","gray","forestgreen"), lty = c(1,2,1), lwd = 2)
