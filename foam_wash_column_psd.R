library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN (PSD retention + impurity wash)
# =====================================================================
# The foam is PREFORMED and PARTICLE-LOADED upstream; this column has NO gas
# sparger. Loaded foam enters at the base, rises as a plug, drains (dries)
# slowly as it goes, and wash liquid added at the top sweeps entrained
# impurity out. The column RETAINS-vs-LOSES what arrives -- no in-column
# attachment, no foam-formation interface.
#
# Two effects added in this version (both observed in practice):
#  (A) DRAINAGE IS TIME-LIMITED. Films drain on a timescale tau_drain that is
#      long compared with the foam residence time, so drainage is INCOMPLETE:
#      the foam stays wet and top-of-column solids content is low (~3-7%),
#      even though the foam is a coherent plug over most of the height. This
#      is the key point: PLUG-FLOW ONSET (bubble jamming at eps_g ~ 0.74) is
#      NOT the same as DRYING (eps_g -> 0.95+). Preformed foam enters already
#      jammed (eps_g_in ~ 0.80), so it is a plug from the base; it just never
#      finishes draining within the column.
#  (B) SNOW-GLOBE vs PLUG is a MATERIAL property. Well-stabilized, lightly/
#      buoyantly-loaded foam jams into a continuous rising plug. Weak films or
#      heavy/dense particle loads thin and rupture films, so the foam cannot
#      hold a continuous network and stays as discrete particle-armored wet
#      flakes ("snow-globe"). Set by film_stability vs. solid loading.
#
# State marched up height z: Js_fine, Js_mid, Js_crs [m/s]; C_imp [%];
#                            eps_l (liquid holdup) [-].
# Solid holdup eps_s and solids content are derived diagnostics.
#
# PLACEHOLDER closures throughout -- tune against data.
# =====================================================================

params <- c(
  H_total = 5.0,

  # --- PREFORMED FOAM FEED (no sparger; foam made upstream) ---
  J_foam   = 0.06,       # superficial foam velocity entering the base [m/s]
  eps_l_in = 0.20,       # liquid holdup of the incoming foam [-] (eps_g = 0.80)
  eps_s_in = 0.020,      # solid volume fraction of incoming foam [-] (for solids%)
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020,  # inlet loading [m/s]
  C_imp_in   = 100,      # inlet entrained impurity [%]

  # --- PREFORMING PRESSURE -> INLET BUBBLE SIZE ---
  P = 150000, P_ref = 150000, d_b_ref = 2.0e-3,

  # --- DRAINAGE (time-limited: foam does NOT fully drain in the column) ---
  eps_l_dry_ref = 0.15,  # fully-drained residual holdup at P_ref [-]
  tau_drain     = 90,    # film drainage time constant [s] (>> residence -> incomplete)
  eps_g_pack    = 0.74,  # packing / collapse threshold gas fraction

  # --- WASH LIQUID (added at top, correlated to the foam feed) ---
  wash_ratio = 0.5,      # wash-water : entrained-feed-liquid ratio

  # --- PSD LOSS: buoyancy-failure detachment + bubble-collapse ---
  k_det_fine = 0.01, k_det_mid = 0.05, k_det_crs = 1.0,
  k_col_fine = 0.05, k_col_mid = 0.30, k_col_crs = 2.0,

  # --- CHANNELING / WASH DISTRIBUTION ---
  k_drainage = 1.5,      # base impurity sweep rate [1/m]
  channel_cv = 0.35,     # CV of Plateau-border widths (0 = ideal plug wash)

  # --- REGIME (snow-globe vs plug) : material stability vs loading ---
  film_stability = 0.76, # foam/film stability (material; higher -> plug)
  load_sens      = 1.0,  # sensitivity of the plug criterion to solid loading
  plug_crit      = 1.0   # plug_index threshold above which foam is a plug
)

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    d_b      <- d_b_ref * (P_ref / P)^(1 / 3)
    eps_l_eq <- eps_l_dry_ref * (d_b_ref / d_b)   # wetter drained foam at high P

    eps_g  <- 1 - eps_l
    U_int  <- J_foam / max(eps_g, 1e-3)           # interstitial foam velocity [m/s]

    # (A) time-limited drainage: dz = U_int*dt, so length scale = U_int*tau_drain
    deps_l <- -(eps_l - eps_l_eq) / (U_int * tau_drain)

    # bubble-collapse rate grows as the foam dries past packing
    coalescence <- max(0, eps_g - eps_g_pack)
    dJs_fine <- -Js_fine * (k_det_fine + k_col_fine * coalescence)
    dJs_mid  <- -Js_mid  * (k_det_mid  + k_col_mid  * coalescence)
    dJs_crs  <- -Js_crs  * (k_det_crs  + k_col_crs  * coalescence)

    # washing: swept by wash liquid, penalized by Plateau-border channeling
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    wash_strength   <- wash_ratio / (wash_ratio + 1)
    dCimp <- -k_drainage * Wash_Efficiency * wash_strength * C_imp

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l))
  })
}

# --- helpers ---------------------------------------------------------------
d_b_of_P    <- function(p) p[["d_b_ref"]] * (p[["P_ref"]] / p[["P"]])^(1/3)
epsleq_of_P <- function(p) p[["eps_l_dry_ref"]] * (p[["d_b_ref"]] / d_b_of_P(p))

# (B) regime map. Jamming (eps_g >= packing) is the GATE for foam at all;
# above it, coherence is a competition between film stability and solid
# loading: plug_index = film_stability / (load_sens * relative loading).
# Stable, lightly-loaded films lock into a continuous PLUG; weak films or
# heavy loading stay as discrete wet SNOW-GLOBE flakes. Loading is highest at
# the base (before coarse sheds), so a marginal material is snow-globe low
# down and consolidates into a plug as it rises, drains and lightens.
regime_of <- function(eps_g, load_rel, p) {
  if (eps_g < p[["eps_g_pack"]]) return("dispersed")
  plug_index <- p[["film_stability"]] / (p[["load_sens"]] * load_rel)
  if (plug_index >= p[["plug_crit"]]) "plug" else "snowglobe"
}

solve_column <- function(p) {
  state_init <- c(Js_fine = p[["Js_fine_in"]], Js_mid = p[["Js_mid_in"]],
                  Js_crs = p[["Js_crs_in"]], C_imp = p[["C_imp_in"]],
                  eps_l = p[["eps_l_in"]])
  z_seq <- seq(0, p[["H_total"]], by = 0.02)
  out <- as.data.frame(ode(y = state_init, times = z_seq,
                           func = column_model_psd, parms = p))
  names(out)[1] <- "z"
  out$eps_g <- 1 - out$eps_l
  # derived: solid holdup scales with retained loading; solids content of the
  # condensed (liquid+solid) phase = what a collapsed sample would assay.
  Js_tot_in <- p[["Js_fine_in"]] + p[["Js_mid_in"]] + p[["Js_crs_in"]]
  load_rel  <- (out$Js_fine + out$Js_mid + out$Js_crs) / Js_tot_in
  out$eps_s <- p[["eps_s_in"]] * load_rel
  out$solids_pct <- 100 * out$eps_s / (out$eps_s + out$eps_l)
  out$regime <- vapply(seq_len(nrow(out)),
                       function(i) regime_of(out$eps_g[i], load_rel[i], p), character(1))
  out
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
U_int <- params[["J_foam"]] / (1 - params[["eps_l_in"]])
t_res <- params[["H_total"]] / U_int
drain_complete <- 1 - exp(-t_res / params[["tau_drain"]])
plug_frac <- mean(out_psd$regime == "plug")
top <- out_psd[nrow(out_psd), ]; ini <- out_psd[1, ]

cat(sprintf("Residence time         = %.0f s ; drainage time tau = %.0f s\n",
            t_res, params[["tau_drain"]]))
cat(sprintf("Drainage completeness  = %.0f%%  (1=fully drained)\n", 100*drain_complete))
cat(sprintf("Solids content: inlet %.1f%% -> top %.1f%%\n", ini$solids_pct, top$solids_pct))
cat(sprintf("Plug-flow height frac  = %.0f%%  (rest: %s)\n", 100*plug_frac,
            paste(unique(out_psd$regime[out_psd$regime != "plug"]), collapse=", ")))
cat(sprintf("Retention: fine=%.0f%% mid=%.0f%% coarse=%.0f%% | impurity %.0f->%.1f%%\n",
            100*top$Js_fine/ini$Js_fine, 100*top$Js_mid/ini$Js_mid,
            100*top$Js_crs/ini$Js_crs, ini$C_imp, top$C_imp))

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# 1. PSD retention
plot(out_psd$z, out_psd$Js_mid, type = "l", col = "green", lwd = 3, ylim = c(0, 0.05),
     ylab = "Solid flux (m/s)", xlab = "Height (m)", main = "PSD retention")
lines(out_psd$z, out_psd$Js_fine, col = "blue", lwd = 2, lty = 2)
lines(out_psd$z, out_psd$Js_crs,  col = "red",  lwd = 2, lty = 3)
legend("right", legend = c("Mid", "Fine", "Coarse"),
       col = c("green", "blue", "red"), lty = c(1, 2, 3), lwd = 2)

# 2. Solids content (drainage) + impurity
plot(out_psd$z, out_psd$solids_pct, type = "l", col = "darkorange", lwd = 3,
     ylim = c(0, 50), ylab = "%", xlab = "Height (m)",
     main = "Solids content & impurity")
lines(out_psd$z, out_psd$C_imp / 2, col = "purple", lwd = 2, lty = 2)  # scaled /2
legend("topright", legend = c("Solids content %", "Impurity %/2"),
       col = c("darkorange", "purple"), lty = c(1, 2), lwd = 2)

# 3. Foam holdup + regime (plug where eps_g jammed & coherent)
plot(out_psd$z, out_psd$eps_g, type = "l", col = "black", lwd = 3, ylim = c(0, 1),
     ylab = "Gas fraction eps_g", xlab = "Height (m)", main = "Holdup & regime")
abline(h = params[["eps_g_pack"]], col = "orange", lty = 2)
plug_z <- out_psd$z[out_psd$regime == "plug"]
if (length(plug_z)) rug(plug_z, col = "forestgreen", lwd = 2)
legend("bottomright", legend = c("eps_g", "packing (0.74)", "plug region"),
       col = c("black", "orange", "forestgreen"), lty = c(1, 2, 1), lwd = 2)
