library(deSolve)

# =====================================================================
# SCRIPT 1: PREFORMED-FOAM WASH COLUMN (PSD retention + impurity wash)
# =====================================================================
# The foam is PREFORMED and PARTICLE-LOADED in an upstream step; this column
# has NO gas sparger of its own. Loaded foam enters at the base, rises as a
# plug, drains (dries) as it goes, and wash liquid added at the top sweeps
# entrained impurity out. The column therefore only RETAINS-vs-LOSES what
# arrives -- there is no in-column attachment and no foam-formation interface.
#
# State marched up height z:
#   Js_fine, Js_mid, Js_crs   attached solid flux per size class [m/s]
#   C_imp                     entrained impurity [%]
#   eps_l                     liquid holdup of the foam [-]
#
# Boundary conditions at z = 0 are the PREFORMED loaded foam (Js_*_in, etc.).
#
# Physics choices (PLACEHOLDER closures, tune later):
#   * Preforming pressure sets the INLET bubble size: d_b ~ (P_ref/P)^(1/3).
#     Smaller bubbles -> tighter Plateau borders (less channeling) but wetter
#     foam (higher drained holdup). Bubble size is treated as constant up the
#     column here (coalescence up the column is a future refinement).
#   * Drainage: eps_l relaxes from the inlet foam wetness toward a drained
#     equilibrium eps_l_eq over a length L_drain (foam dries as it rises).
#   * COLLAPSE / FLOOD: if the foam wetness reaches the packing limit
#     (eps_g <= eps_g_pack) the froth collapses to bubbly liquid. This can be
#     driven by too-wet inlet foam or by a drained equilibrium eps_l_eq that
#     sits above the packing holdup (high preforming pressure ~ 7.8 bar).
#   * COARSE LOSS = two mechanisms: (1) buoyancy-failure / intrinsic
#     detachment k_det (density too high to stay attached) and (2)
#     bubble-collapse loss k_col, gated by foam dryness (coalescence grows as
#     eps_g rises above packing). Fines are buoyant -> both small.
#   * Wash liquid is set by the foam feed's entrained liquid
#     (Jl_wash = wash_ratio * J_foam * eps_l_in) so it scales with feed;
#     channeling from Plateau-border flow maldistribution (Q ~ r^4 => ~4x CV).
# =====================================================================

params <- c(
  H_total = 5.0,

  # --- PREFORMED FOAM FEED (no sparger; foam made upstream) ---
  J_foam   = 0.06,       # superficial foam velocity entering the base [m/s]
  eps_l_in = 0.20,       # liquid holdup of the incoming foam [-] (eps_g = 0.80)
  # inlet solid loading of the preformed foam [m/s] (upstream PSD selectivity)
  Js_fine_in = 0.017, Js_mid_in = 0.040, Js_crs_in = 0.020,
  C_imp_in   = 100,      # inlet entrained impurity [%]

  # --- PREFORMING PRESSURE -> INLET BUBBLE SIZE ---
  P = 150000, P_ref = 150000, d_b_ref = 2.0e-3,

  # --- DRAINAGE (foam dries as it rises) ---
  eps_l_dry_ref = 0.15,  # drained residual holdup at P_ref [-] (eps_g = 0.85)
  L_drain       = 1.0,   # drainage relaxation length [m]
  eps_g_pack    = 0.74,  # packing / collapse (flood) threshold gas fraction

  # --- WASH LIQUID (added at top, correlated to the foam feed) ---
  wash_ratio = 0.5,      # wash-water : entrained-feed-liquid ratio

  # --- PSD LOSS: buoyancy-failure detachment + bubble-collapse ---
  # (1) intrinsic detachment (can't float; density too high) -> coarse worst
  k_det_fine = 0.01, k_det_mid = 0.05, k_det_crs = 1.0,
  # (2) bubble-collapse loss (dumped when films rupture) -> coarse worst
  k_col_fine = 0.05, k_col_mid = 0.30, k_col_crs = 2.0,

  # --- CHANNELING / WASH DISTRIBUTION ---
  k_drainage = 1.5,      # base impurity sweep rate [1/m]
  channel_cv = 0.35      # CV of Plateau-border widths (0 = ideal plug wash)
)

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    # inlet bubble size from preforming pressure -> foam wetness & channeling
    d_b      <- d_b_ref * (P_ref / P)^(1 / 3)
    eps_l_eq <- eps_l_dry_ref * (d_b_ref / d_b)   # wetter drained foam at high P

    eps_g <- 1 - eps_l

    # foam dries as it rises (drainage relaxation toward equilibrium)
    deps_l <- -(eps_l - eps_l_eq) / L_drain

    # bubble-collapse rate grows as the foam dries past packing
    coalescence <- max(0, eps_g - eps_g_pack)

    # size loss: buoyancy-failure detachment + collapse dumping
    dJs_fine <- -Js_fine * (k_det_fine + k_col_fine * coalescence)
    dJs_mid  <- -Js_mid  * (k_det_mid  + k_col_mid  * coalescence)
    dJs_crs  <- -Js_crs  * (k_det_crs  + k_col_crs  * coalescence)

    # washing: swept by wash liquid, penalized by Plateau-border channeling
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    wash_strength   <- wash_ratio / (wash_ratio + 1)   # = Jl_wash/(Jl_wash+Jl_in)
    dCimp <- -k_drainage * Wash_Efficiency * wash_strength * C_imp

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l))
  })
}

# --- helpers ---------------------------------------------------------------
d_b_of_P    <- function(p) p[["d_b_ref"]] * (p[["P_ref"]] / p[["P"]])^(1/3)
epsleq_of_P <- function(p) p[["eps_l_dry_ref"]] * (p[["d_b_ref"]] / d_b_of_P(p))
# foam collapses if the inlet OR the drained-equilibrium holdup can't pack
is_flooded  <- function(p) {
  thr <- 1 - p[["eps_g_pack"]]
  (p[["eps_l_in"]] >= thr) || (epsleq_of_P(p) >= thr)
}

solve_column <- function(p) {
  state_init <- c(Js_fine = p[["Js_fine_in"]], Js_mid = p[["Js_mid_in"]],
                  Js_crs = p[["Js_crs_in"]], C_imp = p[["C_imp_in"]],
                  eps_l = p[["eps_l_in"]])
  z_seq <- seq(0, p[["H_total"]], by = 0.02)
  as.data.frame(ode(y = state_init, times = z_seq,
                    func = column_model_psd, parms = p))
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
names(out_psd)[1] <- "z"
out_psd$eps_g <- 1 - out_psd$eps_l

U_int <- params[["J_foam"]] / (1 - params[["eps_l_in"]])   # interstitial rise [m/s]
cat(sprintf("Inlet bubble size d_b   = %.0f um\n", d_b_of_P(params) * 1e6))
cat(sprintf("Foam rise (interstitial)= %.3f m/s  -> residence %.1f s over %.1f m\n",
            U_int, params[["H_total"]] / U_int, params[["H_total"]]))
cat(sprintf("Drained holdup eps_l_eq = %.3f  (collapse if >= %.2f)\n",
            epsleq_of_P(params), 1 - params[["eps_g_pack"]]))
cat(sprintf("Foam collapses/floods?  = %s\n", ifelse(is_flooded(params), "YES", "no")))
top <- out_psd[nrow(out_psd), ]; ini <- out_psd[1, ]
cat("Recovery (retained fraction of what entered):\n")
cat(sprintf("  fine=%.1f%%  mid=%.1f%%  coarse=%.1f%%  | impurity %.0f%% -> %.1f%%\n",
            100*top$Js_fine/ini$Js_fine, 100*top$Js_mid/ini$Js_mid,
            100*top$Js_crs/ini$Js_crs, ini$C_imp, top$C_imp))

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# 1. PSD retention up the column
plot(out_psd$z, out_psd$Js_mid, type = "l", col = "green", lwd = 3, ylim = c(0, 0.05),
     ylab = "Solid flux (m/s)", xlab = "Height (m)", main = "PSD retention (foam wash)")
lines(out_psd$z, out_psd$Js_fine, col = "blue", lwd = 2, lty = 2)
lines(out_psd$z, out_psd$Js_crs,  col = "red",  lwd = 2, lty = 3)
legend("right", legend = c("Mid", "Fine", "Coarse"),
       col = c("green", "blue", "red"), lty = c(1, 2, 3), lwd = 2)

# 2. Wash impurity removal
plot(out_psd$z, out_psd$C_imp, type = "l", col = "purple", lwd = 3, ylim = c(0, 100),
     ylab = "Impurity (%)", xlab = "Height (m)", main = "Wash impurity removal")

# 3. Foam holdup / stability
plot(out_psd$z, out_psd$eps_g, type = "l", col = "black", lwd = 3, ylim = c(0, 1),
     ylab = "Gas fraction eps_g", xlab = "Height (m)", main = "Foam holdup / stability")
abline(h = params[["eps_g_pack"]], col = "orange", lty = 2)
legend("bottomright", legend = c("eps_g", "collapse (0.74)"),
       col = c("black", "orange"), lty = c(1, 2), lwd = 2)
