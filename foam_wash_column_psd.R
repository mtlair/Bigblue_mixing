library(deSolve)

# =====================================================================
# SCRIPT 1: STEADY-STATE PSD & CHANNELING MODEL (foam wash column)
# =====================================================================
# Marches a lumped steady-state model up the column height z. The liquid
# holdup eps_l(z) is now a STATE VARIABLE (drainage profile), so the
# foam-liquid interface, the foam regimes, and the pressure flood limit are
# all EMERGENT instead of hard-coded.
#
# Regime is set by the local gas fraction eps_g = 1 - eps_l:
#   eps_g <  0.74           bubbly liquid POOL      -> particle ATTACHMENT
#   0.74 <= eps_g < ~0.90   wet ("snow-globe") foam -> detach + wash
#   eps_g >= ~0.90          drier / coalescing foam -> + bubble-collapse loss
# The pool<->foam switch uses a smooth blend (foam_frac) so the ODE RHS is
# continuous and a single solve is valid.
#
# Physics choices (from design discussion) -- PLACEHOLDER closures, tune later:
#   * Pressure sets bubble size: d_b ~ (P_ref/P)^(1/3) (gas compression at
#     formation). Smaller bubbles -> more interfacial area (attachment,
#     capacity) AND wetter foam (higher eps_l,eq) AND tighter Plateau borders
#     (less channeling).
#   * Drainage: eps_l relaxes from the sparger value toward a drained
#     equilibrium eps_l_eq over a length L_drain.
#   * FLOOD LIMIT: if the drained equilibrium holdup eps_l_eq >= (1-eps_g_pack)
#     the foam can never pack into a froth -> the column floods to bubbly
#     liquid (no foam zone). With the reference closure this happens ~7.8 bar.
#   * COARSE LOSS is two mechanisms: (1) buoyancy failure / intrinsic
#     detachment k_det (density too high to stay attached), and (2)
#     bubble-collapse loss k_collapse, gated by foam dryness (coalescence
#     grows as eps_g rises above packing). Fines are buoyant -> both small.
#   * Wash liquid is set by feed (Jl_wash = wash_ratio * Jl_feed); channeling
#     from Plateau-border flow maldistribution (Q ~ r^4 => ~4x CV).
# =====================================================================

params <- c(
  H_total = 5.0,
  Jg = 0.05,             # superficial gas velocity [m/s]

  # --- PRESSURE -> BUBBLE SIZE ---
  P       = 150000,      # operating pressure [Pa]
  P_ref   = 150000,      # reference pressure for the bubble-size closure [Pa]
  d_b_ref = 2.0e-3,      # bubble diameter at P_ref [m]

  # --- HOLDUP / DRAINAGE (drives the emergent interface + regimes) ---
  eps_l_sparger = 0.60,  # liquid holdup at the sparger (bubbly pool) [-]
  eps_l_dry_ref = 0.15,  # drained foam residual holdup at P_ref [-] (eps_g=0.85)
  L_drain       = 1.0,   # drainage relaxation length [m]
  eps_g_pack    = 0.74,  # gas fraction where bubbles pack into foam / flood limit
  blend_w       = 0.02,  # width of the smooth pool<->foam switch [-]

  # --- WASH LIQUID (correlated to feed so impurity can be swept) ---
  Jl_feed    = 0.004,    # feed superficial liquid velocity [m/s]
  wash_ratio = 0.5,      # wash-water : feed ratio -> sets the wash flux

  # --- PSD ATTACHMENT (pool) & CARRYING CAPACITY ---
  k_att_fine = 1.0, k_att_mid = 8.0, k_att_crs = 5.0,
  Js_max_fine = 0.02, Js_max_mid = 0.04, Js_max_crs = 0.02,   # at P_ref

  # --- PSD LOSS (foam): buoyancy-failure detachment + bubble-collapse ---
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

    # --- pressure -> bubble size -> area & foam wetness ---
    d_b         <- d_b_ref * (P_ref / P)^(1 / 3)
    area_factor <- d_b_ref / d_b               # interfacial-area ratio vs ref
    eps_l_eq    <- eps_l_dry_ref * (d_b_ref / d_b)  # wetter foam at high P

    # --- local holdup / regime ---
    eps_g     <- 1 - eps_l
    # smooth pool(0)->foam(1) switch on the packing threshold
    foam_frac <- 1 / (1 + exp(-(eps_g - eps_g_pack) / blend_w))
    U_rise    <- Jg / max(eps_g, 1e-3)         # interstitial foam rise velocity

    # --- holdup drainage (relaxes toward the drained equilibrium) ---
    deps_l <- -(eps_l - eps_l_eq) / L_drain

    # --- POOL attachment (smaller bubbles add area -> faster, more capacity)
    att_fine <- k_att_fine * area_factor * (Js_max_fine * area_factor - Js_fine)
    att_mid  <- k_att_mid  * area_factor * (Js_max_mid  * area_factor - Js_mid)
    att_crs  <- k_att_crs  * area_factor * (Js_max_crs  * area_factor - Js_crs)

    # --- FOAM loss: buoyancy-failure detachment + collapse (dryness-gated) ---
    coalescence <- max(0, eps_g - eps_g_pack)  # drier foam coalesces more
    loss_fine <- -Js_fine * (k_det_fine + k_col_fine * coalescence)
    loss_mid  <- -Js_mid  * (k_det_mid  + k_col_mid  * coalescence)
    loss_crs  <- -Js_crs  * (k_det_crs  + k_col_crs  * coalescence)

    # blend pool and foam physics
    dJs_fine <- (1 - foam_frac) * att_fine + foam_frac * loss_fine
    dJs_mid  <- (1 - foam_frac) * att_mid  + foam_frac * loss_mid
    dJs_crs  <- (1 - foam_frac) * att_crs  + foam_frac * loss_crs

    # --- washing (foam only): swept by wash liquid, penalized by channeling
    flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)  # bigger bubbles worse
    Wash_Efficiency <- 1 / (1 + flow_cv^2)
    wash_strength   <- (wash_ratio * Jl_feed) / (wash_ratio * Jl_feed + Jl_feed)
    dCimp <- -foam_frac * k_drainage * Wash_Efficiency * wash_strength * C_imp

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp, deps_l))
  })
}

# --- helpers ---------------------------------------------------------------
d_b_of_P   <- function(p) p[["d_b_ref"]] * (p[["P_ref"]] / p[["P"]])^(1/3)
epsleq_of_P <- function(p) p[["eps_l_dry_ref"]] * (p[["d_b_ref"]] / d_b_of_P(p))
is_flooded  <- function(p) epsleq_of_P(p) >= (1 - p[["eps_g_pack"]])

# emergent foam-liquid interface: first z where eps_g crosses eps_g_pack
interface_height <- function(out, p) {
  eps_g <- 1 - out[, "eps_l"]
  thr   <- p[["eps_g_pack"]]
  above <- which(eps_g >= thr)
  if (length(above) == 0) return(NA_real_)   # flooded: never becomes foam
  i <- above[1]
  if (i == 1) return(out[1, "time"])
  # linear interpolation between the straddling nodes
  z0 <- out[i-1, "time"]; z1 <- out[i, "time"]
  g0 <- eps_g[i-1];       g1 <- eps_g[i]
  z0 + (thr - g0) * (z1 - z0) / (g1 - g0)
}

solve_column <- function(p) {
  state_init <- c(Js_fine = 0.001, Js_mid = 0.001, Js_crs = 0.001,
                  C_imp = 100, eps_l = p[["eps_l_sparger"]])
  z_seq <- seq(0, p[["H_total"]], by = 0.02)
  out <- ode(y = state_init, times = z_seq, func = column_model_psd, parms = p)
  as.data.frame(out)
}

# --- solve + report --------------------------------------------------------
out_psd <- solve_column(params)
names(out_psd)[1] <- "z"
out_psd$eps_g <- 1 - out_psd$eps_l

H_int <- interface_height(cbind(time = out_psd$z, eps_l = out_psd$eps_l), params)
cat(sprintf("Bubble size d_b        = %.0f um\n", d_b_of_P(params) * 1e6))
cat(sprintf("Drained holdup eps_l_eq= %.3f  (flood if >= %.2f)\n",
            epsleq_of_P(params), 1 - params[["eps_g_pack"]]))
cat(sprintf("Foam flooded?          = %s\n", ifelse(is_flooded(params), "YES", "no")))
cat(sprintf("Computed interface H_int = %s m\n",
            ifelse(is.na(H_int), "-- (flooded)", sprintf("%.2f", H_int))))
top <- out_psd[nrow(out_psd), ]
cat(sprintf("Top recovery: fine=%.4f mid=%.4f coarse=%.4f | impurity=%.1f%%\n",
            top$Js_fine, top$Js_mid, top$Js_crs, top$C_imp))

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# 1. PSD recovery
plot(out_psd$z, out_psd$Js_mid, type = "l", col = "green", lwd = 3, ylim = c(0, 0.05),
     ylab = "Solid flux (m/s)", xlab = "Height (m)", main = "Particle size recovery")
lines(out_psd$z, out_psd$Js_fine, col = "blue", lwd = 2, lty = 2)
lines(out_psd$z, out_psd$Js_crs,  col = "red",  lwd = 2, lty = 3)
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)
legend("topright", legend = c("Mid (optimal)", "Fine", "Coarse"),
       col = c("green", "blue", "red"), lty = c(1, 2, 3), lwd = 2)

# 2. Wash impurity removal
plot(out_psd$z, out_psd$C_imp, type = "l", col = "purple", lwd = 3, ylim = c(0, 100),
     ylab = "Impurity (%)", xlab = "Height (m)", main = "Wash impurity removal")
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)

# 3. Holdup / regime profile (the emergent interface)
plot(out_psd$z, out_psd$eps_g, type = "l", col = "black", lwd = 3, ylim = c(0, 1),
     ylab = "Gas fraction eps_g", xlab = "Height (m)", main = "Holdup / interface")
abline(h = params[["eps_g_pack"]], col = "orange", lty = 2)
if (!is.na(H_int)) abline(v = H_int, col = "gray", lty = 2)
legend("bottomright", legend = c("eps_g", "packing (0.74)"),
       col = c("black", "orange"), lty = c(1, 2), lwd = 2)
