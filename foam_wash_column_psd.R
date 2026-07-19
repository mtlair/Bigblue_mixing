library(deSolve)

# =====================================================================
# SCRIPT 1: STEADY-STATE PSD & CHANNELING MODEL (foam wash column)
# =====================================================================
# Marches a lumped steady-state model up the column height z:
#   Zone 1 (z < H_int): liquid pool, particles attach toward a per-size
#                       carrying capacity.
#   Zone 2 (z >= H_int): foam washing zone, particles detach / entrain and
#                        the impurity is swept out by the added wash liquid.
#
# Notes on the physics choices (from design discussion):
#   * Pressure sets bubble size: d_b ~ (P_ref/P)^(1/3) (gas compression at
#     formation). Smaller bubbles -> (a) more interfacial area, so higher
#     attachment rate and carrying capacity, and (b) smaller, more uniform
#     Plateau borders, so less wash channeling.
#   * channel_cv (was "D_axial") is a placeholder for POORLY DISTRIBUTED
#     WASH FLOW. Plateau-border size follows the packing fraction of the
#     squished bubbles (liquid holdup) and scales with bubble size; drainage
#     sets the *relative* channel width. Because drainage through a Plateau
#     border is ~Poiseuille (Q ~ r^4), a spread channel_cv in border width
#     amplifies ~4x into wash-flow variability -> the widest (largest)
#     borders take most of the water and the rest of the foam is
#     under-washed. This is a flow-maldistribution efficiency, NOT an
#     axial-dispersion (2nd-order) term.
#   * The wash liquid is set by the feed rate (Jl_wash = wash_ratio * Jl_feed)
#     so it can be swept in proportion to what is fed.
#   * Size-dependent loss: heavy COARSE particles settle / drain out with the
#     down-flowing wash liquid and are the dominant product loss; FINE
#     particles are buoyant enough to stay dispersed on the foam, so their
#     loss is small. Loss is coupled to the wash flux (over-washing costs
#     coarse recovery).
# =====================================================================

params <- c(
  H_total = 5.0, H_int = 2.0,
  Jg = 0.05,             # superficial gas velocity [m/s]
  eps_g = 0.85,          # gas holdup in foam (packing fraction of squished bubbles)

  # --- PRESSURE -> BUBBLE SIZE ---
  P       = 150000,      # operating pressure [Pa]
  P_ref   = 150000,      # reference pressure for the bubble-size closure [Pa]
  d_b_ref = 2.0e-3,      # bubble diameter at P_ref [m]

  # --- WASH LIQUID (correlated to feed so impurity can be swept) ---
  Jl_feed    = 0.004,    # feed superficial liquid velocity [m/s]
  wash_ratio = 0.5,      # wash-water : feed ratio -> sets the wash flux

  # --- 1. PARTICLE SIZE DISTRIBUTION (PSD) KINETICS ---
  # sink_* = susceptibility to settling/draining OUT of the foam with the
  # down-flowing wash liquid. Heavy coarse are lost most; buoyant fines least.
  # FINE (< 20 um): Hard to attach (low collision), hard to detach (embedded);
  #                 buoyant, so it rides the foam and is barely lost.
  k_att_fine = 1.0,  k_det_fine = 0.01, sink_fine = 0.02,
  # MID (50-150 um): The sweet spot. High attachment, moderate detachment.
  k_att_mid  = 8.0,  k_det_mid  = 0.1,  sink_mid  = 0.15,
  # COARSE (> 300 um): High collision, but massive detachment (gravity/shear)
  #                    and settles/drains out -> the dominant product loss.
  k_att_crs  = 5.0,  k_det_crs  = 5.0,  sink_crs  = 0.80,

  # Max carrying capacity per size class [m/s] (at P_ref; scales with area)
  Js_max_fine = 0.02, Js_max_mid = 0.04, Js_max_crs = 0.02,

  # --- 2. CHANNELING / WASH DISTRIBUTION ---
  k_drainage = 1.5,      # base impurity sweep rate [1/m]
  channel_cv = 0.35      # CV of Plateau-border widths (0 = ideal plug wash)
)

column_model_psd <- function(z, state, parameters) {
  with(as.list(c(state, parameters)), {

    # Pressure sets bubble size (gas compression at formation): d_b shrinks
    # as P rises. Smaller bubbles -> more interfacial area (area_factor > 1).
    d_b         <- d_b_ref * (P_ref / P)^(1 / 3)
    area_factor <- d_b_ref / d_b            # interfacial-area ratio vs reference

    # Wash liquid flux is set by the feed rate (swept downward through foam).
    Jl_wash <- wash_ratio * Jl_feed        # magnitude of downward wash flux [m/s]

    if (z < H_int) {
      # --- ZONE 1: LIQUID POOL ---
      # Fines struggle to attach; coarse and mid attach rapidly toward capacity.
      # Smaller bubbles (higher P) add interfacial area: faster attachment and
      # higher carrying capacity.
      dJs_fine <- k_att_fine * area_factor * (Js_max_fine * area_factor - Js_fine)
      dJs_mid  <- k_att_mid  * area_factor * (Js_max_mid  * area_factor - Js_mid)
      dJs_crs  <- k_att_crs  * area_factor * (Js_max_crs  * area_factor - Js_crs)

      dCimp_dz <- 0                         # well-mixed pool, no washing yet

    } else {
      # --- ZONE 2: FOAM WASHING ZONE ---
      U_rise <- Jg / eps_g                  # interstitial foam rise velocity [m/s]

      # Drainage ratio: downward wash liquid vs upward foam rise (dimensionless).
      drainage_ratio <- Jl_wash / U_rise

      # Detachment + settling/drainage loss. Heavy coarse settle and drain out
      # with the wash liquid (dominant loss); buoyant fines mostly stay.
      dJs_fine <- -Js_fine * (k_det_fine + sink_fine * drainage_ratio)
      dJs_mid  <- -Js_mid  * (k_det_mid  + sink_mid  * drainage_ratio)
      dJs_crs  <- -Js_crs  * (k_det_crs  + sink_crs  * drainage_ratio)

      # Channeling penalty from Plateau-border flow maldistribution. Border
      # width scales with bubble size, so larger bubbles (lower P) widen the
      # borders and worsen maldistribution. Q ~ r^4 => ~4x CV amplification.
      flow_cv         <- 4 * channel_cv * (d_b / d_b_ref)
      Wash_Efficiency <- 1 / (1 + flow_cv^2)   # effectiveness under maldistribution

      # Impurity swept out by the wash liquid; strength set by wash-vs-feed.
      wash_strength <- Jl_wash / (Jl_wash + Jl_feed)
      dCimp_dz <- -k_drainage * Wash_Efficiency * wash_strength * C_imp
    }

    list(c(dJs_fine, dJs_mid, dJs_crs, dCimp_dz))
  })
}

# --- Solve in two segments so the zone interface is not a discontinuity
#     inside a single adaptive solve. ---
H_int_v   <- params[["H_int"]]
H_total_v <- params[["H_total"]]
state_init <- c(Js_fine = 0.001, Js_mid = 0.001, Js_crs = 0.001, C_imp = 100)

# Zone 1: pool  (0 -> H_int)
z1   <- seq(0, H_int_v, by = 0.05)
out1 <- ode(y = state_init, times = z1, func = column_model_psd, parms = params)

# Zone 2: foam  (H_int -> H_total), starting from the pool-outlet state
state_mid <- out1[nrow(out1), -1]
z2   <- seq(H_int_v, H_total_v, by = 0.05)
out2 <- ode(y = state_mid, times = z2, func = column_model_psd, parms = params)

# Stitch the two segments (drop the duplicated H_int row)
out_psd <- as.data.frame(rbind(out1, out2[-1, ]))
names(out_psd)[1] <- "z"

# =====================================================================
# PLOTS
# =====================================================================
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Plot PSD recovery
plot(out_psd$z, out_psd$Js_mid, type = "l", col = "green", lwd = 3, ylim = c(0, 0.05),
     ylab = "Solid flux (m/s)", xlab = "Height (m)", main = "Particle size recovery")
lines(out_psd$z, out_psd$Js_fine, col = "blue", lwd = 2, lty = 2)
lines(out_psd$z, out_psd$Js_crs,  col = "red",  lwd = 2, lty = 3)
abline(v = H_int_v, col = "gray", lty = 2)
legend("topright", legend = c("Mid (optimal)", "Fine", "Coarse"),
       col = c("green", "blue", "red"), lty = c(1, 2, 3), lwd = 2)

# Plot channeling / wash impact
plot(out_psd$z, out_psd$C_imp, type = "l", col = "purple", lwd = 3, ylim = c(0, 100),
     ylab = "Impurity (%)", xlab = "Height (m)", main = "Wash impurity removal")
abline(v = H_int_v, col = "gray", lty = 2)
legend("topright", legend = "Impurity", col = "purple", lty = 1, lwd = 2)
