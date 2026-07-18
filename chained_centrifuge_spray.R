#!/usr/bin/env Rscript
# =============================================================================
# Chained model: decanter centrifuge -> reslurry (dilution) -> spray dryer.
#
# Sources BOTH model definitions (not their Morris screening runs) into isolated
# environments, maps the centrifuge stream handoff into the spray-dryer feed,
# and runs them end-to-end. Used to decompose WHY the dryer duty must drop to
# hit the same final particle size - i.e. which centrifuge-side driver (less
# entrained gas, higher slurry viscosity, retained plasticizer -> stickiness)
# moves the dryer outlet temperature / DeltaT and the drying difficulty.
#
# Run:  Rscript chained_centrifuge_spray.R
# =============================================================================

# ---- source model DEFINITIONS only, into isolated environments --------------
src_defs <- function(path, stop_marker) {
  L <- readLines(path); cut <- grep(stop_marker, L)[1]
  paste(L[seq_len(cut - 1)], collapse = "\n")
}
cen_env <- new.env(); spray_env <- new.env()
eval(parse(text = src_defs("centrifuge_morris_sensitivity.R", "^# 3\\. SMOKE TEST")), envir = cen_env)
eval(parse(text = src_defs("morris_sensitivity_analysis.R",   "^# 3\\. Morris design")), envir = spray_env)

unified_centrifuge_model <- cen_env$unified_centrifuge_model
centrifuge_to_spray      <- cen_env$centrifuge_to_spray
cen_nominal              <- cen_env$nominal_vec
spray_dry_model          <- spray_env$spray_dry_model
spray_factors            <- spray_env$factors

# ---- which spray inputs come from the centrifuge STREAM vs the dryer ---------
from_centrifuge <- c("rho_L", "C_solid_mass", "alpha_g_0", "sigma", "D_b", "mu_L",
                     "C_surfactant", "C_monomer", "C_plasticizer", "C_binder",
                     "I_strength", "Delta_pH")
# dryer operating / formulation settings: midpoints of the spray model's ranges
sp_mid <- setNames((spray_factors$min + spray_factors$max) / 2, spray_factors$name)

clamp_to_spray <- function(x) {
  for (nm in names(x)) {
    lo <- spray_factors$min[spray_factors$name == nm]
    hi <- spray_factors$max[spray_factors$name == nm]
    if (length(lo)) x[nm] <- min(max(x[nm], lo), hi)
  }
  x
}

# ---- the chained model ------------------------------------------------------
chained_model <- function(cen_row, reslurry_solids = 0.30, spray_op = sp_mid,
                          clamp = TRUE) {
  hand <- centrifuge_to_spray(cen_row, target_solid_mass = reslurry_solids)
  x <- spray_op
  for (nm in from_centrifuge) x[nm] <- hand[[nm]]
  if (clamp) x <- clamp_to_spray(x)
  list(handoff = hand, spray = spray_dry_model(x), spray_in = x)
}

# =============================================================================
# 1. Nominal end-to-end run
# =============================================================================
res <- chained_model(cen_nominal, reslurry_solids = 0.30)
sp  <- res$spray
cat("== Chained centrifuge -> reslurry(30% solids) -> spray dryer (nominal) ==\n")
cat(sprintf("  Feed to dryer:  rho_L %.0f  C_solid %.2f  alpha_g0 %.3f  mu_L %.4f  sigma %.3f  D_b %.0f um\n",
            res$handoff[["rho_L"]], res$handoff[["C_solid_mass"]], res$handoff[["alpha_g_0"]],
            res$handoff[["mu_L"]], res$handoff[["sigma"]], res$handoff[["D_b"]]*1e6))
cat(sprintf("  Dry powder:     D_particle %.1f um  X_moisture %.3f  theta_skin %.3f  Tg_prod %.0f K\n",
            sp[["D_particle_um"]], sp[["X_moisture"]], sp[["theta_skin_z"]], sp[["Tg_product_K"]]))
cat(sprintf("  Dryer state:    T_out %.1f K  DeltaT %.1f K  RH_out %.2f  rho_tapped %.0f\n\n",
            sp[["T_out_K"]], sp[["DeltaT_K"]], sp[["RH_out"]], sp[["rho_tapped"]]))

# =============================================================================
# 2. Driver decomposition: which centrifuge-side lever makes it "hard to dry"?
# =============================================================================
# Perturb one physical driver at a time (via the centrifuge inputs / reslurry)
# and read the spray-side response. "Hard to dry" shows up as larger droplets,
# higher residual moisture, more skin (case hardening), and it changes DeltaT.
metrics <- function(r) c(D_particle = r$spray[["D_particle_um"]],
                         X_moist    = r$spray[["X_moisture"]],
                         theta_skin = r$spray[["theta_skin_z"]],
                         Tg_prod    = r$spray[["Tg_product_K"]],
                         stick_marg = (0.85*r$spray[["T_out_K"]] + 0.15*r$spray_in[["T_feed"]]) -
                                      r$spray[["Tg_product_K"]],   # T_particle - Tg (>0 = sticky)
                         DeltaT     = r$spray[["DeltaT_K"]])

base <- chained_model(cen_nominal, 0.30)

# (a) LESS entrained gas: stronger mechanical degassing (pop) at the centrifuge
low_gas <- cen_nominal; low_gas[["pop_frac"]] <- 0.90; low_gas[["feed_gas_frac"]] <- 0.10
# (b) HIGHER slurry viscosity: reslurry to a higher solids (less dilution)
#     -> higher mu_L and rho_L into the atomizer
hi_visc <- cen_nominal
# (c) MORE retained plasticizer -> stickier powder (case hardening)
hi_plast <- cen_nominal; hi_plast[["plasticizer_frac"]] <- 0.25

rows <- rbind(
  base        = metrics(base),
  less_gas    = metrics(chained_model(low_gas,  0.30)),
  higher_visc = metrics(chained_model(hi_visc,  0.42)),   # 42% solids reslurry
  more_plast  = metrics(chained_model(hi_plast, 0.30))
)
cat("== Driver decomposition (one lever at a time) ==\n")
print(round(rows, 3))
cat("\n  D_droplet/D_particle in um; X_moist & theta_skin in [-]; DeltaT in K.\n")

# =============================================================================
# 3. DeltaT vs dryer feed rate (the "lower the feed" observation)
# =============================================================================
cat("\n== DeltaT vs dryer feed rate mdot_L (at the nominal feed) ==\n")
cat(sprintf("  %-10s %-8s %-10s %-8s\n", "mdot_L", "DeltaT", "D_particle", "X_moist"))
for (m in c(0.004, 0.008, 0.012, 0.016, 0.020)) {
  op <- sp_mid; op[["mdot_L"]] <- m
  r <- chained_model(cen_nominal, 0.30, spray_op = op)
  cat(sprintf("  %-10.3f %-8.1f %-10.1f %-8.3f\n",
              m, r$spray[["DeltaT_K"]], r$spray[["D_particle_um"]], r$spray[["X_moisture"]]))
}
cat("\nLower feed -> lower evaporation load -> smaller DeltaT (outlet runs hotter).\n")
