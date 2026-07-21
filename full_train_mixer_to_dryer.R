#!/usr/bin/env Rscript
# =============================================================================
# FULL PROCESS TRAIN: gassed/templated mixer -> pressurized foam-wash column
#                     -> liquid-solid separator -> reslurry -> nozzle dryer
# =============================================================================
#
#   UP1 mixer            unified/up1_mixer_module.R      (pulled, deSolve ODE)
#     -> UP2 (foam-wash)   foam_wash_column()
#       -> UP3 (separator) centrifuge_morris_sensitivity.R (unified_centrifuge_model)
#         -> reslurry      centrifuge_to_spray()  (dilution to sprayable solids)
#           -> UP4 (dryer) morris_sensitivity_analysis.R (spray_dry_model)
#
# The mixer (UP1) is "the 2nd step prior to the separator"; the foam-wash column (UP2) is
# "the 1st step prior to the separator". Composition, gas holdup, surface
# chemistry, template state, temperature and (critically) the aggregate/floc
# strength are DEFINED ONCE at the mixer and flow downstream through the shared
# `stream` interface (unified/interface_stream.R), so UP3 (separator) and UP4 (dryer)
# see the mixer's transformed values instead of independent knobs.
#
# >>> foam_wash_column() is a deliberately thin PLACEHOLDER <<<
# The real pressurized foam-wash column model is being built separately. This
# stub implements only the column's DEFINING action - it removes entrained gas
# ("washes the foam out") and, being pressurized, collapses/shrinks the residual
# bubbles and sets the discharge pressure. Everything else is pass-through. Drop
# the real model in by replacing this one function; it already conforms to the
# stream contract f(stream, pars) -> stream.
#
# Run:  Rscript full_train_mixer_to_dryer.R
# Deps: deSolve (mixer ODE). Everything else is base R.
# =============================================================================

# ---- upstream: pulled UP1 mixer + shared stream interface --------------------
source("unified/up1_mixer_module.R")
source("unified/interface_stream.R")

# ---- foam-wash column module (algebraic closure from standalone ODE) ---------
source("foam_wash_module.R")

# ---- mixer nominal inputs + pars builder (from unified_model.R definitions) --
# We only need the factor dictionary, nominal_x, equipment and up1_pars_from_x;
# stop before the demo/Morris section so nothing runs on source.
.uL   <- readLines("unified_model.R")
.uCut <- grep("^# 4\\. NOMINAL END-TO-END", .uL)[1]
mx_env <- new.env()
eval(parse(text = paste(.uL[seq_len(.uCut - 1)], collapse = "\n")), envir = mx_env)
mixer_nominal_x  <- mx_env$nominal_x
up1_pars_from_x  <- mx_env$up1_pars_from_x
mx_equipment     <- mx_env$equipment
mx_factors       <- mx_env$factors

# ---- my evolved UP3 (separator) + UP4 (dryer) DEFINITIONS, into isolated envs -----------
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
cen_factors              <- cen_env$factors
spray_dry_model          <- spray_env$spray_dry_model
spray_factors            <- spray_env$factors

# spray operating / formulation setpoints: midpoints, then the plant setpoints
# (same as chained_centrifuge_spray.R).
sp_mid <- setNames((spray_factors$min + spray_factors$max) / 2, spray_factors$name)
sp_mid["T_system"]     <- 423   # dryer inlet ~150 C
sp_mid["T_feed"]       <- 300   # atomizing air ~ ambient
sp_mid["T_sticky_K"]   <- 493   # melt/clumping ~220 C (heat-resistant product)
sp_mid["mdot_gas_dry"] <- 0.25  # lands outlet ~80-85 C
sp_from_centrifuge <- c("rho_L", "C_solid_mass", "alpha_g_0", "sigma", "D_b", "mu_L",
                        "C_surfactant", "C_monomer", "C_plasticizer", "C_binder",
                        "I_strength", "Delta_pH")

# =============================================================================
# NOTE: foam_wash_column() is now sourced from foam_wash_module.R
# It implements hindered settling, bimodal bubble dynamics, surfactant film
# elasticity, and per-class particle loss via algebraic closure.
# =============================================================================
# ADAPTER: mixer-exit stream  ->  UP3 (separator) input factor vector
# =============================================================================
# Overwrites ONLY the UP3 inputs the upstream physically determines;
# the UP3's own operating knobs (rpm, scroll, pond, beach, pop_frac,
# S_base/ceiling, ...) stay at cen_op (nominal, or whatever the sweep sets).
#
# floc_strength_Pa: the aggregates the UP3 tears apart were BUILT by the
# mixer, so the mixer's Bond_Strength drives the UP3's floc strength.
# Anchored so the mixer-nominal Bond_Strength maps to 2000 Pa (the plant-
# calibrated value), then clamped to the UP3 factor range.
.bond_ref <- local({
  r1 <- up1_run_mixer(up1_pars_from_x(mixer_nominal_x), mx_equipment)
  unname(r1$outputs[["Bond_Strength"]])
})
.floc_ref  <- cen_nominal[["floc_strength_Pa"]]              # 2000 Pa (calibrated)
.floc_lo   <- cen_factors$lo[cen_factors$name == "floc_strength_Pa"]
.floc_hi   <- cen_factors$hi[cen_factors$name == "floc_strength_Pa"]

stream_to_centrifuge <- function(stream, cen_op = cen_nominal, couple_floc = TRUE) {
  run <- cen_op
  run["feed_solid_frac"] <- stream$C_solid
  run["feed_gas_frac"]   <- stream$alpha_g                  # already foam-washed
  run["plasticizer_frac"]<- stream$C_plasticizer
  run["C_monomer"]       <- stream$C_monomer
  run["C_binder"]        <- stream$C_binder
  run["I_strength"]      <- stream$ionic_strength
  run["Delta_pH"]        <- stream$Delta_pH
  run["HLB_value"]       <- stream$HLB
  run["T_process"]       <- stream$T_K
  run["P_feed"]          <- stream$P_Pa
  run["surf_dose_kg_m3"] <- stream$C_surfactant * stream$rho_slurry   # wt/wt -> kg/m3
  if (is.finite(stream$D_b_m)) run["D_b"] <- stream$D_b_m
  if (couple_floc && is.finite(.bond_ref) && .bond_ref > 0) {
    fs <- .floc_ref * stream$Bond_Strength / .bond_ref
    run["floc_strength_Pa"] <- min(max(fs, .floc_lo), .floc_hi)
  }
  run
}

# =============================================================================
# THE FULL TRAIN
# =============================================================================
run_full_train <- function(mixer_x = mixer_nominal_x, template_type = 4,
                           wash_pars = list(), cen_op = cen_nominal,
                           reslurry_solids = 0.30, spray_op = sp_mid,
                           verbose = FALSE) {
  eq <- mx_equipment; eq$template_type <- template_type

  # 1. mixer
  r1 <- up1_run_mixer(up1_pars_from_x(mixer_x), eq)
  if (any(is.na(r1$outputs))) return(NULL)
  s <- stream_from_up1(r1, up1_pars_from_x(mixer_x), eq)
  if (verbose) print_stream(s, "1) mixer exit")

  # 2. pressurized foam-wash column
  s <- foam_wash_column(s, wash_pars)
  if (verbose) print_stream(s, "2) UP2 foam-wash exit -> UP3 feed")

  # 3. UP3 (separator) + 4. reslurry handoff (one call: runs the model, dilutes)
  cen_run <- stream_to_centrifuge(s, cen_op)
  cen_out <- unified_centrifuge_model(cen_run)
  hand    <- centrifuge_to_spray(cen_run, target_solid_mass = reslurry_solids)

  # 5. UP4 (dryer)
  x <- spray_op
  for (nm in sp_from_centrifuge) x[nm] <- hand[[nm]]
  sp <- spray_dry_model(x)

  list(mixer = r1$outputs, stream = s, cen_run = cen_run,
       cen_out = cen_out, handoff = hand, spray = sp, spray_in = x)
}

# =============================================================================
# NOMINAL END-TO-END RUN (prints the stream at every interface)
# =============================================================================
cat("=== FULL TRAIN: UP1 -> UP2 (foam-wash) -> UP3 (separator) -> reslurry -> UP4 (dryer) ===\n")
cat("    (nominal inputs, capillary_bridge template)\n\n")
res <- run_full_train(verbose = TRUE)

m <- res$mixer; s <- res$stream; co <- res$cen_out; h <- res$handoff; sp <- res$spray
cat(sprintf("\n[UP1]            C_solid %.3f  alpha_g %.3f  D_agg %.0f um  Bond %.3f  RTF %.3f\n",
            s$C_solid, m[["Blended_Porosity"]], m[["Blended_Size_um"]],
            m[["Bond_Strength"]], m[["Residual_Template_Fraction"]]))
cat(sprintf("[UP2 foam-wash]  alpha_g %.3f -> %.3f (washed)  D_b %.1f um  P %.2f atm\n",
            m[["Blended_Porosity"]], s$alpha_g, s$D_b_m*1e6, s$P_Pa/1.013e5))
cat(sprintf("[UP3 separator]  cake solids %.1f%%  exit dens %.2f g/cc  gas holdup %.3f  floc_used %.0f Pa\n",
            co$Product_Solids_MassFrac*100, co$Exit_Density_kg_m3/1000,
            co$Entrained_Gas_Holdup, res$cen_run[["floc_strength_Pa"]]))
cat(sprintf("[reslurry 30%%]   rho_L %.0f  C_solid %.2f  alpha_g0 %.3f  mu_L %.4f  dil x%.2f\n",
            h[["rho_L"]], h[["C_solid_mass"]], h[["alpha_g_0"]], h[["mu_L"]], h[["dilution_x"]]))
cat(sprintf("[UP4 dryer]      D_particle %.1f um  porosity %.3f  skin %.3f  rho_tap %.0f  X_moist %.3f\n",
            sp[["D_particle_um"]], sp[["phi_porosity_z"]], sp[["theta_skin_z"]],
            sp[["rho_tapped"]], sp[["X_moisture"]]))

cat("\n=== UP2 (foam-wash) sensitivity: how much gas removal changes the train ===\n")
cat(sprintf("  %-10s %-12s %-12s %-12s %-12s\n",
            "eta_gas", "cen_gas_hld", "reslurry_ag0", "D_particle", "porosity"))
for (eg in c(0.0, 0.50, 0.75, 0.95)) {
  r <- run_full_train(wash_pars = list(eta_gas = eg))
  cat(sprintf("  %-10.2f %-12.3f %-12.3f %-12.1f %-12.3f\n", eg,
              r$cen_out$Entrained_Gas_Holdup, r$handoff[["alpha_g_0"]],
              r$spray[["D_particle_um"]], r$spray[["phi_porosity_z"]]))
}
cat("\nMore foam washed out -> less gas into the UP3/UP4 stages -> denser,\n",
    "less porous powder (and, per your plant note, easier to dry).\n", sep = "")
