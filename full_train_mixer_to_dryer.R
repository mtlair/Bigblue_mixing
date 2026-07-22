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

# ---- UP2 foam-wash column: the 491-line ODE (functions only; skip its driver
#      + plots by cutting at the "solve + report" marker) --------------------
fc_env <- new.env()
eval(parse(text = src_defs("foam_wash_column_psd.R", "^# --- solve \\+ report")), envir = fc_env)

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
# UP2 (foam-wash) as the 491-line COLUMN ODE  —  stream-contract wrapper
# =============================================================================
# The ODE (foam_wash_column_psd.R) is a standalone script, not an f(stream)->
# stream stage. This wrapper maps the mixer-exit stream onto the column feed,
# solves the column, and writes back the SAME fields the algebraic placeholder
# touched (alpha_g, D_b_m, P_Pa) — now driven by the ODE instead of a fixed
# eta_gas. Couplings that are unambiguous are taken from the stream; the
# calibrated film-chemistry / geometry knobs keep their column defaults.
#
# Stream -> column feed:
#   rho_liquid <- rho_slurry;  rho_p <- rho_polymer
#   T_col      <- T_K   (T_ref set equal so mu(T_col) == mu_ref, no Andrade offset)
#   P_col      <- P_Pa  (if finite)
#   eps_g_in   <- alpha_g (mixer trapped-gas holdup)         } renormalised
#   eps_s_in   <- phi_s = C_solid*rho_L/rho_polymer          } to sum to 1
#   eps_l_in   <- 1 - eps_g_in - eps_s_in                    }
#   c_surf,cmc <- C_surfactant, CMC  (wt/wt -> mol/m3 via rho_L / MW_surfactant)
#   V_tip      <- mixer tip speed (Hinze born-bubble size), if passed in pars
# Column -> stream (writeback):
#   alpha_g <- alpha_g * (1 - f_slug)  (the coalesced-slug gas — coarse bubbles
#              that grow past the film-rupture size and burst — disengages
#              overhead and is removed; the still-dispersed foam stays entrained.
#              Weaker/less-elastic films burst sooner -> more slug -> stronger
#              wash. f_slug = J_g_slug / J_g_in at the top.  NOTE: the raw ODE
#              *conserves* gas out the top; this slug split is the interpretation
#              used to give the train a gas-removal number, and is the responsive
#              (film-chemistry-coupled) part of the gas balance — the fine mode
#              never bursts, so it carries a fixed share and cannot drive a wash.)
#   D_b_m   <- top FINE-mode bubble diameter (the dispersed population that stays
#              entrained downstream; coarse/slug disengage). This is a foam-column
#              bubble, coarser than the dryer's native D_b factor range — see the
#              caveat in the report.
#   P_Pa    <- column operating pressure
# Diagnostics (per-class retention, impurity, film_stability, derived eta_gas,
# gas-balance closure) are attached as attr(stream,"up2_ode").
foam_wash_column_ode <- function(stream, pars = list()) {
  P <- fc_env$params
  rho_L <- stream$rho_slurry
  P[["rho_liquid"]] <- rho_L
  P[["rho_p"]]      <- stream$rho_polymer
  P[["T_col"]]      <- stream$T_K
  P[["T_ref"]]      <- 298.15                     # let Andrade scale mu to T_col
  if (is.finite(stream$P_Pa) && stream$P_Pa > 0) P[["P_col"]] <- stream$P_Pa

  # Tie the continuous-phase (serum) viscosity to LOCAL composition + temperature:
  # water * binder thickening (matches the dryer's mu_serum = mu_L*(1+10 C_bind)),
  # then the column's Andrade mu(T) scales it to the column temperature. This
  # feeds settling, film drainage and bubble burst; the surfactant conversion
  # below sets film elasticity (sigma / Gibbs E) and the Hinze born-bubble size.
  P[["mu_ref"]] <- 1.0e-3 * (1 + 10 * max(stream$C_binder, 0))

  # holdups (renormalised to close the phase balance)
  C_sol <- stream$C_solid + stream$C_solid_rigid
  phi_s <- max(1e-4, min(0.60, C_sol * rho_L / stream$rho_polymer))
  eg    <- max(1e-4, min(0.90, stream$alpha_g))
  el    <- max(0.05, 1 - eg - phi_s)
  tot   <- eg + phi_s + el
  P[["eps_g_in"]] <- eg / tot
  P[["eps_s_in"]] <- phi_s / tot
  P[["eps_l_in"]] <- el / tot

  # surfactant wt/wt -> mol/m3
  MW_kg <- max(stream$MW_surfactant / 1000, 1e-3)
  P[["c_surf"]] <- max(0,    stream$C_surfactant * rho_L / MW_kg)
  P[["cmc"]]    <- max(1e-6, stream$CMC          * rho_L / MW_kg)

  # UP2 wash-liquid feed. wash_ratio (wash liquid / foam, a native column knob)
  # drives the column's internal impurity wash; wash_carry is the fraction of
  # that wash liquid that leaves OVERHEAD WITH THE FOAM into UP3 (co- vs
  # counter-current), diluting the stream handed downstream.
  wash_carry <- if (!is.null(pars$wash_carry)) max(0, min(1, pars$wash_carry)) else 0.30

  if (!is.null(pars$v_tip)) P[["V_tip"]] <- pars$v_tip
  # explicit overrides (e.g. a sensitivity sweep over a real column knob)
  for (nm in names(pars)) if (nm %in% names(P)) P[[nm]] <- pars[[nm]]

  dpar <- fc_env$derive_state_props(P)
  out  <- fc_env$solve_column(dpar)
  top  <- out[nrow(out), ]
  Jg_in <- dpar[["eps_g_in"]] * dpar[["U_up"]]

  f_slug <- min(max(as.numeric(top$J_g_slug) / Jg_in, 0), 1)

  # wash liquid joining the foam overhead dilutes the slurry going to UP3
  dil_w <- 1 + wash_carry * P[["wash_ratio"]]
  for (f in c("C_solid", "C_binder", "C_monomer", "C_plasticizer",
              "C_surfactant", "ionic_strength"))
    stream[[f]] <- stream[[f]] / dil_w
  stream$rho_slurry <- 1000 + (rho_L - 1000) / dil_w   # blend toward wash water

  stream$alpha_g <- max(0, stream$alpha_g * (1 - f_slug))
  # Entrained-microbubble size handed to the atomizer — NOT the mm foam bubble
  # the column tracks (that is the foam being washed). The fine bubbles that
  # survive in the feed liquid are set by mixer shear (Hinze: d ~ V_tip^-1.2) and
  # film surface tension, calibrated to ~50 um at the reference tip speed. This
  # is the population that gates effervescent atomization at the nozzle, so v_tip
  # (and the film chemistry that sets sigma) reaches the droplet/particle size.
  D_MICRO_REF <- 45e-6
  D_micro <- D_MICRO_REF *
             (P[["V_tip_ref"]] / max(P[["V_tip"]], 0.1))^P[["n_hinze"]] *
             (dpar[["sigma_film"]] / P[["sigma_ref_film"]])^0.6
  stream$D_b_m   <- min(max(D_micro, 5e-6), 5e-4)   # clamp to atomizer range
  stream$P_Pa    <- unname(dpar[["P_col"]])

  attr(stream, "up2_ode") <- list(
    eta_gas        = f_slug,
    d_b_fine_um    = as.numeric(top$d_b_fine)   * 1e6,
    d_b_coarse_um  = as.numeric(top$d_b_coarse) * 1e6,
    film_stability = unname(dpar[["film_stability"]]),
    sigma_mNm      = unname(dpar[["sigma_film"]]) * 1e3,
    mu_serum_mPas  = unname(P[["mu_ref"]]) * 1e3,
    wash_ratio     = unname(P[["wash_ratio"]]),
    dil_wash       = dil_w,
    ret_fine = as.numeric(top$Js_fine / out$Js_fine[1]),
    ret_mid  = as.numeric(top$Js_mid  / out$Js_mid[1]),
    ret_crs  = as.numeric(top$Js_crs  / out$Js_crs[1]),
    impurity_out = as.numeric(top$C_imp),
    gas_closure  = as.numeric((top$J_g_foam_fine + top$J_g_foam_coarse +
                               top$J_g_slug - Jg_in) / Jg_in))
  stream
}

# =============================================================================
# THE FULL TRAIN
# =============================================================================
run_full_train <- function(mixer_x = mixer_nominal_x, template_type = 4,
                           wash_pars = list(), cen_op = cen_nominal,
                           reslurry_add = 0, spray_op = sp_mid,
                           up2 = c("ode", "algebraic"), verbose = FALSE) {
  up2 <- match.arg(up2)
  eq <- mx_equipment; eq$template_type <- template_type
  p1 <- up1_pars_from_x(mixer_x)

  # 1. mixer
  r1 <- up1_run_mixer(p1, eq)
  if (any(is.na(r1$outputs))) return(NULL)
  s <- stream_from_up1(r1, p1, eq)
  if (verbose) print_stream(s, "1) mixer exit")

  # 2. pressurized foam-wash column (UP2)
  if (up2 == "ode") {
    s <- foam_wash_column_ode(s, modifyList(list(v_tip = unname(p1$v_tip)), wash_pars))
  } else {
    s <- foam_wash_column(s, wash_pars)   # algebraic placeholder (foam_wash_module.R)
  }
  if (verbose) print_stream(s, "2) UP2 foam-wash exit -> UP3 feed")

  # 3. UP3 (separator) + 4. reslurry handoff (one call: runs the model, dilutes)
  cen_run <- stream_to_centrifuge(s, cen_op)
  cen_out <- unified_centrifuge_model(cen_run)
  hand    <- centrifuge_to_spray(cen_run, reslurry_add = reslurry_add)

  # 5. UP4 (dryer)
  x <- spray_op
  for (nm in sp_from_centrifuge) x[nm] <- hand[[nm]]
  # UP1-bounded template + feed temperature + v_tip carried through the stream (no
  # pre-heater stage yet, so the dryer feed T is the mixer-exit T; the template
  # emulsion fraction and droplet size are the mixer's, not free dryer knobs).
  x["phi_emulsion"]      <- s$phi_templ_free
  x["D_template"]        <- s$D_template_um * 1e-6
  x["T_feed"]            <- s$T_K
  x["D_primary_exit_um"] <- s$D_primary_exit_um  # UP1 three-regime primary size
  x["D_agg_um"]         <- s$D_agg_um           # UP1 calibrated aggregate d50
  x["D_primary_phys_um"] <- s$D_primary_phys_um # physical colloid bead (200 nm base)
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
d2 <- attr(s, "up2_ode")
cat(sprintf("  (491-ODE)      eta_gas %.2f  film_stab %.2f  sigma %.1f mN/m  retention f/m/c %.0f/%.0f/%.0f%%  imp->%.0f%%  gas-closure %.1e\n",
            d2$eta_gas, d2$film_stability, d2$sigma_mNm,
            100*d2$ret_fine, 100*d2$ret_mid, 100*d2$ret_crs, d2$impurity_out, d2$gas_closure))
cat(sprintf("[UP3 hi-g sep]  cake solids %.1f%%  exit dens %.2f g/cc  gas holdup %.3f  floc_used %.0f Pa\n",
            co$Product_Solids_MassFrac*100, co$Exit_Density_kg_m3/1000,
            co$Entrained_Gas_Holdup, res$cen_run[["floc_strength_Pa"]]))
cat(sprintf("[feed->UP4]      rho_L %.0f  C_solid %.2f (UP3-tied)  alpha_g0 %.3f  mu_L %.4f  reslurry x%.2f\n",
            h[["rho_L"]], h[["C_solid_mass"]], h[["alpha_g_0"]], h[["mu_L"]], h[["dilution_x"]]))
cat(sprintf("[UP4 dryer]      D_particle %.1f um  porosity %.3f  skin %.3f  rho_tap %.0f  X_moist %.3f\n",
            sp[["D_particle_um"]], sp[["phi_porosity_z"]], sp[["theta_skin_z"]],
            sp[["rho_tapped"]], sp[["X_moisture"]]))

cat("\n=== UP2 (491-ODE) sensitivity: surfactant dose -> film stability -> train ===\n")
cat(sprintf("  %-12s %-10s %-12s %-12s %-12s %-12s\n",
            "c_surf", "eta_gas", "cen_gas_hld", "reslurry_ag0", "D_particle", "porosity"))
for (cs in c(1.0, 3.0, 5.0, 8.0)) {
  r  <- run_full_train(wash_pars = list(c_surf = cs))
  d2 <- attr(r$stream, "up2_ode")
  cat(sprintf("  %-12.1f %-10.2f %-12.3f %-12.3f %-12.1f %-12.3f\n", cs,
              d2$eta_gas, r$cen_out$Entrained_Gas_Holdup, r$handoff[["alpha_g_0"]],
              r$spray[["D_particle_um"]], r$spray[["phi_porosity_z"]]))
}
cat("\nMore surfactant -> more-stable films -> less bubble coalescence -> more gas\n",
    "stays finely dispersed (weaker wash, lower eta_gas) -> gassier feed into\n",
    "UP3/UP4 -> more porous powder. Less surfactant washes harder (denser powder).\n", sep = "")
