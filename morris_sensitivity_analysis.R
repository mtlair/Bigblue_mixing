#!/usr/bin/env Rscript
# =============================================================================
# Morris (Elementary Effects) sensitivity analysis of the reduced-order
# two-fluid (internal-mix) atomization + drying model described in
# "deepresearchreport.md", extended with closures from the Latex Coagulation
# Engine Integrated Master Specification (v30.0.0 / R-Engine v47):
#
#   * Flory-Huggins free-volume swelling: residual solvent (C_monomer +
#     C_plasticizer) softens the matrix, Softness = (1+5 C_binder) exp(25 phi_solvent)
#   * Product glass transition: Fox equation on residual solvent -> Tg_eff,
#     driving stickiness / pore collapse / caking above Tg
#   * Surfactant molar stoichiometry: theta_surf from C_surfactant, MW and
#     A_molecule vs colloid + bubble surface area
#   * DLVO electrostatics: E_repulsion = dpH * Debye screening * (1-theta_surf),
#     gating aggregation (skin onset) and fractal openness (D_f -> porosity)
#   * Krieger-Dougherty crowding viscosity from feed solids (phi_s)
#
# Spray-specific extensions beyond the v47 reactor spec:
#   * Pressurized-hold bubble coarsening: Ostwald ripening + coalescence of
#     the entrained foam during the pressurized residence time t_hold before
#     the nozzle (LSW r^3 ~ t kinetics, Henry's-law accelerated by pressure,
#     damped by surfactant film coverage)
#   * Effervescent-airblast hybrid exit plane: nothing atomizes in the feed
#     line (confined flow); the entrained gas acts AT the bi-fluid exit,
#     where its throat voidage thins the annular liquid film (reducing the
#     characteristic length fed to the airblast correlation) and bubbles
#     small enough to ride inside the films flash-shatter the fragments on
#     the final letdown (multiplicative, gated by bubble size vs film
#     thickness). The nozzle carries its own pressure ladder: liquid line
#     P_feed -> mixing chamber P_chamber = 0.8 P_system -> near-choked
#     throat -> ambient
#   * Dryer energy/moisture balance: dryer gas flow rate mdot_gas_dry and
#     inlet absolute humidity Y_in give co-current outlet temperature and
#     relative humidity, which drive the drying Peclet number and the
#     sticky-point state (replaces the former independent RH_gas knob)
#   * Shell permeability closure: residual monomer and plasticizer open
#     free-volume diffusion pathways through the forming skin (permeability
#     rises exponentially with solvent load) while binder films plug
#     inter-particle pores; the relative permeability gates the falling-rate
#     drying resistance and the vapor entrapment (vacuole inflation /
#     blowhole rupture) of skinned droplets
#   * Immiscible template-solvent emulsion: low-boiling (T_bp < T_bp,water)
#     solvent droplets, pre-processed to a target size, ride through the
#     nozzle incompressibly (no expansion work; only superheat flash at the
#     letdown assists atomization) and survive early drying because their
#     escape is solubility-limited. Once the particle crosses the solvent
#     boiling point the fate splits three ways: clean D_e-sized templated
#     pores (permeable shell), balloon inflation (vapor trapped under the
#     skin), or micro-explosion bursting (Clausius-Clapeyron overpressure
#     beats the cake yield stress). The sparged-gas foam chain (holdup,
#     hold-time ripening, effervescent exit, per-mode gas trapping) is
#     retained unchanged and acts in parallel.
#
# NOTE: v_tip does not feed UP4 directly. All colloid effects of tip speed
# (three regimes: unchanged / aggregation onset / milling+attrition) are
# resolved in UP1. The result flows downstream as D_primary_exit_um on the
# stream; UP4 uses that size as-is via Module 0a.
#
# Module chain: Feed Properties -> Pressurized Hold -> Two-Phase
# Conditioning -> Effervescent Stage -> Bi-Fluid Airblast Stage ->
# Secondary Breakup -> Dryer Balance -> Particle Formation.
#
# Outputs screened (nomenclature-sheet symbols where they exist):
#   1. D_particle     - final (dry) particle size                    [um]
#   2. theta_skin_z   - skin network fraction (skin formation)       [-]
#   3. Omega_struct_z - structural memory / sphericity state         [-]
#   4. phi_porosity_z - total porosity (void fraction state)         [-]
#   5. rho_tapped     - powder tapped density (bulk analogue of
#                       rho_colloid_out / SG_out)                    [kg/m3]
#   6. Tg_product     - effective product glass transition T_g,eff   [K]
#
# Method: Morris one-at-a-time screening. Uses the `sensitivity` package
# (morris()) when installed; otherwise falls back to an internal base-R
# implementation of the same OAT trajectory design, so the script has no
# hard package dependencies. The Morris plot (mu* vs sigma) is drawn with
# base graphics.
#
# Run:  Rscript morris_sensitivity_analysis.R
# =============================================================================

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Input factors and ranges
# -----------------------------------------------------------------------------
# One row per factor: name, min, max, log (TRUE = sampled log-uniformly over
# its range, used for wide-ranging positive factors), unit, and description.
# Symbols follow the master nomenclature sheet where they exist; others use
# the spray report's notation. t_hold, mdot_gas_dry and Y_in are spray-line
# additions (pressurized residence before the nozzle; dryer gas flow and
# inlet humidity); k_perm_* are the shell permeability coefficients for
# monomer / plasticizer / binder.
fac <- function(name, min, max, log, unit, desc)
  data.frame(name = name, min = min, max = max, log = log,
             unit = unit, desc = desc, stringsAsFactors = FALSE)

factors <- rbind(
  #   name             min     max     log    unit     description
  fac("ALR",           0.9,    1.8,    FALSE, "-",     "air-liquid mass ratio m_G/m_L (measured up4atom_scfm/up4_feed, visc.xlsx: 0.90-1.78)"),
  fac("P_system",      1.2e5,  1.8e5,  FALSE, "Pa",    "atomizing air supply pressure (measured max_up4atom_psig, visc.xlsx: 2.7-11.4 psig abs)"),
  fac("P_feed",        1.5e5,  1.0e6,  FALSE, "Pa",    "liquid feed line (hold) pressure"),
  fac("mdot_L",        0.013,  0.020,  FALSE, "kg/s",  "liquid feed mass flow (measured up4_feed, visc.xlsx: 107-154 lb/hr = 0.0135-0.0194 kg/s)"),
  fac("sigma",         0.030,  0.070,  FALSE, "N/m",   "liquid surface tension"),
  fac("mu_L",          0.0012, 0.0560, TRUE,  "Pa s",  "serum (continuous phase) viscosity"),
  fac("rho_L",         1000,   1300,   FALSE, "kg/m3", "liquid (slurry) density"),
  fac("alpha_g_0",     0.05,   0.60,   FALSE, "-",     "feed foam quality / entrained gas holdup"),
  fac("D_b",           2.0e-5, 2.0e-4, TRUE,  "m",     "feed bubble diameter (before shear)"),
  fac("C_solid_mass",  0.05,   0.40,   FALSE, "-",     "solid mass fraction in feed (wt/wt)"),
  fac("T_system",      330,    470,    FALSE, "K",     "dryer gas INLET temperature"),
  fac("mdot_gas_dry",  0.10,   1.00,   TRUE,  "kg/s",  "dryer gas mass flow rate"),
  fac("Y_in",          0.001,  0.020,  FALSE, "kg/kg", "dryer gas inlet absolute humidity"),
  fac("T_feed",        280,    330,    FALSE, "K",     "feed / atomizing air temperature"),
  fac("t_hold",        5,      600,    TRUE,  "s",     "pressurized hold time before nozzle"),
  fac("C_monomer",     0.000,  0.020,  FALSE, "-",     "residual monomer concentration (wt/wt)"),
  fac("C_plasticizer", 0.000,  0.050,  FALSE, "-",     "plasticizer concentration (wt/wt)"),
  fac("C_binder",      0.000,  0.050,  FALSE, "-",     "binder concentration (wt/wt)"),
  fac("C_surfactant",  1.0e-4, 2.0e-2, TRUE,  "-",     "formulated surfactant concentration (wt/wt)"),
  fac("Delta_pH",      0.2,    4.0,    FALSE, "pH",    "delta pH vs isoelectric point"),
  fac("I_strength",    1.0e-3, 5.0e-1, TRUE,  "M",     "ionic strength (Debye screening driver)"),
  fac("Tg_polymer",    280,    380,    FALSE, "K",     "dry-polymer glass transition"),
  fac("n_flow",        0.20,   0.67,   FALSE, "-",     "power-law shear-thinning flow index (measured post-UP1 envelope, visc.xlsx: 0.20-0.67, mean 0.44)"),
  fac("k_perm_mono",   5.0,    60.0,   FALSE, "-",     "monomer free-volume permeation coeff."),
  fac("k_perm_plast",  5.0,    60.0,   FALSE, "-",     "plasticizer free-volume permeation coeff."),
  fac("k_perm_bind",   2.0,    40.0,   FALSE, "-",     "binder pore-blocking coefficient"),
  fac("phi_emulsion",  0.00,   0.30,   FALSE, "-",     "template solvent emulsion volume fraction"),
  fac("D_template",    5.0e-7, 1.0e-5, TRUE,  "m",     "pre-processed emulsion droplet diameter"),
  fac("T_bp_solv",     300,    360,    FALSE, "K",     "template solvent boiling point (ambient)")
)
k <- nrow(factors)

# Map the unit-hypercube Morris design onto physical values
scale_design <- function(X01) {
  X <- X01
  for (j in seq_len(k)) {
    if (factors$log[j]) {
      X[, j] <- exp(log(factors$min[j]) +
                    X01[, j] * (log(factors$max[j]) - log(factors$min[j])))
    } else {
      X[, j] <- factors$min[j] + X01[, j] * (factors$max[j] - factors$min[j])
    }
  }
  colnames(X) <- factors$name
  X
}

# -----------------------------------------------------------------------------
# 2. Built-in reduced-order model (one row of physical inputs -> 6 outputs)
# -----------------------------------------------------------------------------
# Fixed constants / nozzle geometry (report worked example + v47 defaults)
P_atm    <- 1.013e5   # ambient pressure                        [Pa]
R_air    <- 287       # gas constant, air                       [J/kg K]
gamma_a  <- 1.4       # heat capacity ratio, air                [-]
A_L      <- 1.0e-6    # liquid passage area                     [m2]
A_G      <- 5.0e-6    # gas passage area                        [m2]
D_h      <- 1.0e-3    # nozzle hydraulic diameter               [m]
rho_s    <- 1400      # dry solid (polymer) density             [kg/m3]
a_prim   <- 1.0e-7    # primary colloid particle radius         [m]
kB       <- 1.380649e-23
N_Av     <- 6.02214e23
MW_surf  <- 0.400     # surfactant molecular weight             [kg/mol]
A_molec  <- 0.5e-18   # surfactant molecule area capacity       [m2]
HLB      <- 12        # surfactant HLB (foam stability driver)  [-]
phi_m    <- 0.63      # Krieger-Dougherty max packing           [-]
Tg_solv  <- 150       # solvent (monomer/plasticizer) Tg, Fox   [K]
d_ratio  <- 0.10      # primary/aggregate size ratio: fallback constant [-]
                      # (overridden inside spray_dry_model when D_agg_um is wired)
C_cham   <- 0.80      # mixing chamber / air supply pressure    [-]
k_rip0   <- 4.0e-15   # Ostwald ripening rate at 1 atm          [m3/s]
cp_gas   <- 1005      # dryer gas heat capacity                 [J/kg K]
h_fg     <- 2.30e6    # latent heat of water evaporation        [J/kg]
gamma_ref <- 100      # rheometer reference shear rate          [1/s]
t_res    <- 8         # dryer residence time                    [s]
Tg_water <- 138       # glass transition of water (Fox)         [K]
sig_y0   <- 4.0e6     # cake yield strength scale (Rumpf-type)  [Pa]
rho_solv <- 750       # template solvent liquid density         [kg/m3]
h_fg_solv <- 3.5e5    # template solvent latent heat            [J/kg]
k_emu0   <- 2.0e-20   # emulsion LSW ripening rate (solubility-
                      # limited; immiscible in water)           [m3/s]

spray_dry_model <- function(x) {
  ALR    <- x[["ALR"]];        P_G    <- x[["P_system"]]
  P_F    <- x[["P_feed"]]
  mdot_L <- x[["mdot_L"]];     sigma  <- x[["sigma"]]
  mu_L   <- x[["mu_L"]];       rho_L  <- x[["rho_L"]]
  alpha0 <- x[["alpha_g_0"]];  D_b    <- x[["D_b"]]
  C_sol  <- x[["C_solid_mass"]]
  T_in   <- x[["T_system"]];   mdot_g <- x[["mdot_gas_dry"]]
  Y_in   <- x[["Y_in"]];       T_feed <- x[["T_feed"]]
  t_hold <- x[["t_hold"]]
  C_mono <- x[["C_monomer"]];  C_plas <- x[["C_plasticizer"]]
  C_bind <- x[["C_binder"]];   C_surf <- x[["C_surfactant"]]
  dpH    <- x[["Delta_pH"]];   I_str  <- x[["I_strength"]]
  Tg_pol <- x[["Tg_polymer"]]; n_fl   <- x[["n_flow"]]
  k_pm   <- x[["k_perm_mono"]]
  k_pp   <- x[["k_perm_plast"]]
  k_pb   <- x[["k_perm_bind"]]
  phi_e  <- x[["phi_emulsion"]]
  D_e    <- x[["D_template"]]
  T_bpS  <- x[["T_bp_solv"]]
  ## --- Module 0a: Primary colloid radius (from UP1 three-regime closure) -----
  # All v_tip effects on the colloid (milling, aggregation) complete in UP1.
  # D_primary_exit_um carries the result: unchanged at low v_tip, attrite past
  # the critical tip speed. UP4 uses it as the physical primary size without
  # any v_tip awareness of its own.
  # Prefer D_primary_phys_um (200 nm physical colloid base) over D_primary_exit_um
  # (ODE-scaled, 1.25 µm base) for packing/surfactant/Pe calculations.
  a_prim_mod <- if ("D_primary_phys_um" %in% names(x) && !is.na(x[["D_primary_phys_um"]]))
                  max(0.01e-6, x[["D_primary_phys_um"]] * 1e-6 / 2)
                else if ("D_primary_exit_um" %in% names(x) && !is.na(x[["D_primary_exit_um"]]))
                  max(0.01e-6, x[["D_primary_exit_um"]] * 1e-6 / 2)
                else a_prim     # fall back to module constant if stream not wired

  # Physical d_ratio: D_primary_exit / D_agg_phys, wired from UP1 stream.
  # Calibration: D_agg ≈ 9.6 µm at v_tip = 8.5 m/s (just above critical 6.81 m/s).
  # Falls back to module constant 0.10 when the stream field is absent.
  d_ratio_stream <- if ("D_agg_um" %in% names(x) && !is.na(x[["D_agg_um"]]) &&
                        x[["D_agg_um"]] > 0)
                      min(max(a_prim_mod * 2e6 / x[["D_agg_um"]], 0.01), 1.0)
                    else d_ratio
  d_ratio <- d_ratio_stream

  ## --- Module 0b: Formulation - Flory-Huggins free-volume swelling ---------
  # (v47 Sect. 4) residual solvent acts as theta-solvent proxy, softening the
  # matrix exponentially and independently of temperature
  phi_solvent <- C_mono + C_plas
  Softness    <- (1 + 5 * C_bind) * exp(25 * phi_solvent)

  # Shell permeability (free-volume theory): residual monomer and
  # plasticizer swell the polymer shell and open diffusive pathways for
  # water vapor - permeability grows exponentially with each solvent load,
  # with species-specific permeation coefficients. Binder films spread over
  # the primary particles and plug the inter-particle pore throats
  # (hydraulic blocking). Perm_shell = 1 for a clean formulation
  # (no monomer, plasticizer, or binder present).
  Perm_shell <- exp(k_pm * C_mono + k_pp * C_plas) / (1 + k_pb * C_bind)

  # Template solvent emulsion: immiscible, incompressible liquid droplets
  # pre-processed to target size D_e; they occupy liquid volume without
  # contributing solids, and their mass rides the feed until the dryer
  rho_lm <- (1 - phi_e) * rho_L + phi_e * rho_solv  # liquid mixture density
  w_e    <- phi_e * rho_solv / rho_lm               # solvent mass fraction

  ## --- Module 0c: Feed rheology (Krieger-Dougherty + power-law thinning) ---
  phi_s      <- C_sol * rho_L / rho_s             # solid phase volume fraction
  phi_disp   <- phi_s * (1 - phi_e) + phi_e       # solids + emulsion droplets
  mu_serum   <- mu_L * (1 + 10 * C_bind)          # binder thickening (K_fluid)
  mu_slurry0 <- mu_serum *
                (1 - min(phi_disp, 0.60) / phi_m)^(-2.5 * phi_m)
  # Override: when the UP1 mixer exit slurry viscosity is wired in, use it
  # directly as the power-law reference, bypassing the serum-based KD build-up.
  # mu_slurry_up1 is the full apparent slurry viscosity at the mixer exit shear
  # rate (~100-1000 s⁻¹), so it already embeds suspension + floc effects.
  if ("mu_slurry_up1" %in% names(x) && !is.na(x[["mu_slurry_up1"]])) {
    mu_slurry0 <- x[["mu_slurry_up1"]]
  }
  # power-law slurry: mu_app(g) = K g^(n-1), anchored so mu_slurry0 is the
  # measured apparent viscosity at the rheometer reference rate gamma_ref;
  # floored at the serum high-shear limit
  mu_app <- function(gdot)
    pmax(1e-3, mu_slurry0 * (gdot / gamma_ref)^(n_fl - 1))

  ## --- Module 0d: Surfactant molar stoichiometry (v47 Sect. 5) -------------
  # monolayer capacity vs total interfacial demand: primary colloid surface
  # (3*phi_s/a per m3, the dominant sink) plus bubble surface (6*alpha/D_b)
  cap_area   <- (C_surf / MW_surf) * rho_L * N_Av * A_molec   # [m2/m3]
  part_area  <- 3 * phi_s * (1 - phi_e) / a_prim_mod          # [m2/m3]
  gas_area   <- 6 * alpha0 / D_b                              # [m2/m3]
  emu_area   <- 6 * phi_e / D_e                               # [m2/m3]
  theta_surf <- min(1, cap_area / (part_area + gas_area + emu_area))
  Foam_Stab  <- 1 + 5 * theta_surf * (HLB / 10)
  # under-stabilized foam sheds entrained gas before the nozzle
  alpha_g    <- alpha0 * (0.3 + 0.7 * theta_surf)

  ## --- Module 0e: Pressurized hold - coalescence / Ostwald ripening --------
  # LSW kinetics D^3 ~ t; Henry's-law gas solubility makes ripening faster
  # under the liquid-line pressure, surfactant coverage retards both
  k_rip <- k_rip0 * (P_F / P_atm) * (1 - 0.8 * theta_surf)
  D_b_h <- (D_b^3 + k_rip * t_hold)^(1/3)         # coarsened bubble size
  # emulsion coarsening is solubility-limited (immiscible solvent, no
  # Henry's-law pressure acceleration) and strongly retarded by surfactant
  # coverage: an under-covered emulsion drifts off its pre-processed
  # target size during the pressurized hold
  k_emu <- k_emu0 * (1 - 0.9 * theta_surf)
  D_e_h <- (D_e^3 + k_emu * t_hold)^(1/3)         # coarsened template size

  ## --- Module 0e: DLVO electrostatics (v47 Sect. 6) ------------------------
  E_rep     <- dpH * (1 / (1 + 10 * I_str)) * (1 - theta_surf)
  stability <- 1 + E_rep + 5 * theta_surf         # electrostatic + steric

  ## --- Module 1: Feed properties (at liquid-line pressure P_feed) ----------
  rho_G0    <- P_F / (R_air * T_feed)             # entrained gas density
  rho_eff0  <- (1 - alpha_g) * rho_lm + alpha_g * rho_G0  # rho_eff,z (HEM)
  mu_eff0   <- mu_slurry0 * (1 + 2.5 * alpha_g)   # low-shear mu_eff,z (bubbly)
  sigma_eff <- sigma * (1 - 0.5 * alpha_g)        # foam-reduced surface tension

  ## --- Module 2: Nozzle pressure ladder & two-phase conditioning -----------
  # liquid line P_feed -> mixing chamber -> near-choked throat -> ambient;
  # alpha_g_0 is measured at line pressure, so all expansions reference P_feed
  P_cham   <- C_cham * P_G                        # internal chamber pressure
  P_throat <- max(P_atm, 0.53 * P_cham)           # near-choked throat pressure
  r_th     <- max(P_F / P_throat, 1)              # feed -> throat expansion
  r_exp    <- max(P_F / P_atm, 1)                 # feed -> ambient expansion
  r_flash  <- max(P_throat / P_atm, 1)            # final letdown (shattering)
  alpha_th <- min(alpha_g * r_th / (1 - alpha_g + alpha_g * r_th), 0.97)
  alpha_e  <- min(alpha_g * r_exp / (1 - alpha_g + alpha_g * r_exp), 0.97)
  D_b_e    <- D_b_h * r_exp^(1/3)                 # bubble growth D_b ~ P^-1/3
  rho_Ge   <- P_atm / (R_air * T_feed)
  rho_eff  <- (1 - alpha_e) * rho_lm + alpha_e * rho_Ge

  ## --- Module 3: Nozzle hydraulics ------------------------------------------
  U_L <- mdot_L / (rho_lm * A_L)                  # liquid exit velocity
  # gas: fully-expanded isentropic velocity (choked upstream if P_G/P_atm>1.89)
  U_G <- 0.9 * sqrt(pmax(2 * gamma_a / (gamma_a - 1) * R_air * T_feed *
                         (1 - (P_atm / P_G)^((gamma_a - 1) / gamma_a)), 0))
  U_rel <- max(U_G - U_L, 5)
  rho_A <- P_atm / (R_air * T_feed)               # ambient air density

  ## --- Module 4: Effervescent-airblast hybrid exit plane -------------------
  # Confined flow cannot atomize: the entrained gas does its work AT the
  # bi-fluid exit orifice, simultaneously with the airblast shear.
  # (a) throat voidage of the entrained gas thins the annular liquid film,
  #     shrinking the characteristic length the airblast correlation sees
  t_film <- (D_h / 2) * (1 - sqrt(alpha_th))      # annular film thickness
  L_c    <- max(2 * t_film, 0.02 * D_h)           # characteristic length
  # atomization sees the power-law apparent viscosity at the film shear rate
  gdot_noz  <- U_rel / max(t_film, 1e-5)          # nozzle shear rate ~1e5-1e7
  mu_eff_hs <- mu_app(gdot_noz) * (1 + 2.5 * alpha_g)
  SMD_ab <- 0.48 * L_c * (sigma_eff / (rho_A * U_rel^2 * L_c))^0.4 *
              (1 + 1 / ALR)^0.4 +
            0.15 * L_c * sqrt(mu_eff_hs^2 / (sigma_eff * rho_lm * L_c)) *
              (1 + 1 / ALR)
  # (b) bubbles small enough to ride inside the films burst on the final
  #     letdown and subdivide the fragments (multiplicative refinement);
  #     coarse (post-hold) bubbles vent to the gas core instead
  D_b_th  <- D_b_h * r_th^(1/3)                   # bubble size at the throat
  g_in    <- 1 / (1 + D_b_th / t_film)            # fraction riding inside film
  F_flash <- (1 + alpha_th * (r_flash - 1))^(-1/3)
  SMD     <- SMD_ab * (1 - g_in * (1 - F_flash))
  # (c) template-solvent flash: the letdown to ambient superheats the
  #     low-boiling emulsion droplets riding in the fragments (T_feed vs
  #     T_bp at 1 atm); partial nucleation flash-shatters them further
  #     (liquid-to-vapor expansion O(200x), damped by the ~20 % fraction
  #     that nucleates on the nozzle timescale)
  S_flashn <- 1 / (1 + exp(-(T_feed - T_bpS) / 5))
  F_fs     <- (1 + 40 * phi_e * S_flashn)^(-1/3)
  SMD      <- SMD * F_fs

  ## --- Module 5a: Secondary breakup (Weber-number correction) --------------
  We_g <- rho_A * U_rel^2 * SMD / sigma_eff       # gas Weber number of drop
  f_sec <- 0.55 + 0.45 / (1 + We_g / 80)          # ~20-45 % diameter reduction
  d_drop <- SMD * f_sec                           # droplet scale entering dryer

  ## --- Module 5b: Droplet size distribution (three-mode volume mixture) ----
  # Main atomized mode plus two tail mechanisms observed with pre-gassed
  # feeds: (i) starved atomizing air leaves unbroken ligaments -> coarse
  # tail (raises d90/d99), amplified by viscous resistance (Ohnesorge);
  # (ii) excess atomizing air shear-strips satellites, and bursting bubble
  # films add fine debris -> fine mode (bimodality, pulls d10 down)
  E_a  <- ALR * U_rel^2                           # specific atomizing energy
  Oh_e <- mu_eff_hs / sqrt(rho_lm * sigma_eff * L_c) # viscous breakup resistance
  w_c  <- min(0.6, 0.5 / (1 + E_a / 1.5e5) * (1 + 3 * Oh_e))   # coarse tail
  w_f  <- min(0.35, 0.30 / (1 + 6e5 / E_a) +      # shear-strip satellites
                    0.15 * g_in * alpha_th +      # + bubble-film debris
                    0.50 * phi_e * S_flashn)      # + solvent-flash satellites
  if ((w_c + w_f) > 0.7) { sc <- 0.7 / (w_c + w_f); w_c <- w_c * sc; w_f <- w_f * sc }
  w_m  <- 1 - w_c - w_f                           # main mode weight

  m_main <- 1.2 * d_drop                          # volume median ~ 1.2 x D32
  modes_w <- c(w_m, w_c, w_f)
  modes_m <- c(m_main, 3.2 * m_main, 0.25 * m_main)
  modes_s <- c(1.9, 1.5, 1.6)                     # geometric std devs

  gr  <- seq(log(0.02 * m_main), log(20 * m_main), length.out = 400)
  cdf <- Reduce(`+`, Map(function(w, m, s)
           w * pnorm((gr - log(m)) / log(s)), modes_w, modes_m, modes_s))
  q_at <- function(p) exp(approx(cdf, gr, xout = p, ties = "ordered")$y)
  d10 <- q_at(0.10); d50 <- q_at(0.50); d90 <- q_at(0.90); d99 <- q_at(0.99)
  span <- (d90 - d10) / d50
  BI   <- 2 * min(w_f, 1 - w_f)                   # fine-mode bimodality index

  ## --- Module 6a: Dryer energy / moisture balance (co-current) -------------
  mdot_w <- mdot_L * (1 - w_e) * (1 - C_sol)      # water evaporation load
  mdot_e <- mdot_L * w_e                          # template solvent load
  T_out  <- T_in - (mdot_w * h_fg + mdot_e * h_fg_solv) / (mdot_g * cp_gas)
  T_out  <- max(T_out, T_feed + 2)                # saturated / underpowered dryer
  Y_out  <- Y_in + mdot_w / mdot_g                # outlet absolute humidity
  p_v    <- Y_out / (Y_out + 0.622) * P_atm       # vapour partial pressure
  p_sat  <- 610.94 * exp(17.625 * (T_out - 273.15) / (T_out - 273.15 + 243.04))
  RH_out <- min(p_v / p_sat, 0.99)                # outlet relative humidity

  ## --- Module 6b: Drying kinetics -------------------------------------------
  # d2-law evaporation coefficient, scaled by outlet-state driving force
  kappa <- 5e-8 * pmax(T_out - T_feed, 1) / 100 * (1 - RH_out)  # [m2/s]
  # Stokes-Einstein diffusivity of primary colloid particles (milling-affected)
  D_diff <- kB * T_feed / (6 * pi * mu_L * a_prim_mod)          # [m2/s]
  Pe <- kappa / (8 * D_diff)                      # drying Peclet number

  # theta_skin,z : skin network fraction. Colloidal instability (low DLVO
  # stability) accelerates shell aggregation; Flory-Huggins softness delays
  # rigid lock-in (soft particles film-form instead of jamming)
  S_skin <- Pe * C_sol * (1 + 1 / stability)
  S_crit <- 500 * (1 + 0.05 * Softness)
  theta_skin <- S_skin / (S_skin + S_crit)

  # Fractal openness (v47 Sect. 9): unstable colloid -> DLCA, open flocs
  # (low D_f, high structural porosity); stable -> compact RLCA packing
  D_f <- 1.8 + 0.7 * min((stability - 1) / 5, 1)
  phi_struct <- 0.30 * theta_skin * (1 - d_ratio^(3 - D_f))

  ## --- Module 6c: Per-mode drying & residual moisture ----------------------
  # d2-law drying time per droplet mode, retarded by the skin (falling-rate
  # period); the skin resistance scales inversely with the shell
  # permeability - plasticized/swollen shells stay open to vapor transport
  # while binder-blocked shells choke it; the coarse tail can leave the
  # tower wet
  tau_dry <- (modes_m^2 / kappa) * (1 + 20 * theta_skin / Perm_shell)
  X_j     <- exp(-t_res / tau_dry)                # residual water per mode
  X_moist <- sum(modes_w * X_j)                   # fraction of initial water

  ## --- Module 7a: Product glass transition (Fox: solvent + moisture) -------
  # in co-current drying particles approach the outlet gas temperature
  T_particle <- 0.85 * T_out + 0.15 * T_feed
  # template solvent boils only if the particle crosses its boiling point;
  # below it the immiscible solvent has no escape route (solubility-limited
  # diffusion through water is negligible) and stays liquid in the matrix
  S_bs   <- 1 / (1 + exp(-(T_particle - T_bpS) / 8))
  R_solv <- 1 - S_bs                              # retained solvent fraction

  # residual solvent (monomer partially evaporates; plasticizer stays;
  # unboiled template solvent adds to the load) plus water carried out by
  # the under-dried coarse tail - moisture plasticizes
  w_res  <- (0.5 * C_mono + C_plas + w_e * R_solv) /
            (0.5 * C_mono + C_plas + w_e * R_solv + C_sol)
  w_wat  <- min(0.12, X_moist * (1 - C_sol) /
                      (X_moist * (1 - C_sol) + C_sol))
  Tg_dry <- 1 / ((1 - w_res) / Tg_pol + w_res / Tg_solv)
  Tg_eff <- 1 / ((1 - w_wat) / Tg_dry + w_wat / Tg_water)

  # stickiness: particle temperature vs Tg_eff + 20 K sticky-point offset
  S_stick <- 1 / (1 + exp(-(T_particle - (Tg_eff + 20)) / 10))

  ## --- Module 7b: Cake (consolidated shell) mechanics -----------------------
  # Rumpf-type yield strength of the particulate cake vs meniscus capillary
  # stress decides collapse vs lock-in. Coagulated (low-stability) cakes
  # bond rigidly; Flory-Huggins softness and above-Tg mobility weaken them
  phi_cake <- 0.45 + 0.19 * min((stability - 1) / 5, 1)
  sigma_y  <- sig_y0 * (phi_cake / 0.64)^4 * (1 + 2 / stability) /
              sqrt(Softness) * (1 - 0.85 * S_stick)
  P_cap    <- 2 * sigma_eff / (0.3 * a_prim_mod)  # meniscus capillary stress
  Pi_col   <- P_cap / max(sigma_y, 1e3)           # collapse severity number
  f_col    <- Pi_col^2 / (1 + Pi_col^2)           # 0 = locks open, 1 = yields

  ## --- Module 7c: Per-mode porosity, sphericity, particle distribution -----
  # Template-solvent fate: vapor generated under the skin must permeate the
  # shell (same f_trap resistance as water); the Clausius-Clapeyron
  # overpressure of the superheated solvent plays against the cake yield
  # stress - a strong shell balloons (inflated pores), a weak one bursts
  # (micro-explosion -> broken hollow shells); a permeable shell vents
  # cleanly, leaving pores at the pre-processed template size
  f_trap_s <- theta_skin / (theta_skin + Perm_shell)
  dP_vap   <- P_atm * pmax(exp(0.025 * (T_particle - T_bpS)) - 1, 0)
  Pi_b     <- dP_vap * f_trap_s * S_bs / max(sigma_y, 1e3)
  f_burst  <- Pi_b^2 / (1 + Pi_b^2)
  B_infl   <- 1.5 * f_trap_s * (1 - f_burst)      # balloon inflation factor
  phi_templ <- min(0.6, phi_e * S_bs * (1 + B_infl) * (1 - 0.5 * f_burst))
  D_pore    <- D_e_h * (1 + B_infl)^(1/3)         # templated pore size

  # gas retention per droplet mode (coarse droplets keep more bubbles)
  a_trap_j <- alpha_e / (1 + D_b_e / modes_m) * (0.5 + 0.5 * theta_surf)
  phi_str_eff <- phi_struct * (1 - 0.7 * f_col)   # yielding crushes pores
  # vapor entrapment: water evaporating beneath a formed skin must permeate
  # the shell; low-permeability shells on superheated droplets build
  # internal vapor pressure and inflate vacuoles (hollow particles)
  f_trap  <- theta_skin / (theta_skin + Perm_shell)
  S_boil  <- 1 / (1 + exp(-(T_particle - 373.15) / 8))  # superheat propensity
  phi_vac <- 0.30 * f_trap * S_boil
  phi_int <- 1 - (1 - phi_str_eff) * (1 - phi_vac) *    # structural + vacuole
                 (1 - phi_templ)                        # + templated pores
  phi_j <- a_trap_j + (1 - a_trap_j) * phi_int
  phi_j <- phi_j * (1 - 0.6 * S_stick * (1 - 0.5 * theta_skin))
  phi_porosity <- sum(modes_w * phi_j)            # mass-weighted mixture

  # Omega_struct,z : pre-shell rounding resisted by sigma and low-shear
  # viscosity; post-shell buckling requires the cake to yield (f_col)
  resist <- (sigma_eff / 0.05) * (mu_eff0 / 0.01)^0.3
  Omega_struct <- 1 - 0.55 * theta_skin * (0.3 + 0.7 * f_col) /
                    ((1 + resist) * (1 + 0.1 * Softness))
  Omega_struct <- Omega_struct + (1 - Omega_struct) * 0.4 * S_stick
  # over-pressurized impermeable shells rupture (blowholes and template
  # micro-explosions), denting the sphericity inflation would restore
  Omega_struct <- Omega_struct * (1 - 0.30 * phi_vac) *
                  (1 - 0.35 * f_burst * min(3 * phi_e, 1))

  # per-mode particle diameters (solids balance) and mixture quantiles
  Dp_j <- modes_m * ((1 - a_trap_j) * (1 - phi_e) * rho_L * C_sol /
                     (rho_s * (1 - phi_j)))^(1/3)
  # Colloid milling coupling: the physical primary size (a_prim_mod from the UP1
  # three-regime closure) sets packing efficiency in the dried aggregate.
  # Larger primaries (less milled) pack more loosely -> larger final particle.
  # Reference: a_prim_ref_min = minimum radius at full milling = a_prim × (1-MILL_MAX)
  #            = 100 nm × 0.30 = 30 nm. At this minimum, factor = 1.0 (baseline size).
  # Exponent PSF_EXP calibrated: at 200 nm primary (v_tip=8.5 m/s, no milling),
  # full chain (UP1->UP2 foam-wash->UP3->UP4) at 20% solid with C_AGG_CAL=0.0797,
  # target UP4 d50 = 15.3 µm → PSF_EXP = 0.537 (factor = 3.333^0.537 = 1.926).
  # (Prior value 0.588 was calibrated at 35% nominal solid, bypassing the foam-wash
  #  column; 0.537 is the full-chain anchor at the 20% solid reference point.)
  a_prim_ref_min <- a_prim * 0.30   # 3e-8 m: minimum primary radius (full milling)
  PSF_EXP        <- 0.537           # calibrated: UP4 d50 = 15.3 µm at 8.5 m/s / 20% solid
  packing_size_factor <- (a_prim_mod / a_prim_ref_min)^PSF_EXP
  Dp_j <- Dp_j * packing_size_factor

  grp  <- seq(log(0.05 * min(Dp_j)), log(10 * max(Dp_j)), length.out = 400)
  cdfp <- Reduce(`+`, Map(function(w, m, s)
            w * pnorm((grp - log(m)) / log(s)), modes_w, Dp_j, modes_s))
  qp_at <- function(p) exp(approx(cdfp, grp, xout = p, ties = "ordered")$y)
  Dp50 <- qp_at(0.50); Dp90 <- qp_at(0.90)

  ## --- Module 7d: Tapped density (v47 Sect. 10 density closure) ------------
  # bimodal blends pack better (fines fill interstices); sticky, plasticized
  # or moist powders cake and pack worse
  rho_env <- sum(modes_w * ((1 - phi_j) * rho_s + phi_j * 1.2))
  f_pack  <- 0.64 * (0.55 + 0.45 * Omega_struct) *
             (Dp50 / (Dp50 + 8e-6)) *
             (1 + 0.5 * modes_w[3] * modes_w[1]) *
             (1 - 0.2 * S_stick) *
             (1 - 0.15 * min(Softness / 25, 1)) *
             (1 - 2 * min(X_moist, 0.2))
  rho_tapped <- rho_env * max(f_pack, 0.05)

  c(d_droplet_um   = d50 * 1e6,
    d10_um         = d10 * 1e6,
    d90_um         = d90 * 1e6,
    d99_um         = d99 * 1e6,
    span           = span,
    BI_bimodal     = BI,
    D_particle_um  = Dp50 * 1e6,
    Dp90_um        = Dp90 * 1e6,
    theta_skin_z   = theta_skin,
    Omega_struct_z = Omega_struct,
    phi_porosity_z = phi_porosity,
    rho_tapped     = rho_tapped,
    Tg_product_K   = Tg_eff,
    X_moisture     = X_moist,
    Perm_shell_rel = Perm_shell,
    D_pore_um      = D_pore * 1e6,
    solv_retained  = R_solv,
    f_burst_solv   = f_burst,
    sigma_y_cake_MPa = sigma_y / 1e6)
}

run_model <- function(X01) {
  Xphys <- scale_design(as.matrix(X01))
  t(apply(Xphys, 1, spray_dry_model))
}

# -----------------------------------------------------------------------------
# 3. Morris design: `sensitivity` package if available, else base-R fallback
# -----------------------------------------------------------------------------
r_traj <- 60   # number of Morris trajectories
levels <- 8    # grid levels
delta  <- levels / (2 * (levels - 1))   # standard Morris step

morris_oat_design <- function(k, r, levels, delta) {
  grid <- seq(0, 1 - delta, length.out = levels - levels / 2)  # feasible bases
  X <- matrix(NA_real_, nrow = r * (k + 1), ncol = k)
  info <- vector("list", r)
  row <- 1
  for (t in seq_len(r)) {
    base <- sample(grid, k, replace = TRUE)
    ord  <- sample.int(k)
    dirs <- numeric(k)
    x <- base
    X[row, ] <- x
    for (s in seq_len(k)) {
      j <- ord[s]
      d <- sample(c(-delta, delta), 1)
      if (x[j] + d < 0 || x[j] + d > 1) d <- -d
      x[j] <- x[j] + d
      dirs[s] <- d
      X[row + s, ] <- x
    }
    info[[t]] <- list(order = ord, dirs = dirs)
    row <- row + k + 1
  }
  list(X = X, info = info)
}

elementary_effects <- function(design, Y, k, r) {
  # Y: matrix, one column per output; returns list of (r x k) EE matrices
  n_out <- ncol(Y)
  ee <- lapply(seq_len(n_out), function(i) matrix(NA_real_, r, k))
  for (t in seq_len(r)) {
    off <- (t - 1) * (k + 1)
    ord <- design$info[[t]]$order
    dirs <- design$info[[t]]$dirs
    for (s in seq_len(k)) {
      j <- ord[s]
      dy <- Y[off + s + 1, ] - Y[off + s, ]
      for (i in seq_len(n_out)) ee[[i]][t, j] <- dy[i] / dirs[s]
    }
  }
  ee
}

use_pkg <- requireNamespace("sensitivity", quietly = TRUE)

if (use_pkg) {
  message("Using sensitivity::morris()")
  mor <- sensitivity::morris(model = NULL, factors = factors$name, r = r_traj,
                             design = list(type = "oat", levels = levels,
                                           grid.jump = levels / 2),
                             binf = 0, bsup = 1)
  Y <- run_model(mor$X)
  ee_list <- lapply(seq_len(ncol(Y)), function(i) {
    m <- mor
    sensitivity::tell(m, Y[, i])
    m$ee
  })
} else {
  message("Package 'sensitivity' not found - using built-in Morris OAT design")
  design <- morris_oat_design(k, r_traj, levels, delta)
  Y <- run_model(design$X)
  ee_list <- elementary_effects(design, Y, k, r_traj)
}

outputs <- colnames(Y)

morris_stats <- function(ee) {
  data.frame(factor  = factors$name,
             mu      = colMeans(ee),
             mu.star = colMeans(abs(ee)),
             sigma   = apply(ee, 2, sd))
}
stats_list <- setNames(lapply(ee_list, morris_stats), outputs)

# -----------------------------------------------------------------------------
# 4. Morris plots (mu* vs sigma), one panel per output
# -----------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)

titles <- c(d_droplet_um   = "Spray droplet size  Dv50 [um]",
            d10_um         = "Fines tail  d10 [um]",
            d90_um         = "Coarse tail  d90 [um]",
            d99_um         = "Extreme coarse tail  d99 [um]",
            span           = "Span  (d90-d10)/d50 [-]",
            BI_bimodal     = "Bimodality index [-]",
            D_particle_um  = "Final particle size  Dp50 [um]",
            Dp90_um        = "Particle coarse tail  Dp90 [um]",
            theta_skin_z   = "Skin formation  theta_skin,z [-]",
            Omega_struct_z = "Sphericity  Omega_struct,z [-]",
            phi_porosity_z = "Porosity  phi_porosity,z [-]",
            rho_tapped     = "Tapped density  rho_tapped [kg/m3]",
            Tg_product_K   = "Product glass transition  Tg_eff [K]",
            X_moisture     = "Residual moisture (frac. of feed water) [-]",
            Perm_shell_rel = "Shell permeability (relative)  Perm_shell [-]",
            D_pore_um      = "Templated pore size  D_pore [um]",
            solv_retained  = "Retained template solvent (frac. of load) [-]",
            f_burst_solv   = "Template micro-explosion severity  f_burst [-]",
            sigma_y_cake_MPa = "Cake yield strength  sigma_y [MPa]")

plot_morris_panel <- function(st, title, n_label = 10) {
  plot(st$mu.star, st$sigma,
       xlim = c(0, max(st$mu.star) * 1.25),
       ylim = c(0, max(st$sigma) * 1.30),
       pch = 21, bg = "steelblue", cex = 1.3,
       xlab = expression(mu * "*  (mean |elementary effect|)"),
       ylab = expression(sigma * "  (sd of elementary effects)"),
       main = title, cex.main = 0.95)
  abline(0, 1, lty = 2, col = "grey55")           # sigma = mu* (non-linear)
  abline(0, 0.1, lty = 3, col = "grey70")         # sigma = 0.1 mu* (~linear)
  # label only the strongest factors to keep panels readable
  top <- order(-st$mu.star)[seq_len(min(n_label, nrow(st)))]
  text(st$mu.star[top], st$sigma[top], labels = st$factor[top], pos = 3,
       cex = 0.70, offset = 0.35, xpd = NA)
}

png(file.path("output", "morris_sensitivity_plots.png"),
    width = 3000, height = 3750, res = 160)
op <- par(mfrow = c(5, 4), mar = c(4.5, 4.5, 2.5, 1), oma = c(2.5, 0, 0, 0))
for (o in outputs) plot_morris_panel(stats_list[[o]], titles[[o]])
mtext(sprintf(paste("Morris screening: %d trajectories, %d levels, %d model",
                    "runs | dashed: sigma = mu* (non-linear/interacting),",
                    "dotted: sigma = 0.1 mu* (near-linear) |",
                    "top 10 factors labelled per panel"),
              r_traj, levels, nrow(Y)),
      side = 1, outer = TRUE, cex = 0.75, line = 1)
par(op); dev.off()

# -----------------------------------------------------------------------------
# 5. Console summary + CSV export
# -----------------------------------------------------------------------------
all_stats <- do.call(rbind, lapply(outputs, function(o)
  cbind(output = o, stats_list[[o]])))
write.csv(all_stats, file.path("output", "morris_indices.csv"),
          row.names = FALSE)

for (o in outputs) {
  st <- stats_list[[o]][order(-stats_list[[o]]$mu.star), ]
  cat("\n==", titles[[o]], "==\n")
  print(st, row.names = FALSE, digits = 3)
}
cat("\nWrote output/morris_sensitivity_plots.png and output/morris_indices.csv\n")
