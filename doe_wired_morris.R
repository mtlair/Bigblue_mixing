# =============================================================================
# DoE screen on the WIRED full chain  UP1 -> UP2 -> UP3 -> UP4
# =============================================================================
# Runs a Morris elementary-effects screen on the fully wired chain (foam wash +
# separator enabled, size-template morphology path ON) to rank which factors
# actually move each product CQA, then proposes experimental bounds for a DoE.
#
# Differs from unified_model.R's screen (which runs UP1->UP4 DIRECT with the two
# intermediate stages muted): here intermediate_stage_1/2 are the real UP2/UP3
# models, so the dryer sees the concentrated ~40% cake, not the raw 25% slurry.
#
# Calibration anchors from SEM (review/sem_morphology_metrics.md):
#   phi_porosity_z  <-> surface open-porosity ~37% (ROBUST)
#   Omega_struct_z  <-> sphericity ~0.91 (ROBUST)
#   theta_skin_z    <-> fusion (DEFERRED: SEM surface_fusion not repeatable per-sample yet)
#
# Run:  Rscript doe_wired_morris.R
# =============================================================================

lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))
source("foam_wash_module.R")
source("up3_separator_module.R")
options(unified.wire_up2 = TRUE, unified.wire_up3 = TRUE)

out_dir <- "unified_output"; dir.create(out_dir, showWarnings = FALSE)
r_traj <- 30; levels <- 8
delta  <- levels / (2 * (levels - 1))

# wired chain run with the size-template morphology path enabled
run_doe <- function(x) {
  x2 <- x; x2["size_template"] <- 1
  run_unified(x2, template_type = TEMPLATE_TYPE)
}

# --- Morris OAT design + elementary effects (same as unified_model.R) ---------
morris_oat_design <- function(k, r, levels, delta) {
  grid <- seq(0, 1 - delta, length.out = levels - levels / 2)
  X <- matrix(NA_real_, nrow = r * (k + 1), ncol = k); info <- vector("list", r); row <- 1
  for (t in seq_len(r)) {
    base <- sample(grid, k, replace = TRUE); ord <- sample.int(k); x <- base
    X[row, ] <- x
    for (s in seq_len(k)) {
      j <- ord[s]; d <- sample(c(-delta, delta), 1)
      if (x[j] + d < 0 || x[j] + d > 1) d <- -d
      x[j] <- x[j] + d; X[row + s, ] <- x
    }
    info[[t]] <- list(order = ord, dirs = rep(delta, k) * sign(diff(c(base[ord[1]], X[row + 1, ord[1]]))))
    # recompute dirs robustly
    dirs <- numeric(k); xb <- base
    for (s in seq_len(k)) { j <- ord[s]; dirs[s] <- X[row + s, j] - X[row + s - 1, j] }
    info[[t]] <- list(order = ord, dirs = dirs)
    row <- row + k + 1
  }
  list(X = X, info = info)
}
elementary_effects <- function(design, Y, k, r) {
  no <- ncol(Y); ee <- lapply(seq_len(no), function(i) matrix(NA_real_, r, k))
  for (t in seq_len(r)) {
    off <- (t - 1) * (k + 1); ord <- design$info[[t]]$order; dirs <- design$info[[t]]$dirs
    for (s in seq_len(k)) {
      j <- ord[s]; dy <- Y[off + s + 1, ] - Y[off + s, ]
      for (i in seq_len(no)) ee[[i]][t, j] <- dy[i] / dirs[s]
    }
  }
  ee
}

set.seed(42)
design <- morris_oat_design(k, r_traj, levels, delta)
Xphys  <- scale_design(design$X)
n_runs <- nrow(Xphys)
cat(sprintf("Wired-chain Morris: %d factors, r=%d -> %d runs (size_template=1)...\n",
            k, r_traj, n_runs))
rows <- lapply(seq_len(n_runs), function(i) Xphys[i, ])
ylist <- if (.Platform$OS.type == "unix")
  parallel::mclapply(rows, run_doe, mc.cores = max(1, parallel::detectCores() - 1)) else
  lapply(rows, run_doe)
Y <- do.call(rbind, ylist)
cat(sprintf("  %d/%d runs valid\n", sum(stats::complete.cases(Y)), n_runs))
saveRDS(list(design = design, Y = Y), file.path(out_dir, "doe_wired_morris.rds"))

ee_list <- setNames(elementary_effects(design, Y, k, r_traj), unified_output_names)
morris_stats <- function(ee) data.frame(
  factor = factors$name, module = factors$module, group = factors$group,
  mu = colMeans(ee, na.rm = TRUE), mu.star = colMeans(abs(ee), na.rm = TRUE),
  sigma = apply(ee, 2, sd, na.rm = TRUE))
stats_list <- setNames(lapply(ee_list, morris_stats), unified_output_names)
all_stats <- do.call(rbind, lapply(unified_output_names, function(o)
  cbind(output = o, stats_list[[o]])))
write.csv(all_stats, file.path(out_dir, "doe_wired_morris_indices.csv"), row.names = FALSE)

# --- CQAs to design against + robust SEM anchor -------------------------------
CQA <- c(up2_D_particle_um = "particle size d50 [um]",
         up2_phi_porosity_z = "porosity  (SEM anchor ~0.37, ROBUST)",
         up2_Omega_struct_z = "sphericity(SEM anchor ~0.91, ROBUST)",
         up2_rho_tapped     = "tapped density [kg/m3]",
         up2_Tg_product_K   = "product Tg [K]",
         up2_X_moisture     = "residual moisture [-]",
         up2_theta_skin_z   = "surface_fusion/fusion (SEM DEFERRED - not repeatable)")

sink(file.path(out_dir, "doe_proposal.txt"))
cat("========================================================================\n")
cat("DoE PROPOSAL - wired chain UP1->UP2->UP3->UP4  (Morris screen, r=30)\n")
cat("Top factors per CQA (mu* = influence, sign of mu = direction), with the\n")
cat("model factor range as the proposed experimental bound.\n")
cat("========================================================================\n\n")
topfac <- list()
for (o in names(CQA)) {
  st <- stats_list[[o]]; st <- st[order(-st$mu.star), ]
  cat(sprintf("### %-18s : %s\n", o, CQA[[o]]))
  cat(sprintf("  %-16s %-6s %10s %10s  %-9s %s\n",
              "factor", "module", "mu*", "sigma", "direction", "range [min,max] unit"))
  for (i in seq_len(6)) {
    fn <- st$factor[i]; fr <- factors[factors$name == fn, ]
    dir <- if (st$mu[i] >= 0) "increase" else "decrease"
    cat(sprintf("  %-16s %-6s %10.4g %10.4g  %-9s [%.3g, %.3g] %s\n",
                fn, st$module[i], st$mu.star[i], st$sigma[i], dir, fr$min, fr$max, fr$unit))
  }
  topfac[[o]] <- head(st$factor, 4)
  cat("\n")
}

# union of top-4 drivers across CQAs -> recommended DoE factor set
uni <- sort(table(unlist(topfac)), decreasing = TRUE)
cat("------------------------------------------------------------------------\n")
cat("RECOMMENDED DoE FACTOR SET (top-4 driver of >=2 CQAs), with bounds:\n")
cat("------------------------------------------------------------------------\n")
cat(sprintf("  %-16s %-6s %-8s %-22s %s\n", "factor", "module", "#CQAs", "range [min,max] unit", "nominal"))
sel <- names(uni)[uni >= 2]
for (fn in sel) {
  fr <- factors[factors$name == fn, ]
  cat(sprintf("  %-16s %-6s %-8d [%.3g, %.3g] %-10s %.3g\n",
              fn, fr$module, uni[[fn]], fr$min, fr$max, fr$unit, nominal_x[[fn]]))
}
sink()
cat("Wrote", file.path(out_dir, "doe_wired_morris_indices.csv"), "and doe_proposal.txt\n")
writeLines(readLines(file.path(out_dir, "doe_proposal.txt")))
