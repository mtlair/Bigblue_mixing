# UP1->UP4-direct Morris screen with the overlays ON (size_template=1,
# morphology_recal=1). Ranks which factors drive the reconciled porous/hollow
# product properties. Writes unified_output/direct_morris_overlays.csv.
lines <- readLines("unified_model.R")
cut <- grep("^cache_file <-", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))

run_ov <- function(x) { x["size_template"] <- 1; x["morphology_recal"] <- 1; run_unified(x) }

set.seed(1)
design <- morris_oat_design(k, r_traj, levels, delta)
Xphys  <- scale_design(design$X)
rows   <- lapply(seq_len(nrow(Xphys)), function(i) Xphys[i, ])
cores  <- max(1, parallel::detectCores() - 1)
ylist  <- parallel::mclapply(rows, run_ov, mc.cores = cores)
Y <- do.call(rbind, ylist)
cat(sprintf("Morris (overlays on): %d factors, %d runs, %d valid\n",
            k, nrow(Xphys), sum(stats::complete.cases(Y))))
ee <- elementary_effects(design, Y, k, r_traj)
names(ee) <- unified_output_names

targets <- c("up2_D_particle_um","up2_phi_porosity_z","up2_Omega_struct_z",
             "up2_theta_skin_z","up2_rho_tapped")
tnice <- c("size Dp50","porosity","sphericity","surface_fusion","tapped density")
allrows <- list()
for (oi in seq_along(targets)) {
  o <- targets[oi]; m <- ee[[o]]
  mustar <- colMeans(abs(m), na.rm = TRUE)
  sig    <- apply(m, 2, sd, na.rm = TRUE)
  ord <- order(-mustar)
  df <- data.frame(output = tnice[oi], factor = factors$name[ord],
                   module = factors$module[ord],
                   mu_star = round(mustar[ord], 4), sigma = round(sig[ord], 4))
  allrows[[o]] <- df
  cat(sprintf("\n=== %s  (top 6 drivers by mu*) ===\n", tnice[oi]))
  for (j in 1:6) cat(sprintf("  %-16s [%-4s]  mu*=%.3g\n",
                             df$factor[j], df$module[j], df$mu_star[j]))
}
out <- do.call(rbind, allrows)
write.csv(out, "unified_output/direct_morris_overlays.csv", row.names = FALSE)
cat("\nWrote unified_output/direct_morris_overlays.csv\n")
