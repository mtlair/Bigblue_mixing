# =============================================================================
# Focused Morris screen — FORMULATION / surface-chemistry variables
# =============================================================================
# Screens the monomer / plasticizer / binder / solvent design space on the wired
# chain (UP1->UP2->UP3->UP4, size-template on), holding the process setpoints at
# nominal, and ranks their effect on the morphology / structure CQAs. Concentrations
# are dosable formulation knobs; the k_perm_* coefficients are each additive's
# intrinsic shell-permeation property (informational -- tells you which additive
# chemistry the response is most sensitive to).
#
# Run:  Rscript doe_formulation_morris.R
# =============================================================================
lines <- readLines("unified_model.R")
cut   <- grep("^# 5\\. MORRIS SCREEN", lines)[1] - 1
eval(parse(text = paste(lines[1:cut], collapse = "\n")))
source("foam_wash_module.R"); source("up3_centrifuge_module.R")
options(unified.wire_up2 = TRUE, unified.wire_up3 = TRUE)
out_dir <- "unified_output"; dir.create(out_dir, showWarnings = FALSE)

# --- formulation subset (name, min, max, log, dosable?) ----------------------
sub <- rbind(
  data.frame(name="C_monomer",    min=0.001, max=0.03, log=FALSE, kind="dose:monomer"),
  data.frame(name="C_plasticizer",min=0.005, max=0.04, log=FALSE, kind="dose:plasticizer"),
  data.frame(name="C_binder",     min=0.01,  max=0.15, log=FALSE, kind="dose:binder"),
  data.frame(name="template_dose",min=0.02,  max=1.10, log=FALSE, kind="dose:solvent(amount)"),
  data.frame(name="T_bp_solv",    min=300,   max=360,  log=FALSE, kind="solvent(bp/type)"),
  data.frame(name="k_perm_mono",  min=5,     max=60,   log=FALSE, kind="material:monomer perm"),
  data.frame(name="k_perm_plast", min=5,     max=60,   log=FALSE, kind="material:plast perm"),
  data.frame(name="k_perm_bind",  min=2,     max=40,   log=FALSE, kind="material:binder perm"),
  stringsAsFactors = FALSE)
ks <- nrow(sub); r_traj <- 40; levels <- 8
delta <- levels / (2 * (levels - 1))

scale_sub <- function(u) sub$min + u * (sub$max - sub$min)   # all linear here

run_form <- function(u) {
  x <- nominal_x
  x[sub$name] <- scale_sub(u)
  x["size_template"] <- 1
  run_unified(x, template_type = TEMPLATE_TYPE)
}

# OAT design over the subset
set.seed(7)
grid <- seq(0, 1 - delta, length.out = levels - levels/2)
X <- matrix(NA_real_, r_traj*(ks+1), ks); info <- vector("list", r_traj); row <- 1
for (t in seq_len(r_traj)) {
  base <- sample(grid, ks, TRUE); ord <- sample.int(ks); x <- base; X[row,] <- x
  for (s in seq_len(ks)) { j<-ord[s]; d<-sample(c(-delta,delta),1)
    if (x[j]+d<0||x[j]+d>1) d<--d; x[j]<-x[j]+d; X[row+s,]<-x }
  dirs <- numeric(ks); for (s in seq_len(ks)){j<-ord[s]; dirs[s]<-X[row+s,j]-X[row+s-1,j]}
  info[[t]] <- list(order=ord, dirs=dirs); row <- row+ks+1
}
cat(sprintf("Formulation Morris: %d factors, r=%d -> %d wired runs...\n", ks, r_traj, nrow(X)))
Y <- do.call(rbind, parallel::mclapply(seq_len(nrow(X)),
       function(i) run_form(X[i,]), mc.cores = max(1, parallel::detectCores()-1)))
cat(sprintf("  %d/%d valid\n", sum(stats::complete.cases(Y)), nrow(Y)))
colnames(Y) <- unified_output_names

ee <- lapply(seq_len(ncol(Y)), function(o) {
  m <- matrix(NA_real_, r_traj, ks)
  for (t in seq_len(r_traj)) { off<-(t-1)*(ks+1); ord<-info[[t]]$order; dr<-info[[t]]$dirs
    for (s in seq_len(ks)){ j<-ord[s]; m[t,j] <- (Y[off+s+1,o]-Y[off+s,o])/dr[s] } }
  m })
names(ee) <- unified_output_names
stat <- function(m) data.frame(factor=sub$name, kind=sub$kind,
  mu=colMeans(m,na.rm=TRUE), mu.star=colMeans(abs(m),na.rm=TRUE), sigma=apply(m,2,sd,na.rm=TRUE))
sl <- setNames(lapply(ee, stat), unified_output_names)
write.csv(do.call(rbind, lapply(unified_output_names, function(o) cbind(output=o, sl[[o]]))),
          file.path(out_dir, "doe_formulation_morris_indices.csv"), row.names=FALSE)

CQA <- c(up2_Tg_product_K="product Tg [K]", up2_theta_skin_z="skin/fusion",
  up2_phi_porosity_z="porosity", up2_Omega_struct_z="sphericity",
  up2_Perm_shell_rel="shell permeability", up2_D_pore_um="pore size [um]",
  up2_solv_retained="retained solvent", up2_f_burst_solv="micro-explosion risk",
  up2_sigma_y_cake_MPa="cake strength [MPa]", up2_D_particle_um="d50 [um]")

sink(file.path(out_dir, "doe_formulation_proposal.txt"))
cat("=====================================================================\n")
cat("FORMULATION DoE - monomer/plasticizer/binder/solvent (wired chain)\n")
cat("Morris r=40; process setpoints held at nominal. mu*=influence,\n")
cat("sign(mu)=direction. Dosable knobs get a proposed range; k_perm_* are\n")
cat("intrinsic material properties (which additive chemistry matters).\n")
cat("=====================================================================\n\n")
for (o in names(CQA)) {
  st <- sl[[o]][order(-sl[[o]]$mu.star),]
  cat(sprintf("### %-20s %s\n", o, CQA[[o]]))
  for (i in 1:4) cat(sprintf("   %-14s %-20s mu*=%.4g  %s\n",
      st$factor[i], st$kind[i], st$mu.star[i], ifelse(st$mu[i]>=0,"increase","decrease")))
  cat("\n")
}
cat("---------------------------------------------------------------------\n")
cat("DOSABLE FORMULATION KNOBS - proposed DoE ranges:\n")
cat("---------------------------------------------------------------------\n")
for (i in which(grepl("^dose|solvent\\(", sub$kind))) {
  fn <- sub$name[i]
  # rank across CQAs: count top-3 appearances
  hits <- sum(sapply(names(CQA), function(o){
    st<-sl[[o]][order(-sl[[o]]$mu.star),]; fn %in% st$factor[1:3] }))
  cat(sprintf("  %-14s %-20s [%.4g, %.4g]  nominal %.4g   (top-3 for %d/%d CQAs)\n",
      fn, sub$kind[i], sub$min[i], sub$max[i], nominal_x[[fn]], hits, length(CQA)))
}
sink()
cat("Wrote doe_formulation_morris_indices.csv and doe_formulation_proposal.txt\n")
writeLines(readLines(file.path(out_dir, "doe_formulation_proposal.txt")))
