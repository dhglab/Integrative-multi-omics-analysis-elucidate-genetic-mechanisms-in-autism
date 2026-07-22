#!/usr/bin/env Rscript
# ==============================================================================
# run_demo.R
#
# This is the "Reproduction instructions" counterpart to demo/run_demo.R.
# demo/run_demo.R (simulated data) satisfies the Nature Portfolio checklist's
# required "Demo" section on its own -- this script is for reproducing the
# actual manuscript results from the real processed data, which is optional
# per the checklist but useful for readers with approved data access.
#
# Adapted from:
#   scripts/01_chromatin_accessibility/
#     5_4_v7_PCA_build_linear_model_SignificantCorWithPC_VIF2p5_lmeSubject_organized.R
#
# Expected inputs (repo-relative paths, assuming script is run from repo root):
#   results/5_DiffATAC/RPKM_cqn.rda        - loads object `RPKM.cqn`
#   results/3_QC/BiolTechCov_and_QC.rda    - loads object `df_wAllCov`
#
# Expected outputs (written to results/5_DiffATAC/):
#   Corrplot_topPCofCQNnormedRPKM_Covariates.pdf
#   p_magnitude_lme<model>.rda
#   HistPval_lme<model>.pdf
#
# Expected runtime: under 20 minutes on a normal desktop computer.
# ==============================================================================
 
rm(list = ls())
options(stringsAsFactors = FALSE)
library(tidyverse)
library(nlme)
library(corrplot)
 
repo_root <- "."  # run from repo root, or set explicitly
results_dir <- file.path(repo_root, "results", "5_DiffATAC")
qc_dir      <- file.path(repo_root, "results", "3_QC")
 
## ---- Load data --------------------------------------------------------------
load(file.path(results_dir, "RPKM_cqn.rda"))              # expects: RPKM.cqn
load(file.path(qc_dir, "BiolTechCov_and_QC.rda"))         # expects: df_wAllCov
 
Covariates <- df_wAllCov
Covariates <- Covariates[match(colnames(RPKM.cqn), Covariates$sample_name), ]
 
idx_rm <- which(grepl("_quality", colnames(Covariates)))
if (length(idx_rm) > 0) Covariates <- Covariates[, -idx_rm]
idx_rm <- which(colnames(Covariates) %in%
                  c("QC_summary", "Subject_ID", "Subject_ID_2", "sample_name"))
if (length(idx_rm) > 0) Covariates <- Covariates[, -idx_rm]
 
Covariates$PMI <- as.numeric(Covariates$PMI)
Covariates$PMI[is.na(Covariates$PMI)] <- mean(Covariates$PMI, na.rm = TRUE)
Covariates$Diagnosis <- factor(Covariates$Diagnosis, levels = c("CTL", "ASD"))
Covariates$Sex <- factor(Covariates$Sex, levels = c("M", "F"))
Covariates$Cortex <- factor(Covariates$Cortex, levels = c("Parietal", "Frontal", "Temporal"))
Covariates$Age2 <- Covariates$Age^2
 
colnames(Covariates) <- gsub("_", ".", colnames(Covariates))
colnames(Covariates)[which(colnames(Covariates) == "subject")] <- "Subject"
 
Covariates <- Covariates[, c("Subject", "Diagnosis", "Age", "Age2", "PMI", "Sex",
                              "BrainBank", "Batch", "Cortex", "Region",
                              "tssenrich.score", "FRiP",
                              "picard.alignment.PCT.READS.ALIGNED.IN.PAIRS",
                              "picard.duplication.PERCENT.DUPLICATION")]
colnames(Covariates)[(ncol(Covariates) - 1):ncol(Covariates)] <-
  c("PCT.READS.ALIGNED.IN.PAIRS", "PERCENT.DUPLICATION")
 
## ---- Remove technical replicates for PCA/covariate screening ---------------
TechRep_rm <- c("B4334-1", "B4721A", "B5000B", "B5718D", "B5813B", "CQ56-2")
RPKM.cqn2 <- RPKM.cqn[, which(!colnames(RPKM.cqn) %in% TechRep_rm)]
Covariates2 <- Covariates[which(!rownames(Covariates) %in% TechRep_rm), ]
 
stopifnot(rownames(Covariates) == colnames(RPKM.cqn))
stopifnot(rownames(Covariates2) == colnames(RPKM.cqn2))
 
## ---- Top PCs -----------------------------------------------------------------
norm <- t(scale(t(RPKM.cqn2), scale = FALSE))
pc <- prcomp(norm, center = FALSE)
varexp <- (pc$sdev)^2 / sum(pc$sdev^2)
n_pcs_80 <- which(cumsum(varexp) >= 0.80)[1]
top_pc <- pc$rotation[, 1:n_pcs_80]
cat(sprintf("Top %d PCs explain %.1f%% of total variance.\n",
            n_pcs_80, 100 * sum(varexp[1:n_pcs_80])))
 
## ---- PC-covariate correlation screening --------------------------------------
mod_mat_expr <- paste(colnames(Covariates2)[-1], collapse = " + ")
mod_mat_expr <- paste0("~ ", mod_mat_expr)
mod_mat <- model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates2)[, -1]
mod_mat_withPC <- cbind(top_pc, mod_mat)
 
cor_mat <- cor(mod_mat_withPC)
cor_spearman <- cor(mod_mat_withPC, method = "spearman")
idx_spearman <- ncol(top_pc) + c(1, 5:15)  # TODO: confirm these indices still
                                            # correspond to the intended
                                            # non-parametric covariates for
                                            # your current covariate set
cor_mat[idx_spearman, ] <- cor_spearman[idx_spearman, ]
cor_mat[, idx_spearman] <- cor_spearman[, idx_spearman]
 
cor_sig <- matrix(nrow = nrow(cor_mat), ncol = ncol(cor_mat),
                   dimnames = list(colnames(cor_mat), colnames(cor_mat)))
for (i in seq_len(ncol(mod_mat_withPC))) {
  for (j in seq_len(ncol(mod_mat_withPC))) {
    cor_sig[i, j] <- suppressWarnings(
      cor.test(mod_mat_withPC[, i], mod_mat_withPC[, j])$p.value
    )
  }
}
 
cor_fdr <- matrix(p.adjust(cor_sig, method = "fdr"), nrow = nrow(cor_sig),
                   dimnames = dimnames(cor_sig))
diag(cor_fdr) <- 1
 
flagged <- colnames(cor_fdr)[
  apply(cor_fdr[1:ncol(top_pc), , drop = FALSE], 2, function(x) any(x < 0.2))
]
cat("Covariates flagged at FDR < 0.2:\n  ", paste(flagged, collapse = ", "), "\n")
 
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(results_dir, "Corrplot_topPCofCQNnormedRPKM_Covariates.pdf"),
    height = 35, width = 45)
corrplot(cor_mat, method = "ellipse", tl.pos = "lt", tl.col = "black",
         tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
dev.off()
 
## ---- Final peak-wise linear mixed model --------------------------------------
expression_model <- "lme(logRPKM ~ Diagnosis + Age + Batch + RegionBA38 + RegionBA44_45 + tssenrich.score + FRiP, random=~1|Subject, data = cur_data)"
 
Covariates$Subject <- factor(Covariates$Subject)
Covariates$RegionBA38 <- ifelse(Covariates$Region == "BA38", 1, 0)
Covariates$RegionBA44_45 <- ifelse(Covariates$Region == "BA44_45", 1, 0)
 
n_terms <- 7  # number of fixed-effect coefficients (excl. intercept) -- update
              # if you change expression_model above
p_mat <- magnitude_mat <- matrix(nrow = nrow(RPKM.cqn), ncol = n_terms)
 
runlme <- function(thisdat, expression) {
  lm1 <- eval(parse(text = expression))
  s <- summary(lm1)
  list(s$coefficients$fixed, anova(lm1))
}
 
for (i in seq_len(nrow(RPKM.cqn))) {
  if (i %% 5000 == 0) cat(sprintf("  ...processed %d/%d peaks\n", i, nrow(RPKM.cqn)))
  cur_data <- as.data.frame(cbind(logRPKM = RPKM.cqn[i, ], Covariates))
  out <- tryCatch(runlme(cur_data, expression_model), error = function(e) NULL)
  if (is.null(out)) {
    cat(sprintf("  Error fitting peak %d (%s) -- setting NA\n", i, rownames(RPKM.cqn)[i]))
    next
  }
  magnitude_mat[i, ] <- out[[1]][-1]
  p_mat[i, ] <- out[[2]][-1, "p-value"]
}
 
n_na <- sum(is.na(p_mat[, 1]))
cat(sprintf("Fitted %d/%d peaks (%d failed to converge).\n",
            nrow(RPKM.cqn) - n_na, nrow(RPKM.cqn), n_na))
 
save(p_mat, magnitude_mat,
     file = file.path(results_dir, "p_magnitude_lmeDiagnosisASD_ASD_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.rda"))
 
pdf(file.path(results_dir, "HistPval_lmeDiagnosisASD_Age_Batch_RegionBA38_RegionBA44.45_tssenrich.score_FRiP.pdf"), height = 5, width = 8)
hist(p_mat[, 1], breaks = 50, xlab = "p-value (Diagnosis)",
     main = "Reproduction: Diagnosis term p-value distribution")
dev.off()
 
fdr <- p.adjust(p_mat[, 1], method = "fdr")
cat(sprintf("Peaks at FDR<0.05: %d | FDR<0.1: %d | FDR<0.2: %d\n",
            sum(fdr < 0.05, na.rm = TRUE),
            sum(fdr < 0.1, na.rm = TRUE),
            sum(fdr < 0.2, na.rm = TRUE)))
 
cat("\nReproduction run complete.")
