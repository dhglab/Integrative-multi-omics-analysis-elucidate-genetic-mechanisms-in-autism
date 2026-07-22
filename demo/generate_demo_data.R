setwd("/Volumes/DataTransferBwMac/Working_Dir/Documents/Geschwind_lab/LAB/Projects/Project_ASD/ATAC-seq/8_ASDandCTL_batch2/data_analysis_20220727_startover/results/5_DiffATAC/")
lnames = load("CQN.rda"); rm(list = setdiff(ls(), "RPKM.cqn"))
str(RPKM.cqn) #  num [1:128036, 1:38] 2.007 0.434 0.868 -1.019 -0.549 ...
save(RPKM.cqn, file = "RPKM_cqn.rda", compress = "xz")

## check
rm(list = ls())
lnames = load("RPKM_cqn.rda")
str(RPKM.cqn) #  num [1:128036, 1:38] 2.007 0.434 0.868 -1.019 -0.549 ...
# Great.
