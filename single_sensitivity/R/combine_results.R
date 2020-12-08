#! /usr/bin/Rscript

install.packages("here")
library(SuperLearner)
library(cvAUC)

# in docker container
setwd("/home/output")

library(here)

folder_roots <- read.table(here("bash_output", "all_nabs.txt"))[,1]

out <- vector(mode = "list", length = length(folder_roots))
names(out) <- folder_roots
ct <- 0
for(folder_root in folder_roots){
	ct <- ct + 1
	print(folder_root)
	orig_files <- list.files(here("bash_output", folder_root))

	# assumes only a single csv is present
	orig_data_csv <- orig_files[grep(".csv", orig_files)]
	orig_data <- read.csv(here("bash_output", folder_root, orig_data_csv))

	# load learner
	rds_files <- orig_files[grep(".rds", orig_files)]
	sens_files <- rds_files[grep("sens", rds_files)]
	# assumes two files, one starts with cv
	sens_fit <- readRDS(here("bash_output", folder_root, sens_files[2]))
	stopifnot(class(sens_fit) == "SuperLearner")

	# predictions and evaluations
	sens_pred <- predict(sens_fit)
	sens_sl_pred <- as.numeric(sens_pred[[1]])

	# auc + ci
	est_sens_auc <- ci.cvAUC(predictions = sens_sl_pred, labels = orig_data[, "sens"])

	out[[ct]] <- est_sens_auc
}

bnab <- AUC <- cil <- ciu <- NULL
for(i in folder_roots){
	bnab <- c(bnab, i)
	AUC <- c(AUC, out[[i]]$cvAUC)
	cil <- c(cil, out[[i]]$ci[1])
	ciu <- c(ciu, out[[i]]$ci[2])
}

rslt_df <- data.frame(bnab, AUC, cil, ciu)
save(rslt_df, file = here("R_output", "rslt_df.RData"))
