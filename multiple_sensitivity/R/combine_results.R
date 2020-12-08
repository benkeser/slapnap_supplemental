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
	folder_root2 <- paste0(folder_root, "_2")
	orig_files <- list.files(here("bash_output", folder_root))
	new_files <- list.files(here("bash_output", folder_root2))

	# assumes only a single csv is present
	orig_data_csv <- orig_files[grep(".csv", orig_files)]
	orig_data <- read.csv(here("bash_output", folder_root, orig_data_csv))

	new_data_csv <- new_files[grep(".csv", new_files)]
	new_data <- read.csv(here("bash_output", folder_root2, new_data_csv))

	if(all(new_data$sens == 0) | all(new_data$sens == 1)){
		tmp_out <- list(cvAUC = NA, ci = c(NA, NA))
		out[[ct]] <- list(est_sens = tmp_out, mult_sens = tmp_out)
	}else{
		orig_min_col_idx <- min(grep("geographic", colnames(orig_data)))
		orig_ft_names <- colnames(orig_data)[orig_min_col_idx:ncol(orig_data)]

		new_min_col_idx <- min(grep("geographic", colnames(new_data)))
		new_ft_names <- colnames(new_data)[new_min_col_idx:ncol(new_data)]

		shared_col <- intersect(orig_ft_names, new_ft_names)
		missing_col_new_dat_idx <- which(!(orig_ft_names %in% new_ft_names))

		new_ft <- new_data[ , shared_col]
		n_new <- nrow(new_data)
		for(i in missing_col_new_dat_idx){
			new_ft[[orig_ft_names[i]]] <- rep(0, n_new)
		}
		new_ft_ordered <- new_ft[,orig_ft_names]

		# load learner
		rds_files <- orig_files[grep(".rds", orig_files)]
		est_sens_files <- rds_files[grep("estsens", rds_files)]
		mult_sens_files <- rds_files[grep("multsens", rds_files)]
		# assumes two files, one starts with cv
		est_sens_fit <- readRDS(here("bash_output", folder_root, est_sens_files[2]))
		stopifnot(class(est_sens_fit) == "SuperLearner")
		mult_sens_fit <- readRDS(here("bash_output", folder_root, mult_sens_files[2]))
		stopifnot(class(mult_sens_fit) == "SuperLearner")

		# predictions and evaluations
		est_sens_pred <- predict(est_sens_fit, newdata = new_ft_ordered)
		est_sens_sl_pred <- as.numeric(est_sens_pred[[1]])

		mult_sens_pred <- predict(mult_sens_fit, newdata = new_ft_ordered)
		mult_sens_sl_pred <- as.numeric(mult_sens_pred[[1]])

		# auc + ci
		est_sens_auc <- ci.cvAUC(predictions = est_sens_sl_pred, labels = new_data$sens)
		mult_sens_auc <- ci.cvAUC(predictions = mult_sens_sl_pred, labels = new_data$sens)

		out[[ct]] <- list(est_sens = est_sens_auc, mult_sens = mult_sens_auc)
	}
}

nab <- sens <- auc <- auc_lower <- auc_upper <- NULL
for(i in folder_roots){
	nab <- c(nab, rep(i, 2))
	sens <- c(sens, c("est", "mult"))
	auc <- c(auc, out[[i]]$est_sens$cvAUC, out[[i]]$mult_sens$cvAUC)
	auc_lower <- c(auc_lower, out[[i]]$est_sens$ci[1], out[[i]]$mult_sens$ci[1])
	auc_upper <- c(auc_upper, out[[i]]$est_sens$ci[2], out[[i]]$mult_sens$ci[2])
}

rslt_df <- data.frame(nab, sens, auc, auc_lower, auc_upper)
save(rslt_df, file = here("R_output", "rslt_df.RData"))