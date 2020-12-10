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
	orig_files <- list.files(here("docker_output", folder_root))

	# assumes only a single csv is present
	orig_data_csv <- orig_files[grep(".csv", orig_files)]
	orig_data <- read.csv(here("docker_output", folder_root, orig_data_csv))

	# load learner
	rds_files <- orig_files[grep(".rds", orig_files)]
	sens_files <- rds_files[grep("sens", rds_files)]
	# assumes two files, first one starts with 'cv'
	sens_fit <- readRDS(here("docker_output", folder_root, sens_files[1]))
	stopifnot(class(sens_fit) == "CV.SuperLearner")

	# predictions and evaluations
	sens_sl_pred <- as.numeric(sens_fit$SL.predict)

	# auc + ci
    ests <- vector("list", length = sens_fit$V)
    ses <- vector("list", length = sens_fit$V)
    for (v in 1:sens_fit$V) {
        est_sens_auc <- ci.cvAUC(predictions = sens_sl_pred[sens_fit$folds[[v]]], labels = sens_fit$Y[sens_fit$folds[[v]]])
        ests[[v]] <- est_sens_auc$cvAUC
        ses[[v]] <- est_sens_auc$se
    }
    est_vec <- unlist(ests)
    se_vec <- unlist(ses)
    mn_est <- mean(est_vec)
    mn_se <- mean(se_vec)
    grad <- 1 / (mn_est - mn_est ^ 2)
    logit_se <- sqrt(mn_se ^ 2 * grad ^ 2)
    cil <- plogis(qlogis(mn_est) - qnorm(0.975) * logit_se)
    ciu <- plogis(qlogis(mn_est) + qnorm(0.975) * logit_se)

	out[[ct]] <- list(AUC = mn_est, ci = c(cil, ciu))
}

bnab <- AUC <- cil <- ciu <- NULL
for(i in folder_roots){
	bnab <- c(bnab, i)
	AUC <- c(AUC, out[[i]]$AUC)
	cil <- c(cil, out[[i]]$ci[1])
	ciu <- c(ciu, out[[i]]$ci[2])
}

rslt_df <- data.frame(bnab, AUC, cil, ciu)
save(rslt_df, file = here("R_output", "rslt_df.RData"))
