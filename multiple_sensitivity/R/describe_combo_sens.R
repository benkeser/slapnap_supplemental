library(SuperLearner)
library(cvAUC)
library(here)

folder_roots <- read.table(here("bash_output", "all_nabs.txt"))[,1]

out <- vector(mode = "list", length = length(folder_roots))
names(out) <- folder_roots
ct <- 0
for(folder_root in folder_roots){
	ct <- ct + 1
	print(folder_root)
	folder_root2 <- paste0(folder_root, "_2")
	new_files <- list.files(here("bash_output", folder_root2))

	new_data_csv <- new_files[grep(".csv", new_files)]
	new_data <- read.csv(here("bash_output", folder_root2, new_data_csv))
	out[[ct]] <- list(sens = sum(new_data$sens), res = sum(1 - new_data$sens))
}

nab <- sens <- n_sens <- n_res <- NULL
for(i in folder_roots){
	nab <- c(nab, rep(i, 2))
	sens <- c(sens, c("est", "mult"))
	n_sens <- c(n_sens, rep(out[[i]]$sens, 2))
	n_res <- c(n_res, rep(out[[i]]$res, 2))
}

n_sens_df <- data.frame(nab, sens, n_sens, n_res)
save(n_sens_df, file = here("R_output", "n_sens_df.RData"))