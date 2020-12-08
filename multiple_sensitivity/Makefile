## auc_fig                  : a figure of auc results across all combo nabs
.PHONY: auc_fig
auc_fig: R/make_auc_fig.R R_output/rslt_df.RData R_output/n_sens_df.RData
	Rscript R/make_auc_fig.R

## R_output/rslt_df.RData   : a data.frame of auc (95% CI) across all nabs
R_output/rslt_df.RData: R_output/all_nabs.txt R/combine_results.R slapnap_output
	chmod +x R/combine_results.R && \
	docker run --entrypoint /usr/bin/Rscript -v $$(pwd):/home/output slapnap/slapnap /home/output/R/combine_results.R

## R_output/n_sens_df.RData : data.frame of # sequences sens and res for each combo
R_output/n_sens_df.RData: R_output/all_nabs.txt R/combine_results.R
	Rscript R/describe_combo_sens.R

## slapnap_output           : populates output directory with slapnap results for all nabs.
.PHONY: slapnap_output
slapnap_output: bash/run_docker_containers.sh
	chmod +x bash/run_docker_containers.sh && ./bash/run_docker_containers.sh

## R_output/all_nabs.txt    : a string of output folder names read in by R
R_output/all_nabs.txt: bash/get_nabs_for_R.sh
	chmod +x bash/get_nabs_for_R.sh && ./bash/get_nabs_for_R.sh

.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<