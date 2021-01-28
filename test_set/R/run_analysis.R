# do the test-set analysis

# -----------------------------------------------------------------------------
# load required functions and packages
# -----------------------------------------------------------------------------
library("here")
library("tibble")
library("dplyr")
library("tidyr")
library("SuperLearner")

# -----------------------------------------------------------------------------
# read in the dataset (after running pull_outcome_data.R)
# and read in the CV learner
# -----------------------------------------------------------------------------
analysis_dataset <- readRDS(here("test_set", "data", "analysis_dataset.rds"))
cc_data <- analysis_dataset %>% 
    mutate(cc = complete.cases(analysis_dataset)) %>% 
    filter(cc == 1) %>% 
    select(-cc)
cc_covars <- cc_data %>% 
    select(-seq.id.lanl, -seq.id.catnap, -ic50, -ic80, -estsens, -multsens)

# get the correct names for the cv learners (based on date that you ran SLAPNAP)
all_mab_files <- list.files(here("test_set", "docker_output", 
                                   "10-1074_3bnc117"))
learner_files <- all_mab_files[grepl("learner", all_mab_files) & 
                                   !grepl("cv", all_mab_files)]
learner_ic50 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117",
         learner_files[grepl("ic50", learner_files)])
)
learner_estsens <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117",
         learner_files[grepl("estsens", learner_files)])
)
learner_multsens <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117",
         learner_files[grepl("multsens", learner_files)])
)

# ------------------------------------------------------------------------------
# estimate test-set performance (using AUC, accuracy for binary outcomes; R-squared for IC50)
# ------------------------------------------------------------------------------
set.seed(4747)
preds_ic50 <- predict(learner_ic50, newdata = cc_covars)$pred
set.seed(1234)
preds_estsens <- predict(learner_estsens, newdata = cc_covars)$pred
set.seed(5678)
preds_multsens <- predict(learner_multsens, newdata = cc_covars)$pred

perf_ic50 <- vimp::measure_r_squared(fitted_values = preds_ic50, 
                                     y = cc_data$ic50,
                                     scale = "logit")
perf_estsens <- vimp::measure_auc(fitted_values = preds_estsens,
                                  y = cc_data$estsens,
                                  scale = "logit")
perf_multsens <- vimp::measure_auc(fitted_values = preds_multsens,
                                   y = cc_data$multsens,
                                   scale = "logit")
point_est_ic50 <- perf_ic50$point_est
point_est_estsens <- perf_estsens$point_est
point_est_multsens <- perf_multsens$point_est
ci_ic50 <- vimp::vimp_ci(point_est_ic50,
                         se = vimp::vimp_se(point_est_ic50, perf_ic50$eif),
                         scale = "logit")
ci_estsens <- vimp::vimp_ci(point_est_estsens,
                         se = vimp::vimp_se(point_est_estsens, 
                                            perf_estsens$eif),
                         scale = "logit")
ci_multsens <- vimp::vimp_ci(point_est_multsens,
                         se = vimp::vimp_se(point_est_multsens, 
                                            perf_multsens$eif),
                         scale = "logit")
perf_text_ic50 <- sprintf("%.2f [%.2f, %.2f]", max(perf_ic50$point_est, 0), 
                          ci_ic50[1], 
                          ci_ic50[2])
perf_text_estsens <- sprintf("%.2f [%.2f, %.2f]", point_est_estsens, ci_estsens[1], 
                          ci_estsens[2])
perf_text_multsens <- sprintf("%.2f [%.2f, %.2f]", point_est_multsens, ci_multsens[1], 
                          ci_multsens[2])
output_tib <- tibble::tibble(Outcome = c("IC$_{50}$", "Estimated sensitivity", 
                           "Multiple sensitivity"),
               Performance = c(perf_text_ic50, perf_text_estsens, 
                               perf_text_multsens))
if (!dir.exists(here("test_set", "R_output"))) {
    dir.create(here("test_set", "R_output"))
}
saveRDS(output_tib, here("test_set", "R_output", "test_set_perf.rds"))
