# do the test-set analysis

# -----------------------------------------------------------------------------
# load required functions and packages
# -----------------------------------------------------------------------------
library("here")
library("tibble")
library("dplyr")
library("tidyr")
library("SuperLearner")
library("ROCR")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

source(here("test_set", "R", "utils.R"))
source(here("test_set", "R", "00_utils.R"))
source(here("test_set", "R", "05_plotting_functions.R"))

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
    select(-seq.id.lanl, -seq.id.catnap, -ic50, -ic80, 
           -starts_with("estsens"), -starts_with("multsens"))

# get the correct names for the cv learners (based on date that you ran SLAPNAP)
# note that the final piece determines the analysis threshold for binary variables
# (2 matches the Mendoza et al. and Bar-On et al. papers; 1 matches other SLAPNAP analyses)
all_mab_files_2 <- list.files(here("test_set", "docker_output", 
                                   "10-1074_3bnc117_2"))
training_data_2 <- readr::read_csv(
    here("test_set", "docker_output", "10-1074_3bnc117_2", 
         all_mab_files_2[grepl("slapnap", all_mab_files_2)])
)
learner_files_2 <- all_mab_files_2[grepl("learner", all_mab_files_2) & 
                                   !grepl("cv", all_mab_files_2)]
all_mab_files_1 <- list.files(here("test_set", "docker_output", 
                                   "10-1074_3bnc117_1"))
training_data_1 <- readr::read_csv(
    here("test_set", "docker_output", "10-1074_3bnc117_1", 
         all_mab_files_1[grepl("slapnap", all_mab_files_1)])
)
learner_files_1 <- all_mab_files_1[grepl("learner", all_mab_files_1) & 
                                       !grepl("cv", all_mab_files_1)]
cv_learner_files_2 <- all_mab_files_2[grepl("cvlearner", all_mab_files_2)]
learner_ic50 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         learner_files_2[grepl("ic50", learner_files_2)])
)
cvlearner_ic50 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         cv_learner_files_2[grepl("ic50", learner_files_2)])
)
learner_estsens_2 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         learner_files_2[grepl("estsens", learner_files_2)])
)
cvlearner_estsens_2 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         cv_learner_files_2[grepl("estsens", cv_learner_files_2)])
)
learner_multsens_2 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         learner_files_2[grepl("multsens", learner_files_2)])
)
cvlearner_multsens_2 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_2",
         cv_learner_files_2[grepl("multsens", cv_learner_files_2)])
)
learner_estsens_1 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_1",
         learner_files_1[grepl("estsens", learner_files_1)])
)
learner_multsens_1 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117_1",
         learner_files_1[grepl("multsens", learner_files_1)])
)

# read in cross-validated fits
fit_list_out <- load_cv_fits(opts = list(outcomes = c("ic50", "sens1",
                                                      "sens2"),
                                         cvperf = TRUE, cvtune = TRUE),
                             run_sls = c(TRUE, TRUE, TRUE),
                             code_dir = paste0(here("test_set", "docker_output",
                                             "10-1074_3bnc117_2"), "/"))

if (!dir.exists(here("test_set", "R_output"))) {
    dir.create(here("test_set", "R_output"))
}

# summaries of CV superlearners
summary(cvlearner_estsens_2, method = "method.AUC", 
        opts = list(learners = c("rf", "lasso", "xgboost")))

# ------------------------------------------------------------------------------
# estimate test-set performance (using AUC, accuracy for binary outcomes; R-squared for IC50)
# ------------------------------------------------------------------------------
preds_ic50 <- predict(learner_ic50, newdata = cc_covars)$pred
preds_estsens_2 <- predict(learner_estsens_2, newdata = cc_covars)$pred
preds_multsens_2 <- predict(learner_multsens_2, newdata = cc_covars)$pred
preds_estsens_1 <- predict(learner_estsens_1, newdata = cc_covars)$pred
preds_multsens_1 <- predict(learner_multsens_1, newdata = cc_covars)$pred

# plot predictions vs outcomes for continuous; plot roc curves for binary
continuous_tib <- tibble::tibble(y = cc_data$ic50, yhat = preds_ic50)
ic50_cor <- cor(cc_data$ic50, preds_ic50)
continuous_plot <- continuous_tib %>% 
    ggplot(aes(x = yhat, y = y)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    ylim(c(-3.5, 2)) +
    xlim(c(-3.5, 2)) +
    labs(y = "Observed log10 IC-50", x = "Predicted log10 IC-50") +
    annotate(geom = "text", x = -2, y = 1, label = sprintf("r = %.2f", ic50_cor),
             color = "red")
ggsave(here("test_set", "R_output", "pred_ic50.png"),
       continuous_plot, 
       width = 10, height = 10, units = "cm")

plot_cv_predictions(fit_list_out$out$ic50, outcome_name = "IC-50")

# predicted probability of sensitivity: threshold 2
predprob_estsens_2 <- tibble::tibble(probs = preds_estsens_2, 
                                   sensitive = factor(
                                       case_when(
                                           cc_data$estsens_2 == 1 ~ "Sensitive",
                                           cc_data$estsens_2 == 0 ~ "Resistant"
                                       )
                                   ))
predprob_multsens_2 <- tibble::tibble(probs = preds_multsens_2, 
                                   sensitive = factor(
                                       case_when(
                                           cc_data$multsens_2 == 1 ~ "Sensitive",
                                           cc_data$multsens_2 == 0 ~ "Resistant"
                                       )
                                   ))

estsens_plot_2 <- predprob_estsens_2 %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of estimated sensitivity",
         x = "") +
    guides(fill = FALSE) +
    ggtitle("Sensitivity threshold: 2")

multsens_plot_2 <- predprob_multsens_2 %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of multiple sensitivity",
         x = "") +
    guides(fill = FALSE) +
    ggtitle("Sensitivity threshold: at least 1 mAb sensitive (at threshold 2)")
ggsave(here("test_set", "R_output", "pred_sens_2.png"),
       plot_grid(estsens_plot_2, multsens_plot_2, labels = "AUTO"), 
       width = 40, height = 15, units = "cm")

# predicted probability of sensitivity: threshold 1
predprob_estsens_1 <- tibble::tibble(probs = preds_estsens_1, 
                                     sensitive = factor(
                                         case_when(
                                             cc_data$estsens_1 == 1 ~ "Sensitive",
                                             cc_data$estsens_1 == 0 ~ "Resistant"
                                         )
                                     ))
predprob_multsens_1 <- tibble::tibble(probs = preds_multsens_1, 
                                      sensitive = factor(
                                          case_when(
                                              cc_data$multsens_1 == 1 ~ "Sensitive",
                                              cc_data$multsens_1 == 0 ~ "Resistant"
                                          )
                                      ))

estsens_plot_1 <- predprob_estsens_1 %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of estimated sensitivity",
         x = "") +
    guides(fill = FALSE) +
    ggtitle("Sensitivity threshold: 1")

multsens_plot_1 <- predprob_multsens_1 %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of multiple sensitivity",
         x = "") +
    guides(fill = FALSE) +
    ggtitle("Sensitivity threshold: at least 1 mAb sensitive (at threshold 1)")
ggsave(here("test_set", "R_output", "pred_sens_1.png"),
       plot_grid(estsens_plot_1, multsens_plot_1, labels = "AUTO"), 
       width = 40, height = 15, units = "cm")

# in-sample predicted prob of resistant, threshold 2
summ_estsens_2 <- subset.summary.myCV.SuperLearner(object = fit_list_out$out$sens1,
                                 method = "method.AUC", 
                                 opts = list(learners = c("rf", "xgboost", "lasso")),
                                 subset_vec = training_data_2$subtype.is.B)
summ_estsens_2$Table 
summ_multsens_2 <- subset.summary.myCV.SuperLearner(object = fit_list_out$out$sens2,
                                                   method = "method.AUC", 
                                                   opts = list(learners = c("rf", "xgboost", "lasso")),
                                                   subset_vec = training_data_2$subtype.is.B)
summ_multsens_2$Table 

cv_estsens_plot_2 <- plot_predicted_prob_boxplots_subset(fit_list_out$out$sens1, 
                                    opts = list(learners = c("rf", "xgboost", "lasso")),
                                    subset_vec = factor(training_data_2$subtype.is.B,
                                                        levels = c(0, 1),
                                                        labels = c("Not subtype B",
                                                                   "Subtype B")),
                                    shape_lab = "Virus subtype")
cv_multsens_plot_2 <- plot_predicted_prob_boxplots_subset(fit_list_out$out$sens2, 
                                                         opts = list(learners = c("rf", "xgboost", "lasso")),
                                                         subset_vec = factor(training_data_2$subtype.is.B,
                                                                             levels = c(0, 1),
                                                                             labels = c("Not subtype B",
                                                                                        "Subtype B")),
                                                         shape_lab = "Virus subtype")
ggsave(here("test_set", "R_output", "cv_sens_2.png"),
       plot_grid(cv_estsens_plot_2, cv_multsens_plot_2, labels = "AUTO"),
       width = 60, height = 20, units = "cm")

# cv performance metrics
perf_ic50 <- vimp::measure_r_squared(fitted_values = preds_ic50, 
                                     y = cc_data$ic50)
mse_ic50 <- mean((preds_ic50 - cc_data$ic50) ^ 2)
var_ic50 <- mean((cc_data$ic50 - mean(cc_data$ic50)) ^ 2)
mse_train_ic50 <- mean((mean(training_data_2$ic50) - cc_data$ic50) ^ 2)
1 - mse_ic50 / var_ic50
1 - mse_ic50 / mse_train_ic50
perf_estsens_2 <- vimp::measure_auc(fitted_values = preds_estsens_2,
                                  y = cc_data$estsens_2,
                                  scale = "logit")
perf_estsens_1 <- vimp::measure_auc(fitted_values = preds_estsens_1,
                                    y = cc_data$estsens_1,
                                    scale = "logit")
perf_multsens_2 <- vimp::measure_auc(fitted_values = preds_multsens_2,
                                   y = cc_data$multsens_2,
                                   scale = "logit")
perf_multsens_1 <- vimp::measure_auc(fitted_values = preds_multsens_1,
                                     y = cc_data$multsens_1,
                                     scale = "logit")
point_est_ic50 <- perf_ic50$point_est
point_est_estsens_2 <- perf_estsens_2$point_est
point_est_multsens_2 <- perf_multsens_2$point_est
point_est_estsens_1 <- perf_estsens_1$point_est
point_est_multsens_1 <- perf_multsens_1$point_est
ci_ic50 <- vimp::vimp_ci(point_est_ic50,
                         se = vimp::vimp_se(point_est_ic50, perf_ic50$eif))
ci_estsens_2 <- vimp::vimp_ci(point_est_estsens_2,
                         se = vimp::vimp_se(point_est_estsens_2, 
                                            perf_estsens_2$eif),
                         scale = "logit")
ci_multsens_2 <- vimp::vimp_ci(point_est_multsens_2,
                         se = vimp::vimp_se(point_est_multsens_2, 
                                            perf_multsens_2$eif),
                         scale = "logit")
ci_estsens_1 <- vimp::vimp_ci(point_est_estsens_1,
                              se = vimp::vimp_se(point_est_estsens_1, 
                                                 perf_estsens_1$eif),
                              scale = "logit")
ci_multsens_1 <- vimp::vimp_ci(point_est_multsens_1,
                               se = vimp::vimp_se(point_est_multsens_1, 
                                                  perf_multsens_1$eif),
                               scale = "logit")
perf_text_ic50 <- sprintf("%.2f [%.2f, %.2f]", max(perf_ic50$point_est, 0), 
                          ci_ic50[1], 
                          ci_ic50[2])
perf_text_estsens_2 <- sprintf("%.2f [%.2f, %.2f]", point_est_estsens_2, 
                               ci_estsens_2[1], ci_estsens_2[2])
perf_text_multsens_2 <- sprintf("%.2f [%.2f, %.2f]", point_est_multsens_2, 
                                ci_multsens_2[1], ci_multsens_2[2])
perf_text_estsens_1 <- sprintf("%.2f [%.2f, %.2f]", point_est_estsens_1, 
                               ci_estsens_1[1], ci_estsens_1[2])
perf_text_multsens_1 <- sprintf("%.2f [%.2f, %.2f]", point_est_multsens_1, 
                                ci_multsens_1[1], ci_multsens_1[2])
output_tib <- tibble::tibble(Outcome = c("IC$_{50}$", "Estimated sensitivity (2)", 
                           "Multiple sensitivity (2)", "Estimated sensitivity (1)",
                           "Multiple sensitivity (1)"),
               Performance = c(perf_text_ic50, perf_text_estsens_2, 
                               perf_text_multsens_2, perf_text_estsens_1,
                               perf_text_multsens_1))
saveRDS(output_tib, here("test_set", "R_output", "test_set_perf.rds"))

# check consensus
all_vars <- apply(cc_covars, 2, var)
sum(all_vars == 0)
sum(all_vars == 0) / ncol(cc_covars)
training_covars <- training_data_1 %>% 
    select(-seq.id.lanl, -seq.id.catnap, -ic50, -ic80, -iip, -estsens, -multsens)
all_training_vars <- apply(training_covars, 2, var)
sum(all_training_vars == 0)
all_one <- colMeans(cc_covars)
sum(all_one == 1)
