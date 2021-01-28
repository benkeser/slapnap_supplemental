# do the test-set analysis

# -----------------------------------------------------------------------------
# load required functions and packages
# -----------------------------------------------------------------------------
library("here")
library("tibble")
library("dplyr")
library("tidyr")
library("SuperLearner")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

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
cv_learner_files <- all_mab_files[grepl("cvlearner", all_mab_files)]
learner_ic50 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117",
         learner_files[grepl("ic50", learner_files)])
)
cvlearner_ic50 <- readRDS(
    here("test_set", "docker_output", "10-1074_3bnc117",
         cv_learner_files[grepl("ic50", learner_files)])
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

predprob_estsens <- tibble::tibble(probs = preds_estsens, 
                                   sensitive = factor(
                                       case_when(
                                           cc_data$estsens == 1 ~ "Sensitive",
                                           cc_data$estsens == 0 ~ "Resistant"
                                       )
                                   ))
predprob_multsens <- tibble::tibble(probs = preds_multsens, 
                                   sensitive = factor(
                                       case_when(
                                           cc_data$multsens == 1 ~ "Sensitive",
                                           cc_data$multsens == 0 ~ "Resistant"
                                       )
                                   ))

estsens_plot <- predprob_estsens %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of sensitivity",
         x = "") +
    guides(fill = FALSE)
ggsave(here("test_set", "R_output", "pred_estsens.png"),
       estsens_plot, 
       width = 10, height = 10, units = "cm")

multsens_plot <- predprob_multsens %>% 
    ggplot(aes(x = sensitive, y = probs, fill = sensitive)) +
    geom_boxplot() + 
    geom_point(position = position_jitterdodge()) +
    labs(y = "Predicted probability of sensitivity",
         x = "") +
    guides(fill = FALSE)
ggsave(here("test_set", "R_output", "pred_multsens.png"),
       multsens_plot, 
       width = 10, height = 10, units = "cm")

roc_perf_estsens <- ROCR::performance(
    ROCR::prediction(preds_estsens, cc_data$estsens), measure = c("tpr", "fpr")
    )
roc_tib_estsens <- tibble::tibble(xval = roc_perf_estsens@x.values[[1]],
                                  yval = roc_perf_estsens@y.values[[1]])
roc_tib_estsens %>% 
    ggplot(aes(x = xval, y = yval)) +
    geom_step(lwd = 2) + 
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    labs(x = "False positive rate", y = "True positive rate")

roc_perf_multsens <- ROCR::performance(
    ROCR::prediction(preds_multsens, cc_data$multsens), measure = c("tpr", "fpr")
)
roc_tib_multsens <- tibble::tibble(xval = roc_perf_multsens@x.values[[1]],
                                  yval = roc_perf_multsens@y.values[[1]])
roc_tib_multsens %>% 
    ggplot(aes(x = xval, y = yval)) +
    geom_step(lwd = 2) + 
    ylim(c(0, 1)) +
    xlim(c(0, 1)) +
    labs(x = "False positive rate", y = "True positive rate")

perf_ic50 <- vimp::measure_r_squared(fitted_values = preds_ic50, 
                                     y = cc_data$ic50)
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
                         se = vimp::vimp_se(point_est_ic50, perf_ic50$eif))
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
