# utility functions for use throughout SLAPNAP

# -------------------------------------------------------------------
# get system variables
# -------------------------------------------------------------------
# read in a system variable
# use boolean for boolean options
# use !boolean for semi-colon-separated list options
get_sys_var <- function(option = "nab", boolean = FALSE){
    read_string <- Sys.getenv(option)
    if (boolean) {
        out <- read_string == "TRUE"
    } else {
        out <- strsplit(read_string, split = ";")[[1]]
        if (length(out) == 0) {
            out <- ""
        }
    }
    return(out)
}

# read in permanent options
get_global_options <- function(options = c("nab", "outcomes", "learners", "cvtune", "cvperf", "nfolds", "combination_method", "binary_outcomes",
                                           "importance_grp", "importance_ind", "report_name", "return",
                                           "sens_thresh", "multsens_nab", "same_subset", "var_thresh"),
                               options_boolean = c(FALSE, FALSE, FALSE, TRUE,
                                                   TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                                                   FALSE, FALSE, TRUE, FALSE)){
    out <- mapply(option = options, boolean = options_boolean,
                  FUN = get_sys_var, SIMPLIFY = FALSE)
    # replace sensitivity with sens1/2 labels
    if (length(out$nab) == 1){
        out$outcomes <- gsub("sens", "sens1", out$outcomes)
    } else {
        out$outcomes <- gsub("estsens", "sens1", out$outcomes)
        out$outcomes <- gsub("multsens", "sens2", out$outcomes)
    }
    return(out)
}

# -------------------------------------------------------------------
# Figure and table captions, printing nice outcomes,
# and text throughout the report
# -------------------------------------------------------------------
# ---------------------------------------
# General-purpose text
# ---------------------------------------
make_nice_outcome <- function(outcome = "ic50") {
    if (outcome == "ic50") {
        if (length(opts$nab) > 1) {
            "estimated IC$_{50}$"
        } else {
            "IC$_{50}$"
        }
    } else if (outcome == "ic80") {
        if (length(opts$nab) > 1) {
            "estimated IC$_{80}$"
        } else {
            "IC$_{80}$"
        }
    } else if (outcome == "iip") {
        "IIP"
    } else if (outcome == "sens1") {
        if(length(opts$nab) > 1){
          "estimated sensitivity"
        }else{
            "sensitivity"
        }
    } else if (outcome == "sens2") {
        "multiple sensitivity"
    } else {
        ""
    }
}
make_xlab <- function(outcome = "ic50", nab = "VRC01", transformation = "none") {
    if (outcome == "ic50") {
        if (transformation == "none") {
            return(bquote(IC[50]~.(nab)))
        } else {
            return(bquote(log[10]*"("*IC[50]~.(nab)*")"))
        }
    } else if (outcome == "ic80") {
        if (transformation == "none") {
            return(bquote(IC[80]~.(nab)))
        } else {
            return(bquote(log[10]*"("*IC[80]~.(nab)*")"))
        }
    } else {
        if (transformation == "none") {
            return(bquote("IIP"~.(nab)))
        } else {
            return(bquote(log[10*"(IIP"~.(nab)*")"]))
        }
    }
}
get_complete_data_text <- function(opts = list(same_subset = TRUE, outcomes = c("ic50", "ic80", "iip", "sens")), ncomplete_ic50 = NA,
                                   ncomplete_ic80 = NA,
                                   ncomplete_ic5080 = NA,
                                   iip_undef = NA){
    ic50_pres <- ("ic50" %in% opts$outcomes) || ((opts$binary_outcomes == "ic50") && (any(grepl("sens", opts$outcomes))))
    ic80_pres <- ("ic80" %in% opts$outcomes) || ((opts$binary_outcomes == "ic80")  && (any(grepl("sens", opts$outcomes))))
    iip_pres <- "iip" %in% opts$outcomes
    if(!opts$same_subset | !(("ic80" %in% opts$outcomes | "iip" %in% opts$outcomes) & length(opts$outcomes) > 1)){
        tmp <- ""
        # if ic50 was present or the binary outcomes are defined using ic50, put that in
        if (ic50_pres) {
            tmp <- paste0(tmp, ncomplete_ic50, " of these sequences had measured IC$_{50}$")
        }
        if (ic80_pres) {
            txt <- paste0(ncomplete_ic80, " of these sequences had measured IC$_{80}$")
            tmp <- paste0(tmp, ifelse(ic50_pres, paste0(" and ", txt), txt))
        }
        if (iip_pres) {
            tmp <- paste0(tmp, ifelse(!ic50_pres && !ic80_pres, "", ", and "), ncomplete_ic5080, ifelse(ic50_pres & ic80_pres, " had both", " had both IC$_{50}$ and IC$_{80}$"), " measured")
        }
    }else{
        tmp <- paste0(ncomplete_ic5080, " of these sequences had IC$_{50}$ and IC$_{80}$ measured.")
        if(iip_undef > 0){
            tmp <- paste0(tmp, " Amongst these sequences, ", iip_undef, " had the same IC$_{50}$ and IC$_{80}$ and thus were excluded from the analysis.",
                          " Since `same_subset` was set to TRUE in the call to `slapnap`, there were a total of ", ncomplete_ic5080 - iip_undef, " sequences included in the analysis of each endpoint")
        }else{
            tmp <- paste0(tmp, " Since `same_subset` was set to TRUE in the call to `slapnap`, this was the total number of sequences included in the analysis of each endpoint")
        }
    }
    return(tmp)
}
# @param outcomes: can specify to only return text for a single outcome
get_num_obs_text <- function(opts, num_obs_fulls = NA, num_obs_reds = NA, n_row_now = NA, outcomes = "all") {
    if (length(unique(n_row_now)) == 1) {
        complete_obs_txt <- paste0("n = ", n_row_now[1])
        full_obs_txt <- paste0("n = ", num_obs_fulls[1])
        redu_obs_txt <- paste0("n = ", num_obs_reds[1])
    } else {
        outcomes_txt <- gsub(".", "", get_comma_sep_outcomes(opts), fixed = TRUE)
        plural_txt <- ", respectively"
        if (length(n_row_now) == 1) {
            ntot <- n_row_now
            tot_postlim <- " for all outcomes "
        } else {
            ntot <- paste0(paste0(n_row_now[-length(n_row_now)], collapse = ", "), " and ", n_row_now[length(n_row_now)])
            tot_postlim <- paste0(" for ", outcomes_txt, plural_txt)
        }
        if (outcomes == "all") {
            if (length(unique(num_obs_fulls)) == 1) {
                full_prelim <- num_obs_fulls[1]
                full_postlim <- " for all outcomes"
            } else {
                full_prelim <- paste0(paste0(num_obs_fulls[-length(num_obs_fulls)], collapse = ", "), " and ", num_obs_fulls[length(num_obs_fulls)])
                full_postlim <- paste0(" for ", outcomes_txt, plural_txt)
            }
            if (length(unique(num_obs_reds)) == 1) {
                redu_prelim <- num_obs_reds[1]
                redu_postlim <- " for all outcomes"
            } else {
                redu_prelim <- paste0(paste0(num_obs_reds[-length(num_obs_reds)], collapse = ", "), " and ", num_obs_reds[length(num_obs_reds)])
                redu_postlim <- paste0(" for ", outcomes_txt, plural_txt)
            }
        } else {
            full_prelim <- num_obs_fulls[grepl(outcomes, names(num_obs_fulls))]
            redu_prelim <- num_obs_reds[grepl(outcomes, names(num_obs_reds))]
            outcomes_txt <- make_nice_outcome(outcomes)
            plural_txt <- ""
            full_postlim <- redu_postlim <- ""
        }
        complete_obs_txt <- paste0("overall n = ", ntot, tot_postlim)
        full_obs_txt <- paste0("n = ", full_prelim, full_postlim)
        redu_obs_txt <- paste0("n = ", redu_prelim, redu_postlim)
    }
    return(list(complete = complete_obs_txt, full = full_obs_txt, redu = redu_obs_txt))
}

# check whether or not vimp or sls were run for estimated/multiple sensitivity
#' @param opts options
#' @param ran_sl did we run sl for the given outcome?
#' @param ran_vimp did we run vimp for the given outcome?
#' @param outcome_nm the outcome of interest ("sens1" or "sens2")
check_sl_vimp_bin <- function(opts, ran_sl = TRUE, ran_vimp = FALSE, outcome_nm = "sens1") {
    if (outcome_nm == "sens1") {
        ifelse(!ran_sl, paste0(". There were too few observations in at least one class for results to be reliable, and thus ", ifelse(length(opts$nab) > 1, "estimated ", ""), "sensitivity is not included in any learning or biological variable importance analyses"), ifelse(!ran_vimp, paste0(". There were too few observations in at least one class for variable importance results to be reliable, and thus ", ifelse(length(opts$nab) > 1, "estimated ", ""), "sensitivity is not included in any biological variable importance analyses"), ""))
    } else {
        ifelse(!ran_sl, ". There were too few observations in at least one class for results to be reliable, and thus multiple sensitivity is not included in any learning or biological variable importance analyses", ifelse(!ran_vimp, ". There were too few observations in at least one class for variable importance results to be reliable, and thus multiple sensitivity is not included in any biological variable importance analyses", "."))
    }
}

# text for estimated IC-50
get_combo_text <- function(opts, ic50_pres = FALSE, ic80_pres = FALSE) {
    if ("additive" %in% opts$combination_method) {
        txt <- paste0("computed based on the additive model of @wagh2016optimal; ",
                       "for $J$ bnAbs, it is computed as \\[ \\mbox{estimated IC} = \\left( \\sum_{j=1}^J \\mbox{IC}_j^{-1} \\right)^{-1} \\ , \\]",
                       " where $\\mbox{IC}_j$ denotes the measured ",
                       paste0(ifelse(ic50_pres, "IC$_{50}$ ", ""),
                              ifelse(ic50_pres & ic80_pres, "or ", ""),
                              ifelse(ic80_pres, "IC$_{80}$ ", ""), collapse = ""),
                       "for bnAb $j$. ")
    } else {
        txt <- paste0("computed based on the Bliss-Hill model of @wagh2016optimal; ",
                       "for $J$ bnAbs, it is computed using Brent's algorithm [@brent1971] as the concentration value $c$ that minimizes \\[ \\lvert f_J(c) - k \\rvert \\ , \\]",
                       " where $k$ denotes the desired neutralization fraction (",
                       paste0(ifelse(ic50_pres, "50%", ""),
                              ifelse(ic50_pres & ic80_pres, " or ", ""),
                              ifelse(ic80_pres, "80%", ""), collapse = ""), "), ",
                          "\\[f_J(c) = 1 - \\prod_{j=1}^J \\{1 - f_j(c, c \\ / \\ J)\\} \\ , \\]",
                      " $f_j(c, c_j) = (c^m) / (\\mbox{IC}_{50,j}^{m} + c_j^m)$,",
                  " $m = \\mbox{log}_{10}(4) / (\\mbox{log}_{10}(\\mbox{IC}_{80,j}) - \\mbox{log}_{10}(\\mbox{IC}_{50,j}))$,",
                  " and $\\mbox{IC}_{50,j}$ and $\\mbox{IC}_{80,j}$ denote the measured IC$_{50}$ and IC$_{80}$ for bnAb $j$, respectively. ")
    }
    txt
}

get_iip_text <- function(opts) {
    est_txt <- ifelse(length(opts$nab) > 1, "estimated ", "")
    suffix <- ifelse(length(opts$nab) > 1, " and estimated IC$_{50}$ and IC$_{80}$ are computed as described above. ", ". ")
    txt <- paste0("IIP [@shen2008dose; @wagh2016optimal] is calculated as ",
                  "\\[ \\frac{10^m}{\\mbox{", est_txt, "IC$_{50}$}^m + 10^m} \ , \\]", "where $m = \\mbox{log}_{10}(4) / (\\mbox{log}_{10}(\\mbox{", est_txt, "IC}_{80}) - \\mbox{log}_{10}(\\mbox{", est_txt, "IC}_{50}))$",
                  suffix,
                  collapse = "")
    txt
}
# get whether or not binary outcomes are defined using IC-50 or IC-80
get_binary_outcome_text <- function(opts) {
    txt <- paste0("IC$_{", ifelse(opts$binary_outcomes == "ic50", 50, 80) ,"}$ < ", opts$sens_thresh)
    txt
}
# get descriptions of outcomes
get_outcome_descriptions <- function(opts, collapse = TRUE){
    # first describe IC-50, IC-80 if present
    tmp_text <- NULL
    ic80_pres <- "ic80" %in% opts$outcomes
    iip_pres <- "iip" %in% opts$outcomes
    sens1_pres <- "sens1" %in% opts$outcomes
    sens2_pres <- "sens2" %in% opts$outcomes
    ic50_pres <- "ic50" %in% opts$outcomes
    binary_outcomes_txt <- get_binary_outcome_text(opts)
    if(length(opts$nab) > 1){
        if(ic50_pres | ic80_pres | iip_pres | sens1_pres){
            tmp <- paste0("Estimated ",
                          ifelse(ic50_pres | iip_pres | ((sens1_pres || sens2_pres) && opts$binary_outcomes == "ic50"), "IC$_{50}$ ", ""),
                          ifelse((ic50_pres & ic80_pres) | iip_pres, "and ", ""),
                          ifelse(ic80_pres | iip_pres | ((sens1_pres || sens2_pres) && opts$binary_outcomes == "ic80"), "IC$_{80}$ ", ""), collapse = "")
            tmp1_5 <- ifelse((ic50_pres & ic80_pres) | iip_pres, "were ", "was ")
            tmp2 <- get_combo_text(opts, ic50_pres || ((sens1_pres || sens2_pres) && opts$binary_outcomes == "ic50"), ic80_pres || ((sens1_pres || sens2_pres) && opts$binary_outcomes == "ic80"))
            tmp_text <- c(tmp_text, paste0(tmp, tmp1_5, tmp2, collapse = ""))
        }
        if(iip_pres){
            tmp <- get_iip_text(opts)
            tmp_text <- c(tmp_text, tmp)
        }
        if(sens1_pres){
            tmp_text <- c(tmp_text, "Estimated sensitivity is defined by the binary indicator that estimated ", binary_outcomes_txt, ". ")
        }
        if(sens2_pres){
            min_num <- min(c(length(opts$nab), opts$multsens_nab))
            tmp_text <- c(tmp_text, "Multiple sensitivity is defined as the binary indicator of having measured ", binary_outcomes_txt, " for at least ", min_num, " bnAb", ifelse(min_num == 1, ".", "s."))
        }
    } else {
        if(iip_pres){
            tmp <- get_iip_text(opts)
            tmp_text <- c(tmp_text, tmp)
        }
        if(sens1_pres | sens2_pres){
            tmp_text <- c(tmp_text, "Sensitivity is defined by the binary indicator that ", binary_outcomes_txt, ".")
        }
    }
    if (collapse) {
        return(paste0(tmp_text, collapse = ""))
    } else {
        return(tmp_text)
    }
}

# get a comma separated list of outcomes for report
get_comma_sep_outcomes <- function(opts){
    if(length(opts$outcomes) > 2){
        tmp <- paste0(paste0(opts$outcomes[-length(opts$outcomes)], collapse = ", "), ", and ", opts$outcomes[length(opts$outcomes)])
    }else if(length(opts$outcomes) == 2){
        tmp <- paste0(opts$outcomes, collapse = " and ")
    }else{
        tmp <- opts$outcomes
    }
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    all_labels <- c("IC$_{50}$", "IC$_{80}$", "IIP", ifelse(length(opts$nab) == 1, "sensitivity", "estimated sensitivity"), "multiple sensitivity")
    for(i in seq_along(all_outcomes)){
        tmp <- gsub(all_outcomes[i], all_labels[i], tmp)
    }
    tmp <- paste0(tmp, ".")
    return(tmp)
}

# ---------------------------------------
# SL text
# ---------------------------------------
# SL learner descriptions
get_learner_descriptions <- function(opts, n_total_ft = NA, n_ft_screen = NA){

    if(length(opts$learners) == 1){
        tmp <- if(opts$learners[1] == "rf"){
            "random forest [@breiman2001]"
        }else if(opts$learners[1] == "xgboost"){
            "extreme gradient boosting [@chen2016]"
        }else if(opts$learners[1] == "lasso"){
            "elastic net regression [@zou2005]"
        }
        if(all(opts$var_thresh == 0) & !opts$cvtune){
            tmp <- paste0(tmp, " with tuning parameters set to their 'default' values.")
        }else if(all(opts$var_thresh == 0) & opts$cvtune){
            tmp <- paste0(tmp, " with tuning parameters selected using a limited grid search and cross-validation.")
        }else if(!all(opts$var_thresh == 0) & !opts$cvtune & opts$learner != "xgboost"){
            if(length(opts$var_thresh) == 1){
                tmp <- paste0(tmp, " with tuning parameters set to their 'default' values.",
                              " Variable pre-screening was applied to ensure all binary features had at least ", opts$var_thresh, " minority variants.",
                              " This constituted a total of ", n_ft_screen, "/", n_total_ft, " features.")
            }else{
                tmp <- paste0(tmp, " with tuning parameters set to their 'default' values.",
                              " Variable pre-screening procedures were applied that ensured that all binary features had at least ", paste0(opts$var_thresh, collapse = ", "), " minority variants.",
                              " This constituted a total of ", paste0(n_ft_screen, "/", n_total_ft, collapse = ", "), " features, respectively.",
                              " The optimal screening approach was selected using cross-validation.")
            }
        }else if(!all(opts$var_thresh == 0) & opts$cvtune & opts$learner != "xgboost"){
            if(length(opts$var_thresh) == 1){
                tmp <- paste0(tmp, " with tuning parameters selected using a limited grid search and cross-validation.",
                              " Variable pre-screening was applied to ensure all binary features had at least ", opts$var_thresh, " minority variants.",
                              " This constituted a total of ", n_ft_screen, "/", n_total_ft, " features.")
            }else{
                tmp <- paste0(tmp, " with tuning parameters selected using a limited grid search and cross-validation.",
                              " Variable pre-screening procedures were also applied that ensured that all binary features had at least ", paste0(opts$var_thresh, collapse = ", "), " minority variants.",
                              " This constituted a total of ", paste0(n_ft_screen, "/", n_total_ft, collapse = ", "), " features, respectively.",
                              " The optimal screening approach was also selected using cross-validation.")
            }
        }else{
            tmp <- paste0(tmp, ".")
        }
    } else {
        lib_label <- NULL
        if("rf" %in% opts$learners){
            lib_label <- c(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "random forest", ifelse(opts$cvtune, "s [@breiman2001] with varied tuning parameters", " [@breiman2001]")))
        }
        if("xgboost" %in% opts$learners){
            if("rf" %in% opts$learners){
                if(!("lasso" %in% opts$learners)){
                    lib_label <- paste0(lib_label, " and ")
                }else{
                    lib_label <- paste0(lib_label, ", ")
                }
            }
            lib_label <- paste0(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "gradient boosted tree", ifelse(opts$cvtune, "s [@chen2016] with varied tuning parameters", " [@chen2016]")))
        }
        if("lasso" %in% opts$learners){
            if("rf" %in% opts$learners | "xgboost" %in% opts$learners){
                lib_label <- paste0(lib_label, " and ")
            }
            lib_label <- paste0(lib_label, paste0(ifelse(opts$cvtune, "several ", ""), "elastic net regression", ifelse(opts$cvtune, "s [@zou2005] with varied tuning parameters", " [@zou2005]")))
        }
        tmp <- paste0("a super learner ensemble [@vanderlaan2007] of ", lib_label, " and intercept-only regression.")
        if(!all(opts$var_thresh == 0)){
            if(length(opts$var_thresh) == 1){
                tmp <- paste0(tmp, " Each algorithm ", ifelse("xgboost" %in% opts$learners, "(excepting xgboost) ", ""), "included a variable pre-screening to ensure all binary features had at least ", opts$var_thresh, " minority variants.",
                              " This constituted a total of ", n_ft_screen, "/", n_total_ft, " features.")
            }else{
                tmp <- paste0(tmp, " Each algorithm ", ifelse("xgboost" %in% opts$learners, "(excepting xgboost) ", ""), "was additionally implemented in combination with variable pre-screening procedures to ensure that all binary features had at least ", paste0(opts$var_thresh, collapse = ", "), " minority variants.",
                              " This constituted a total of ", paste0(n_ft_screen, "/", n_total_ft, collapse = ", "), " features, respectively.")
            }
        }
    }
    return(tmp)
}
# given options and n_row_now (for captions) load appropriate SuperLearner/
# CV.SuperLearner fits and create a list of CV results for all continuous
# valued outcomes ([[1]] of output) and dichotomous outcomes ([[2]] of output)

get_cont_table_cap <- function(opts, V = 2, n_row_ic50 = NA, n_row_ic80 = NA, n_row_iip = NA){
    tmp <- paste0("Estimates of ", V, "-fold cross-validated $R^2$ for predictions of ")
    if ("ic50" %in% opts$outcomes) {
        tmp <- paste0(tmp, "IC$_{50}$")
    }
    if (any(c("ic80", "iip") %in% opts$outcomes)) {
        if (n_row_ic50 != n_row_ic80) {
            tmp <- paste0(tmp, " (n = ", n_row_ic50, ")")
        } else {
            # tmp <- paste0(tmp, ", ")
        }
        if ("ic80" %in% opts$outcomes) {
            tmp <- paste0(tmp, ifelse(!("ic50" %in% opts$outcomes), "IC$_{80}$", ifelse("iip" %in% opts$outcomes, ", IC$_{80}$", " and IC$_{80}$")))
            if(n_row_ic80 != n_row_iip | !("iip" %in% opts$outcomes)){
                tmp <- paste0(tmp, " (n = ", n_row_ic80, ")")
            }
        }
        if("iip" %in% opts$outcomes){
            tmp <- paste0(tmp, ", and IIP (n = ", n_row_iip, ").")
        } else {
            tmp <- paste0(tmp, ".")
        }
    } else {
        tmp <- paste0(tmp, " (n = ", n_row_ic50, ").")
    }
    return(tmp)
}
# table caption for SL library
get_sllibtab_caption <- function(opts = list(learners = "rf", var_thresh = 0)){
    tmp <- NULL
    if(length(opts$learners) == 1){
        tmp <- paste0(tmp, "Algorithms used in the analysis")
    }else{
        tmp <- paste0(tmp, "Algorithms used in the super learner library")
    }
    if(!all(opts$var_thresh == 0)){
        if(length(opts$var_thresh) == 1){
            tmp <- paste0(tmp, ". Variable pre-screening was applied to each algorithm to ensure all binary features had at least ", opts$var_thresh, " minority variants.")
        }else{
            if(length(opts$learners) == 1){
                tmp <- paste0(tmp, ". Variable pre-screening procedures were applied to each algorithm to ensure that all binary features had at least ", paste0(opts$var_thresh, collapse = ", "), " minority variants.",
                              " The optimal screening approach was selected using cross-validation.")
            }else{
                tmp <- paste0(tmp, ". Each algorithm ", ifelse("xgboost" %in% opts$learners, "(excepting xgboost) ", ""), "was additionally implemented in combination with variable pre-screening procedures to ensure that all binary features had at least ", paste0(opts$var_thresh, collapse = ", "), " minority variants.")
            }
        }
    }
    return(tmp)
}
# nice labels for the SL library
relabel_library <- function(library_names = "SL.ranger", opts = list(var_thresh = 0)){
    labels <- list(
      c("SL.ranger.reg", "rf_default"),
      c("SL.ranger.small", "rf_tune1"),
      c("SL.ranger.large", "rf_tune2"),
      c("SL.xgboost.4", "xgboost_default"),
      c("SL.xgboost.2", "xgboost_tune1"),
      c("SL.xgboost.6", "xgboost_tune2"),
      c("SL.xgboost.8", "xgboost_tune3"),
      c("SL.xgboost.12", "xgboost_tune4"),
      c("SL.glmnet.0", "lasso_default"),
      c("SL.glmnet.25", "lasso_tune1"),
      c("SL.glmnet.50", "lasso_tune2"),
      c("SL.glmnet.75", "lasso_tune3"),
      c("SL.mean", "mean")
    )
    screens_included <- any(grepl("_All", library_names)) | any(grepl("_var_thresh", library_names))
    for(i in seq_along(labels)){
        library_names <- gsub(paste0(labels[[i]][1], ifelse(screens_included, "_All", "")),
                              labels[[i]][2], library_names)
        for(j in opts$var_thresh){
            library_names <- gsub(paste0(labels[[i]][1], ifelse(screens_included, paste0("_var_thresh_", j), "")),
                              paste0(labels[[i]][2],"_screen", j), library_names)
        }
    }
    if(!all(opts$var_thresh == 0) & !screens_included){
        tmp <- NULL
        library_names_minus_mean <- library_names[!(library_names %in% c("mean", "Discrete SL", "Super Learner"))]
        library_names_minus_xgboost <- library_names_minus_mean[!grepl("xgboost", library_names_minus_mean)]
        for(i in opts$var_thresh){
            tmp <- c(tmp, paste0(library_names_minus_xgboost, "_screen", i))
        }
        if(any(grepl("xgboost", library_names))){
            tmp <- c(tmp, library_names[grepl("xgboost", library_names)])
        }
        if(any(library_names == "mean")){
            tmp <- c(tmp, "mean")
        }
        library_names <- tmp
    }
    library_names <- gsub("Super Learner", "super learner", library_names)
    library_names <- gsub("Discrete SL", "cv selector", library_names)

    return(library_names)
}
# get the text for prediction importance
#' @param opts options
#' @param imp_df importance data.frame
#' @param n_ft number of features shown
get_importance_text <- function(opts = list(learners = "rf", cvtune = FALSE, cvperf = FALSE), imp_df = NULL, n_ft = 20){
    if (missing(imp_df)) {
        return("")
    }
    algo_with_highest_wt <- imp_df$algo[1]

    # check if super learner
    is_sl <- (length(opts$learners) > 1) | (opts$cvtune | opts$cvperf)
    # check if tuning parameters varied
    is_tuned <- opts$cvtune

    # if random forest
    if(!is_sl & "rf" %in% opts$learners){
        text_out <- paste0("Specifically, random forest impurity importance for a given feature is computed by taking a normalized sum of the decrease in impurity (i.e., Gini index for binary outcomes; mean squared-error for continuous outcomes) over all nodes in the forest at which a split on that feature has been conducted.")
        if(is_tuned){
            text_out <- paste0(text_out, " These measures are shown for the choice of tuning parameters with the best model fit, as chosen by cross-validation.")
        }
    }else if(!is_sl & "xgboost" %in% opts$learners){
        text_out <- paste0("Specifically, xgboost gain importance measures were computed and the top ", n_ft, " features are shown. Gain measures the improvement in accuracy brought by a given feature to the tree branches on which it appears. The essential idea is that before adding a split on a given feature to the branch, there may be some observations that were poorly predicted, while after adding an additional split on this feature, and each resultant branch is more accurate. Gain measures this change in accuracy.")
        if(is_tuned){
            text_out <- paste0(text_out, " These measures are shown for the choice of tuning parameters with the best model fit, as chosen by cross-validation.")
        }
    }else if(!is_sl & "lasso" %in% opts$learners){
        text_out <- paste0("Specifically, lasso variable importance is taken to be the magnitude of the coefficient for the model with $\\lambda$ chosen via cross-validation, and the top ", n_ft, " are shown.")
        if(is_tuned){
            text_out <- paste0(text_out, " These ranks are shown for the choice of alpha that resulted in the best model fit, as chosen by cross-validation.")
        }
        text_out <- paste0(text_out, " Overall, there were ", sum(abs(imp_df$Importance) > 0), " features that had non-zero coefficient in the final fit.")
    }else{
        text_out <- "Specifically, the algorithm with the largest weight in the super learner ensemble was selected and associated variable importance metrics for this algorithm are shown."
        text_out <- paste0(text_out, " In this case, the highest weight was assigned to a `", algo_with_highest_wt, "` algorithm, and thus the variable importance measures presented correspond to ")
        if(algo_with_highest_wt == "rf"){
            text_out <- paste0(text_out, "random forest impurity importance measures. Impurity is computed by taking a normalized sum of the decrease in impurity (i.e., Gini index for binary outcomes; mean squared-error for continuous outcomes) over all nodes in the forest at which a split on that feature has been conducted.")
        }else if(algo_with_highest_wt == "lasso"){
            text_out <- paste0(text_out, "the magnitude of the coefficient for the model with $\\lambda$ chosen via cross-validation.")
            text_out <- paste0(text_out, " Overall, there were ", sum(abs(imp_df$Importance) > 0), " features that had non-zero coefficient in the final lasso fit.")
        }else if(algo_with_highest_wt == "xgboost"){
            text_out <- paste0(text_out, "xgboost gain importance measures were computed and are shown by their rank. Gain measures the improvement in accuracy brought by a given feature to the tree branches on which it appears. The essential idea is that before adding a split on a given feature to the branch, there may be some observations that are poorly predicted, while after adding an additional split on this feature, and each resultant branch is more accurate. Gain measures this change in accuracy.")
        }
    }
    return(text_out)
}
# ---------------------------------------
# Biological importance text
# ---------------------------------------
#' get biological importance table text
#' @param opts options
#' @param cont_nms nice names of continuous outcomes
#' @param bin_nms nice names of binary outcomes
#' @param num_obs_fulls number of total obs for "full" regression (may differ for each outcome)
#' @param num_obs_reds number of total obs for reduced regression (may differ for each outcome)
#' @param n_row_now (may differ for each outcome)
#' @param importance_type "marginal" or "conditional"
get_biological_importance_table_description <- function(opts, any_cont = TRUE, any_dich = TRUE, cont_nms = c("IC$_{50}$", "IC$_{80}$", "IIP"), bin_nms = c("Sensitivity"), num_obs_fulls = 100, num_obs_reds = 100, n_row_now = 100, vimp_threshold = 0.05, importance_type = "marg", any_signif = FALSE) {
    rel_txt <- ifelse(importance_type == "marginal", "the group of geographic confounders", "the remaining features")
    full_func_txt <- ifelse(importance_type == "marginal", "the feature group of interest", "all available features")
    redu_func_txt <- ifelse(importance_type == "marginal", "the group of geographic confounders", "the reduced set of features (defined by removing the feature group of interest)")
    all_obs_txt <- get_num_obs_text(opts, num_obs_fulls, num_obs_reds, n_row_now)
    complete_obs_txt <- all_obs_txt$complete
    full_obs_txt <- all_obs_txt$full
    redu_obs_txt <- all_obs_txt$redu
    correct_outcomes <- ifelse(length(opts$outcomes) == 1, get_comma_sep_outcomes(opts), "each outcome.")
    cont_nms <- switch((any_cont) + 1, NULL, cont_nms)
    bin_nms <- switch((any_dich) + 1, NULL, bin_nms)
    if (any_cont) {
        cont_nms <- ifelse(length(opts$nab) == 1, cont_nms, paste0("estimated ", cont_nms))
    }
    cont_nm_descr <- ifelse(length(cont_nms) > 0, ifelse(length(cont_nms) == 1, paste0("$R^2$ for ", cont_nms), paste0("$R^2$ for continuous outcomes (", paste(cont_nms, collapse = ", "), ")")), "")
    bin_nm_descr <- ifelse(length(bin_nms) > 0, ifelse(length(bin_nms) == 1, paste0("AUC for ", tolower(bin_nms)), paste0("AUC for binary outcomes (", paste(bin_nms, collapse = ", "), ")")), "")
    correct_description <- paste0(" Importance is measured via ", cont_nm_descr, ifelse(length(cont_nms) > 0 & length(bin_nms) > 0, " and ", ""), bin_nm_descr, ".")
    signif_txt <- paste0(" Stars next to ranks denote groups with p-value less than ", vimp_threshold, " from a hypothesis test with null hypothesis of zero importance.")
    descr <- paste0("Ranked ", importance_type, " variable importance of groups relative to ", rel_txt, " for predicting ", correct_outcomes, correct_description, ifelse(any_signif, signif_txt, ""), " (", complete_obs_txt, "; for estimating the prediction functions based on ", full_func_txt, ", ", full_obs_txt, "; for estimating the prediction functions based on ", redu_func_txt, ", ", redu_obs_txt, ")")
    return(descr)
}
# get biological importance text
#' @param opts options
#' @param grp whether or not this is a group description
#' @return character vector with importance text
get_biological_importance_plot_description <- function(opts, grp = TRUE) {
    if (grp) {
        these_opts <- opts$importance_grp
        this_text <- "group"
    } else {
        these_opts <- opts$importance_ind
        this_text <- "feature"
    }
    if (("marg" %in% these_opts) & ("cond" %in% these_opts)) {
        return(paste0("The left-hand plot shows the marginal biological importance of the ", this_text, " relative to the null model with geographic confounders only. The right-hand plot shows the conditional importance of the ", this_text, " relative to all other ", this_text, "s."))
    } else if ("marg" %in% these_opts) {
        return(paste0("The plot shows the marginal biological importance of the ", this_text, " relative to the null model with geographic confounders only."))
    } else if ("cond" %in% these_opts){
        return(paste0("The plot shows the conditional biological importance of the ", this_text, " relative to all other ", this_text, "s."))
    } else {
        return("")
    }
}
# return figure caption for biological importance
#' @param ncomplete number of complete cases
#' @param num_obs_full number of obs used in the "full" regression
#' @param num_obs_red number of obs used in the "reduced" regression
#' @param outcome the outcome (e.g., "ic50")
#' @param grp whether or not this is group importance
biological_importance_figure_caption <- function(ncomplete = NA, num_obs_full = NA, num_obs_red = NA, outcome = "all", grp = TRUE, marg = TRUE, cond = TRUE, opts, vimp_threshold = 0.05, any_signif = list(grp_conditional = FALSE, grp_marginal = FALSE, ind_conditional = FALSE, ind_marginal = FALSE)) {
    outcome_text <- paste0(make_nice_outcome(outcome), ".")
    if (grp) {
        outer_descr <- "Group"
        inner_descr <- "feature group"
        signif_check <- ifelse(marg & cond, any(any_signif[grepl("grp", names(any_signif))]), ifelse(marg, any(any_signif$grp_marginal), any(any_signif$grp_conditional)))
    } else {
        outer_descr <- "Individual"
        inner_descr <- "feature"
        signif_check <- ifelse(marg & cond, any(any_signif[grepl("ind", names(any_signif))]), ifelse(marg, any(any_signif$ind_marginal), any(any_signif$ind_conditional)))
    }
    all_obs_txt <- get_num_obs_text(opts, num_obs_fulls, num_obs_reds, ncomplete, outcomes = outcome)
    complete_obs_txt <- all_obs_txt$complete
    full_obs_txt <- all_obs_txt$full
    redu_obs_txt <- all_obs_txt$redu
    full_func_txt <- ifelse(marg & cond, "s based on all available features and geographic confounders only", ifelse(marg, " based on geographic confounders only", " based on all available features"))
    redu_func_txt <- ifelse(marg & cond, paste0("s based on the reduced set of features (defined by removing the ", inner_descr, " of interest) and the ", inner_descr, " of interest plus geographic confounders"), ifelse(marg, paste0(" based on the ", inner_descr, " of interest plus geographic confounders"), paste0(" based on the reduced set of features (defined by removing the ", inner_descr, " of interest)")))
    signif_txt <- paste0(" and stars denoting p-values less than ", vimp_threshold)
    ci_txt <- paste0(" ", (1 - vimp_threshold) * 100, "\\% confidence intervals", ifelse(signif_check, signif_txt, ""), " are displayed in blue.")
    cap <- paste0(outer_descr, " biological variable importance for predicting ", outcome_text, ci_txt, " (", complete_obs_txt, "; for estimating the prediction function", full_func_txt, ", ", full_obs_txt, "; for estimating the prediction function", redu_func_txt, ", ", redu_obs_txt, ")")
    return(cap)
}
# ------------------------------------------------------------------------------
# Generate plots and tables of CV performance
# ------------------------------------------------------------------------------
# each entry in the output list is a kable that should be properly labeled.
get_cv_outcomes_tables <- function(fit_list_out = NULL, run_sls = NULL, run_sls2 = NULL, opts){
    fit_list <- fit_list_out$out
    V <- fit_list_out$V
    table_list <- lapply(fit_list[!is.na(names(fit_list))], summary.myCV.SuperLearner, opts = opts)
    all_possible_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")
    if (any(is.na(names(fit_list))) || length(fit_list) < length(all_possible_outcomes)) {
        na_sum <- sum(is.na(names(fit_list)))
        na_list <- rep(list(NA), length(all_possible_outcomes) - (length(fit_list) - na_sum))
        name_check_mat <- sapply(names(fit_list), function(y) sapply(all_possible_outcomes, function(x) grepl(y, x)))
        names(na_list) <- all_possible_outcomes[rowSums(name_check_mat, na.rm = TRUE) == 0]
        table_list2 <- c(table_list, na_list)
        table_list3 <- list(ic50 = table_list2$ic50, ic80 = table_list2$ic80, iip = table_list2$iip, sens1 = table_list2$sens1, sens2 = table_list2$sens2)
        table_list <- table_list3
    }

    # re-label
    all_outcomes <- all_possible_outcomes
    all_labels <- c("IC$_{50}$", "IC$_{80}$", "IIP", ifelse(length(opts$nab) == 1, "Sensitivity", "Estimated sensitivity"), "Multiple sensitivity")
    tmp <- all_outcomes
    for(i in seq_along(all_outcomes)){
        tmp <- gsub(all_outcomes[i], all_labels[i], tmp)
    }
    # now format continuous outcomes table
    sls_run <- (all_possible_outcomes %in% opts$outcomes) & run_sls
    cont_idx <- which((opts$outcomes %in% c("ic50", "ic80", "iip")) & run_sls2)
    rsq_kab <- NULL

    if(length(cont_idx) > 0){
        list_rows <- sapply(cont_idx, get_est_and_ci, fit_list = table_list[!is.na(table_list)], Rsquared = TRUE, simplify = FALSE)
        rsqtab <- Reduce(rbind, lapply(list_rows, unlist, use.names = FALSE))
        if(is.null(dim(rsqtab))) rsqtab <- matrix(rsqtab, nrow = 1)
        row.names(rsqtab) <- tmp[!is.na(table_list)][cont_idx]
        rsq_kab <- knitr::kable(rsqtab, col.names = c("CV-R$^2$", "Lower 95% CI", "Upper 95% CI"),
              digits = 3, row.names = TRUE,
              caption = get_cont_table_cap(opts, V, fit_list_out$n_row_ic50, fit_list_out$n_row_ic80, fit_list_out$n_row_iip))
    }
    # now format dichotomous outcomes table
    dich_idx <- which((opts$outcomes %in% c("sens1", "sens2")) & run_sls2)
    auc_kab <- NULL
    if(length(dich_idx) > 0){
        list_rows <- sapply(dich_idx, get_est_and_ci, fit_list = table_list[!is.na(table_list)], Rsquared = FALSE, simplify = FALSE)
        auctab <- Reduce(rbind, lapply(list_rows, unlist, use.names = FALSE))
        if(is.null(dim(auctab))) auctab <- matrix(auctab, nrow = 1)
        row.names(auctab) <- tmp[!is.na(table_list)][dich_idx]
        auc_kab <- knitr::kable(auctab, col.names = c("CV-AUC", "Lower 95% CI", "Upper 95% CI"),
              digits = 3, row.names = TRUE,
              caption = paste0("Estimates of ", V, "-fold cross-validated AUC for predictions of ", ifelse(length(dich_idx) == 1, tolower(tmp[!is.na(table_list)][dich_idx]), "the binary-valued outcomes"), " (n = ", fit_list_out$n_row_ic50, ")."))
    }
    return(list(r2 = rsq_kab, auc = auc_kab))
}

# load cv_fits for given set of opts, needed since the naming convention
# is different if length(opts$learners) == 1 and opts$cvtune == FALSE
load_cv_fits <- function(opts, run_sls = NULL, code_dir = "/home/lib/"){
    if(!opts$cvtune & !opts$cvperf){
        stop("no cross-validated fit for these options")
    }
    out_list <- vector(mode = "list", length = length(opts$outcomes))
    all_outcomes <- c("ic50", "ic80", "iip", "sens1", "sens2")[run_sls]
    all_file_labels <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")[run_sls]
    # if only one learner used with no cvtuning and only one variable screen, then super learner was fit to assess CV performance
    if(length(opts$learners) == 1 & length(opts$var_thresh) == 1 & !opts$cvtune){
        ct <- 0
        for(i in seq_along(all_outcomes)){
            if(all_outcomes[i] %in% opts$outcomes){
                ct <- ct + 1
                out_list[[ct]] <- readRDS(paste0(code_dir, "fit_", all_file_labels[i], ".rds"))
                class(out_list[[ct]]) <- c("myCV.SuperLearner", class(out_list[[ct]]))
                names(out_list)[ct] <- all_outcomes[i]
            }
        }
    # other wise cv super learner was used to assess CV performance
    }else{
        ct <- 0
        for(i in seq_along(all_outcomes)){
            if(all_outcomes[i] %in% opts$outcomes){
                ct <- ct + 1
                out_list[[ct]] <- readRDS(paste0(code_dir, "cvfit_", all_file_labels[i], ".rds"))
                class(out_list[[ct]]) <- c("myCV.SuperLearner", class(out_list[[ct]]))
                names(out_list)[ct] <- all_outcomes[i]
            }
        }
    }
    if("SuperLearner" %in% class(out_list[[1]])){
        V <- length(out_list[[1]]$validRows)
    }else{
        V <- out_list[[1]]$V
    }
    # figure out number
    n_row_ic50 <- if("ic50" %in% opts$outcomes){
        length(out_list$ic50$Y)
    } else if ("sens1" %in% opts$outcomes) {
        length(out_list$sens1$Y)
    } else if ("sens2" %in% opts$outcomes) {
        length(out_list$sens2$Y)
    } else {
        NULL
    }
    n_row_ic80 <- if("ic80" %in% opts$outcomes){
        length(out_list$ic80$Y)
    }else{
        n_row_ic50
    }
    n_row_iip <- if("iip" %in% opts$outcomes){
        length(out_list$iip$Y)
    }else{
        n_row_ic80
    }

    return(list(out = out_list, V = V, n_row_ic50 = n_row_ic50,
                n_row_ic80 = n_row_ic80, n_row_iip = n_row_iip))
}
# for a given outcome make a panel histogram of the individual
# nabs and a summary table
get_individual_nab_summaries <- function(outcome = "ic50", opts, dat = NULL){
    out_hist <- list()
    out_summary <- list()
    # re-label and grab correct column
    if (outcome == "ic50") {
        outcome_label <- "IC$_{50}$"
        name_prefix <- "nab_"
        name_postfix <- ".ic50.imputed"
    } else if (outcome == "ic80") {
        outcome_label <- "IC$_{80}"
        name_prefix <- "nab_"
        name_postfix <- ".ic80.imputed"
    } else {
        outcome_label <- "IIP"
        name_prefix <- ""
        name_postfix <- ""
    }
    ct <- 0
    for(i in seq_along(opts$nab)){
        ct <- ct + 1
        this_name <- gsub("-", ".", paste0(name_prefix, opts$nab[i], name_postfix))
        out_hist[[ct]] <- make_hist_plot(dat, var_name = this_name,
                                          x_lab = make_xlab(outcome, opts$nab[i], "none"),
                                          y_lab = "Density")
        tmp_sum <- summary(dat[, this_name])[1:6] # to ignore NA columns
        tmp_sum <- c(tmp_sum[1:3], 10^mean(log10(dat[, this_name])), tmp_sum[4:6])
        names(tmp_sum)[4] <- "Geom. Mean"
        out_summary[[i]] <- tmp_sum
        ct <- ct+1
        dat[,paste0("log10_",this_name)] <- log10(dat[, this_name])
        out_hist[[ct]] <- make_hist_plot(dat, var_name = paste0("log10_",this_name), x_lab = make_xlab(outcome, opts$nab[i], "log10"),
                                          y_lab = "")
    }
    return(list(hist = out_hist, summary = out_summary))
}
# -----------------------------------------------------------------------------
# Summarize CV.SuperLearner objects and extract point estimates, CIs
# -----------------------------------------------------------------------------
# this function takes as input a fitted object EITHER of class SuperLearner OR
# of class CV.SuperLearner and computes a summary table of performance.
# if all object$Y are 0/1, it will compute AUC; otherwise it computes R^2
# if SuperLearner is used to evaluate CV performance of a single algorithm,
# a single row table is returned with that algorithms performance.
# if CV.SuperLearner is used to evaluate CV performance of a CV-tuned single algorithm,
# a table with a row for each choice of tuning parameters and for the cv-selected tuning
# parameters is returned.
# if CV.SuperLearner is used to evaluated CV performance of SuperLearner, then an
# additional row is added that describes the performance of SuperLearner.
summary.myCV.SuperLearner <- function (object, obsWeights = NULL, method = NULL, opts, ...) {
    if ("env" %in% names(object)) {
        env = object$env
    }else {
        env = parent.frame()
    }

    is_sl <- "SuperLearner" %in% class(object)
    is_cvsl <- "CV.SuperLearner" %in% class(object)
    if(is_sl | is_cvsl){
      library.names <- object$libraryNames
      if(is_cvsl){
        V <- object$V
      }else{
        V <- length(object$validRows)
      }
      n <- length(object$SL.predict)
      if (is.null(obsWeights)) {
        obsWeights <- rep(1, length(object$Y))
      }

      if(is_cvsl){
        folds <- object$folds
      }else if(is_sl){
        folds <- object$validRows
      }

      if(is_cvsl){
        # only will use this if multiple learners selected
        SL.predict <- object$SL.predict
        # this will be only output if single learner used and opts$cvtune
        discreteSL.predict <- object$discreteSL.predict
        # only will use this if multiple learners selected
        library.predict <- object$library.predict
      }else if(is_sl){
        # in this case a single "default" learner was requested
        # so we can pull Z out from the object
        SL.predict <- object$Z[,1]
      }

      Y <- object$Y
      Risk.SL <- rep(NA, length = V)
      se.SL <- rep(NA, length = V)
      if(is_cvsl){
        Risk.dSL <- rep(NA, length = V)
        se.dSL <- rep(NA, length = V)
        Risk.library <- matrix(NA, nrow = length(library.names),
            ncol = V)
        se.library <- matrix(NA, nrow = length(library.names),
            ncol = V)
        rownames(Risk.library) <- library.names
      }
      if (!(all(Y %in% c(0,1)))) {
          for (ii in seq_len(V)) {
              Risk.SL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                  SL.predict[folds[[ii]]])^2)
              if(is_cvsl){
                Risk.dSL[ii] <- mean(obsWeights[folds[[ii]]] * (Y[folds[[ii]]] -
                    discreteSL.predict[folds[[ii]]])^2)
                Risk.library[, ii] <- apply(library.predict[folds[[ii]],
                    , drop = FALSE], 2, function(x) mean(obsWeights[folds[[ii]]] *
                    (Y[folds[[ii]]] - x)^2))
              }
          }
          if_sl <- (Y - SL.predict)^2 - mean((Y - SL.predict)^2)
          if(is_cvsl){
            if_dsl <- (Y - discreteSL.predict)^2 - mean((Y - discreteSL.predict)^2)
            if_library <- apply(library.predict, 2, function(x){ (Y - x)^2 - mean((Y - x)^2) })
          }
          if_varY <- (Y - mean(Y))^2 - mean((Y - mean(Y))^2)
          get_log_se <- function(if_risk, if_varY, risk, varY,
                                 n = length(if_risk)){
              grad <- matrix(c(1 / risk, - 1 /varY), nrow = 2)
              Sig <- cov(cbind(if_risk, if_varY))
              se_log <- t(grad) %*% Sig %*% grad
              return(se_log)
          }

          if(is_cvsl){
            if(length(opts$learners) > 1){
              se <- (1/sqrt(n)) * c(
                get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y)),
                get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
                mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                       function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
              )
              Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
              library.names), Ave = c(1 - mean(Risk.SL)/var(Y), 1 - mean(Risk.dSL)/var(Y),
              apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se, Min = c(min(1 - Risk.SL/var(Y)),
              min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})), Max = c(max(1 - Risk.SL/var(Y)),
              max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))

            }else{
              se <- (1/sqrt(n)) * c(
                get_log_se(if_risk = if_dsl, if_varY = if_varY, risk = mean(Risk.dSL), varY = var(Y)),
                mapply(if1 = split(if_library, col(if_library)), risk = split(Risk.library, row(Risk.library)),
                       function(if1, risk){ get_log_se(if_risk = if1, if_varY = if_varY, risk = mean(risk), varY = var(Y))})
              )

              Table <- data.frame(Algorithm = c("Discrete SL",
              library.names), Ave = c(1 - mean(Risk.dSL)/var(Y),
              apply(Risk.library, 1, function(x){ 1 - mean(x)/var(Y) })), log_se = se,
              Min = c(min(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ min(1 - mean(x)/var(Y))})),
              Max = c(max(1 - Risk.dSL/var(Y)), apply(Risk.library, 1, function(x){ max(1 - mean(x)/var(Y)) })))
            }
          }else{
            se <- (1/sqrt(n)) * get_log_se(if_risk = if_sl, if_varY = if_varY, risk = mean(Risk.SL), varY = var(Y))
              Table <- data.frame(Algorithm = c(library.names[1]), Ave = c(1 - mean(Risk.SL)/var(Y)),
                                  log_se = se,
                                  Min = c(min(1 - Risk.SL/var(Y))),
                                  Max = c(max(1 - Risk.SL/var(Y))))
          }
      }else {
          requireNamespace("cvAUC")
          for (ii in seq_len(V)) {
            sl_auc <- cvAUC::ci.cvAUC(predictions = SL.predict[folds[[ii]]],
                labels = Y[folds[[ii]]], folds = NULL)
            Risk.SL[ii] <- sl_auc$cvAUC
            se.SL[ii] <- sl_auc$se
            if(is_cvsl){
              dsl_auc <- cvAUC::ci.cvAUC(predictions = discreteSL.predict[folds[[ii]]],
                  labels = Y[folds[[ii]]], folds = NULL)
              Risk.dSL[ii] <- dsl_auc$cvAUC
              se.dSL[ii] <- dsl_auc$se
              library_auc <- apply(library.predict[folds[[ii]], , drop = FALSE], 2, function(x){
                  tmp <- cvAUC::ci.cvAUC(predictions = x, labels = Y[folds[[ii]]], folds = NULL)
                  return(c(tmp$cvAUC, tmp$se))
                })
              Risk.library[,ii] <- library_auc[1,]
              se.library[,ii] <- library_auc[2,]
            }
          }
          if(is_cvsl){
            if(length(opts$learners) > 1){
              se <- c(mean(se.SL, na.rm = TRUE), mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
              Table <- data.frame(Algorithm = c("Super Learner", "Discrete SL",
                        library.names), Ave = c(mean(Risk.SL), mean(Risk.dSL),
                        apply(Risk.library, 1, mean)), se = se, Min = c(min(Risk.SL),
                        min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(max(Risk.SL),
                        max(Risk.dSL), apply(Risk.library, 1, max)))
            }else{
              se <- c(mean(se.dSL, na.rm = TRUE),
                  rowMeans(se.library, na.rm = TRUE))
              Table <- data.frame(Algorithm = c("Discrete SL",
                        library.names), Ave = c(mean(Risk.dSL),
                        apply(Risk.library, 1, mean)), se = se, Min = c(
                        min(Risk.dSL), apply(Risk.library, 1, min)), Max = c(
                        max(Risk.dSL), apply(Risk.library, 1, max)))
            }
          }else{
            se <- c(mean(se.SL, na.rm = TRUE))
              Table <- data.frame(Algorithm = c(library.names[1]),
                                  Ave = c(mean(Risk.SL)), se = se,
                                  Min = c(min(Risk.SL)),
                                  Max = c(max(Risk.SL)))
          }
      }
      out <- list(call = object$call, method = method, V = V, Table = Table)
    }
    class(out) <- "summary.myCV.SuperLearner"
    return(out)
}
# get point estimate and CI
get_est_and_ci <- function(idx = 1, fit_list = NULL, Rsquared = FALSE, constant = qnorm(0.975)){
  cv_fit_table <- fit_list[[idx]]
  Mean <- cv_fit_table$Table$Ave
  if(Rsquared){
    se <- cv_fit_table$Table$log_se
    Lower <- 1 - exp( log(-Mean + 1) + constant * se)
    Upper <- 1 - exp( log(-Mean + 1) - constant * se)
  }else{
    se <- cv_fit_table$Table$se
    # put AUC CI on logit scale
    grad <- 1 / (Mean - Mean^2)
    logit_se <- sqrt(se^2 * grad^2)
    Lower <- plogis(qlogis(Mean) - constant * logit_se); Upper <- plogis(qlogis(Mean) + constant * logit_se)
  }
  return(list(est = Mean[1], ci = c(Lower[1], Upper[1])))
}

# -----------------------------------------------------------------------------
# Check to see whether there are enough outcomes in a given class for
# reliable biological variable importance
# -----------------------------------------------------------------------------
# check dichotomous outcomes to make sure there are enough observations in each class to do SLs and vimp
check_outcome <- function(dat, outcome_nm, V) {
    if (grepl("dichot", outcome_nm)) {
        outcome_tbl <- table(dat[, outcome_nm])
        run_sl <- TRUE
        run_vimp <- TRUE
        if (any(outcome_tbl <= V)) { # need to be able to split into outer folds with at least one of each class per fold
            run_sl <- FALSE
        }
        if (any(outcome_tbl <= 3 * V)) { # need to be able to split into outer folds (half of the data in each, for hypothesis testing); then split into outer folds with at least V of each class per fold; then split into inner folds with at least one of each class per fold
            run_vimp <- FALSE
        }
    } else {
        run_sl <- TRUE
        run_vimp <- TRUE
        outcome_tbl <- length(dat[, outcome_nm])
    }
    return(list(run_sl = run_sl, run_vimp = run_vimp, num_obs = outcome_tbl))
}
check_outcomes <- function(dat, outcome_names, V) {
    checked_outcomes <- lapply(as.list(outcome_names), function(x) check_outcome(dat, x, V))
    names(checked_outcomes) <- outcome_names
    all_outcome_names <- c("log10.pc.ic50", "log10.pc.ic80", "iip", "dichotomous.1", "dichotomous.2")
    all_other_outcomes <- all_outcome_names[!(all_outcome_names %in% outcome_names)]
    checked_other_outcomes <- sapply(all_other_outcomes,
                                     FUN = function(x) list(run_sl = FALSE, run_vimp = FALSE, num_obs = NA),
                                     simplify = FALSE)
    names(checked_other_outcomes) <- all_other_outcomes
    run_sls <- unlist(lapply(c(checked_outcomes, checked_other_outcomes), function(x) x[1]))[c("log10.pc.ic50.run_sl", "log10.pc.ic80.run_sl", "iip.run_sl", "dichotomous.1.run_sl", "dichotomous.2.run_sl")]
    run_vimps <- unlist(lapply(c(checked_outcomes, checked_other_outcomes), function(x) x[2]))[c("log10.pc.ic50.run_vimp", "log10.pc.ic80.run_vimp", "iip.run_vimp", "dichotomous.1.run_vimp", "dichotomous.2.run_vimp")]
    return(list(run_sl = run_sls, run_vimp = run_vimps))
}
# -----------------------------------------------------------------------------
# Get SL- or vimp-specific options from the full options list
# -----------------------------------------------------------------------------
# determine SL options based on outcome name
get_sl_options <- function(outcome_name = "ic50", V = 2) {
    if (grepl("dichot", outcome_name)) {
        sl_fam <- binomial()
        cv_ctrl_lst <- list(V = V, stratifyCV = TRUE)
        sl_method <- "tmp_method.CC_nloglik"
    } else {
        sl_fam <- gaussian()
        cv_ctrl_lst <- list(V = V)
        sl_method <- "tmp_method.CC_LS"
    }
    return(list(fam = sl_fam, ctrl = cv_ctrl_lst, method = sl_method))
}
# Determine vimp options based on outcome name
get_vimp_options <- function(outcome_name = "ic50") {
    if (grepl("dichot", outcome_name)) {
        vimp_measure <- "auc"
    } else {
        vimp_measure <- "r_squared"
    }
    return(list(vimp_measure = vimp_measure))
}

# -----------------------------------------------------------------------------
# Functions for cross-fitted variable importance
# -----------------------------------------------------------------------------
# make outer folds for VIM hypothesis test (based on sample splitting)
make_folds <- function(y = NULL, V = 2, stratified = TRUE) {
  if (stratified) {
    y_1 <- y == 1
    y_0 <- y == 0
    folds_1 <- rep(seq_len(V), length = sum(y_1))
    folds_1 <- sample(folds_1)
    folds_0 <- rep(seq_len(V), length = sum(y_0))
    folds_0 <- sample(folds_0)
    folds <- vector("numeric", length(y))
    folds[y_1] <- folds_1
    folds[y_0] <- folds_0
  } else {
    folds <- rep(seq_len(V), length = length(y))
    folds <- sample(folds)
  }
  return(folds)
}
# get cv folds from a list created by CV.SuperLearner
get_cv_folds <- function(folds_lst = NULL) {
    V <- length(folds_lst)
    v_lst <- sapply(1:V, function(s) rep(s, length(folds_lst[[s]])), simplify = FALSE)
    joint_lst <- mapply(list, v_lst, folds_lst, SIMPLIFY = FALSE)
    folds_mat <- do.call(rbind, lapply(joint_lst, function(x) cbind(x[[1]], x[[2]])))
    folds <- folds_mat[order(folds_mat[, 2]), 1]
    return(folds)
}
# Make list of vimp objects
make_vimp_list <- function(var_groups = NULL, var_inds = NULL) {
    list_names <- c("grp_conditional", "grp_marginal", "ind_conditional", "ind_marginal")
    lst <- sapply(list_names, function(x) NULL, simplify = FALSE)
    return(lst)
}
# Make lists of cv objects
make_cv_lists <- function(folds_lst = NULL, full_vec = NULL, redu_vec = NULL) {
    folds <- get_cv_folds(folds_lst)
    # make lists of the fitted values
    full_lst <- lapply(as.list(1:length(unique(folds))), function(x) full_vec[folds == x])
    redu_lst <- lapply(as.list(1:length(unique(folds))), function(x) redu_vec[folds == x])
    return(list(folds = folds, full_lst = full_lst, redu_lst = redu_lst))
}

# Nice group names for vimp
vimp_nice_group_names <- function(nm_vec = NULL) {
    nice_names <- c("Cysteine counts", "Viral geometry", "Region-specific counts of PNG sites", "gp120 CD4 binding sites", "gp120 V2", "gp120 V3", "gp41 MPER")
    reference_nm_vec <- c("cysteines", "geometry", "glyco", "gp120_cd4bs", "gp120_v2", "gp120_v3", "gp41_mper")
    reference_positions <- apply(as.matrix(nm_vec), 1, function(x) grep(x, reference_nm_vec))
    return(nice_names[reference_positions])
}
vimp_nice_ind_names <- function(nm_vec = NULL) {
    no_hxb2 <- gsub("hxb2.", "", nm_vec, fixed = TRUE)
    no_1mer <- gsub(".1mer", "", no_hxb2, fixed = TRUE)
    return(no_1mer)
}
# nice plotting names
vimp_plot_name <- function(vimp_str = NULL, one_nab = TRUE) {
    plot_nms <- rep(NA, length(vimp_str))
    plot_nms[grepl("iip", vimp_str)] <- "IIP"
    plot_nms[grepl("pc.ic50", vimp_str)] <- paste0(ifelse(one_nab, "", "Estimated "), "IC$_{50}$")
    plot_nms[grepl("pc.ic80", vimp_str)] <- paste0(ifelse(one_nab, "", "Estimated "), "IC$_{80}$")
    plot_nms[grepl("dichotomous.1", vimp_str)] <- ifelse(one_nab, "Sensitivity", "Estimated sensitivity")
    plot_nms[grepl("dichotomous.2", vimp_str)] <- "Multiple sensitivity"
    return(plot_nms)
}
vimp_plot_name_expr <- function(vimp_str = NULL, one_nab = TRUE) {
    est_fillin <- ifelse(one_nab, "", "Estimated")
    if (grepl("iip", vimp_str)) {
        return(bquote("IIP: "))
    } else if (grepl("pc.ic50", vimp_str)) {
        return(bquote(.(est_fillin) ~ IC[50]*": "))
    } else if (grepl("pc.ic80", vimp_str)) {
        return(bquote(.(est_fillin) ~ IC[80]*": "))
    } else if (grepl("dichotomous.1", vimp_str)) {
        if (one_nab) {
            return(bquote("Sensitivity: "))
        } else {
            return(bquote("Estimated sensitivity: "))
        }
    } else {
        return(bquote(.("Multiple sensitivity: ")))
    }
}
vimp_plot_type_expr <- function(str = NULL) {
    grp_nm <- rev(unlist(strsplit(str, "_", fixed = TRUE)))
    nice_grp <- gsub("grp", "Group Variable Importance", grp_nm)
    nice_grp_ind <- gsub("ind", "Individual Variable Importance", nice_grp)
    nice_grp_ind_marg <- gsub("marginal", "Marginal ", nice_grp_ind)
    nice_grp_ind_marg_cond <- gsub("conditional", "Conditional ", nice_grp_ind_marg)
    return(bquote(.( paste(nice_grp_ind_marg_cond, collapse = ""))))
}
vimp_nice_rownames <- function(vimp_obj = NULL, cv = FALSE) {
    mat_s <- vimp_obj$mat$s
    lst_s <- vimp_obj$s
    indx_mat <- sapply(1:length(mat_s), function(x) which(mat_s[x] == lst_s))
    paste_ind <- 3
    if (cv) {
        paste_ind <- 4
    }
    tmp_nms <- unlist(lapply(strsplit(names(lst_s), "_", fixed = TRUE), function(x) paste(x[paste_ind:length(x)], collapse = "_")))
    return(tmp_nms[indx_mat])
}

# -----------------------------------------------------------------------------
# Compute Bliss-Hill predictions of combination neutralization
# -----------------------------------------------------------------------------

# compute individual mAb neutralization curve based on concentration and IC-50
# @param conc the concentration
# @param ic50 the IC-50 value
# @param m the slope (in our analyses, taken to be log(4) / [log(IC-80) - log(IC-50)])
individual_mab_neutralization <- function(conc, mab_conc, ic50, m) {
  f_mab_c <- (conc ^ m) / (ic50 ^ m + mab_conc ^ m)
  f_mab_c
}
# compute neutralization curves using the Bliss-Hill model (with independence assumption)
# @param conc the concentration
# @param ic50 vector or matrix of IC-50 values (columns are mAbs)
# @param ic80 vector or matrix of IC-80 values (columns are mAbs)
bliss_hill_predictions <- function(conc, ic50, ic80) {
  m <- log10(4) / (log10(ic80) - log10(ic50))
  if (!is.null(ncol(ic50))) {
    num_indices <- ncol(ic50)
    this_function <- function(indx, conc) {
      individual_mab_neutralization(conc = conc, mab_conc = conc / ncol(ic50), ic50[, indx], m[, indx])
    }
  } else {
    num_indices <- length(ic50)
    this_function <- function(indx, conc) {
      individual_mab_neutralization(conc = conc, mab_conc = conc / length(ic50), ic50[indx], m[indx])
    }
  }
  one_minus_fs <- do.call(cbind, sapply(1:num_indices, function(indx) 1 - this_function(indx, conc),
                         simplify = FALSE))
  f_c <- 1 - apply(one_minus_fs, 1, prod)
  f_c
}

# back-solve to obtain combination IC-50 or IC-80 from the Bliss-Hill model (with independence assumption)
# @param conc the concentration of interest
# @param ic50 the IC-50 values for all mAbs of interest (a vector)
# @param ic80 the IC-80 values for all mAbs of interest (a vector)
# @return the predicted IC (IC-50 if conc = 0.5, IC-80 if conc = 0.8)
predict_bh_concentration <- function(conc = 0.5, ic50, ic80) {
  optim_func <- function(another_conc, ic50, ic80, conc) {
    abs(bliss_hill_predictions(another_conc, ic50, ic80) - conc)
  }
  if (conc == 0.5) {
    init_par <- 1 / sum(1 / ic50)
    upper_lim <- switch(all(is.na(ic50)) + 1, min(ic50, na.rm = TRUE), NA)
  } else if (conc == 0.8) {
    init_par <- 1 / sum(1 / ic80)
    upper_lim <- switch(all(is.na(ic80)) + 1, min(ic80, na.rm = TRUE), NA)
  } else {
    init_par <- 0
    upper_lim <- 200
  }
  # if any individual one is NA, the whole thing should be
  if (any(is.na(c(ic50, ic80))) | is.na(upper_lim)) {
    ret <- NA
  } else { # do the optimization
    suppressWarnings(optimized <- optim(init_par, fn = optim_func, ic50 = ic50, ic80 = ic80, conc = conc,
                                        method = "Brent", lower = 0, upper = upper_lim,
                                        control = list(abstol = 1e-3, reltol = 1e-5)))
    ret <- optimized$par
  }
  ret
}

# ------------------------------------------------------------------------------
# Functions for creating the metadata
# ------------------------------------------------------------------------------
# create metadata table
# @param dataset the clean analysis dataset
# @param opts the global options
create_metadata <- function(dataset, opts) {
    # set up the names
    all_nms <- names(dataset)
    id_nms <- all_nms[1:2]
    geog_idx <- grepl("geog", all_nms)
    outcome_nms <- all_nms[3:(which(geog_idx)[1] - 1)]
    geog_nms <- all_nms[geog_idx]
    subtype_idx <- grepl("subtype", all_nms)
    subtype_nms <- all_nms[subtype_idx]
    aa_idx <- grepl("hxb2", all_nms)
    aa_var_nms <- all_nms[aa_idx]
    geom_nms <- all_nms[grepl("length", all_nms) | grepl("num", all_nms)]
    # set up the tibble
    metadata_tib <- tibble(Variable = all_nms, Name = NA, Description = NA)
    # make nice names
    nice_ids <- str_to_title(gsub(".", " ", id_nms, fixed = TRUE))
    nice_outcomes <- make_nice_outcomes(outcome_nms)
    nice_geog <- str_to_title(gsub(".", " ", geog_nms, fixed = TRUE))
    nice_subtype <- str_to_title(gsub(".", " ", subtype_nms, fixed = TRUE))
    nice_aas <- str_to_title(gsub(".", " ", vimp_nice_ind_names(aa_var_nms), fixed = TRUE))
    nice_geom <- str_to_title(gsub(".", " ", geom_nms, fixed = TRUE))
    metadata_tib$Name <- c(nice_ids, nice_outcomes, nice_geog, nice_subtype, nice_aas, nice_geom)
    # make descriptions
    id_descriptions <- apply(matrix(id_nms), 1, describe_id_var)
    outcome_descriptions <- apply(matrix(outcome_nms), 1, describe_outcome_var, opts = opts)
    geog_descriptions <- apply(matrix(geog_nms), 1, describe_geog_var)
    subtype_descriptions <- apply(matrix(subtype_nms), 1, describe_subtype_var)
    aa_descriptions <- apply(matrix(aa_var_nms), 1, describe_aa_var)
    geom_descriptions <- apply(matrix(geom_nms), 1, describe_geom_var)
    metadata_tib$Description <- c(id_descriptions, outcome_descriptions, geog_descriptions, subtype_descriptions, aa_descriptions, geom_descriptions)
    return(metadata_tib)
}
make_nice_outcomes <- function(outcome_nms) {
    gsub("multsens", "Multiple Sensitivity", gsub("estsens", "Estimated Sensitivity", gsub("sens", "Sensitivity", gsub("iip", "IIP", gsub("ic80", "IC$_{80}$", gsub("ic50", "IC$_{50}$", outcome_nms))))))
}
# describe id variables for the metadata
describe_id_var <- function(var) {
    if (grepl("id", var)) {
        descr <- "Subject identifier"
        if (grepl("lanl", var)) {
            descr <- paste0(descr, " from LANL viral sequences database [@yoon2015catnap].")
        } else if (grepl("catnap", var)) {
            descr <- paste0(descr, " from CATNAP viral sequence database [@yoon2015catnap].")
        }
    } else {
        descr <- NA
    }
    return(descr)
}
# describe outcomes for the metadata
describe_outcome_var <- function(var, opts) {
    predicted_text_ic50 <- get_combo_text(opts, ic50_pres = TRUE, ic80_pres = FALSE)
    predicted_text_ic80 <- get_combo_text(opts, ic50_pres = FALSE, ic80_pres = TRUE)
    binary_outcomes_txt <- get_binary_outcome_text(opts)
    if (grepl("ic50", var)) {
        descr <- paste0("Outcome variable: $\\log_{10}$ IC$_{50}$ (50% inhibitory concentration)", ifelse(length(opts$nab) == 1, ".", predicted_text_ic50))
    } else if (grepl("ic80", var)) {
        descr <- paste0("Outcome variable: $\\log_{10}$ IC$_{80}$ (80% inhibitory concentration)", ifelse(length(opts$nab) == 1, ".", predicted_text_ic80))
    } else if (grepl("iip", var)) {
        descr <- get_iip_text(opts)
    } else if (var == "sens" | var == "estsens") {
        descr <- paste0("Outcome variable: ", ifelse(length(opts$nab) == 1, "", "estimated "), "sensitivity. Defined as the binary indicator that ", ifelse(length(opts$nab) == 1, "", "estimated "), binary_outcomes_txt, ". Note that in the dataset, 1 denotes sensitive (i.e., ", ifelse(length(opts$nab) == 1, "", "estimated "), binary_outcomes_txt, ") while 0 denotes resistant")
    } else if (var == "multsens") {
        descr <- paste0("Outcome variable: multiple sensitivity. Defined as the binary indicator of having measured ", binary_outcomes_txt, " for at least ", min(c(length(opts$nab), opts$multsens_nab)) ," bnAbs. note that in the dataset, 1 denotes multiple sensitivity (i.e., measured ", binary_outcomes_txt, " for $\\ge$ ", min(c(length(opts$nab), opts$multsens_nab)) ," bnAbs).")
    }
    return(descr)
}
# geographic region descriptions
describe_geog_var <- function(var) {
    split_str <- unlist(strsplit(var, ".", fixed = TRUE))
    descr <- paste0("Binary indicator variable describing the geographic region of origin for the given virus. 1 denotes that the virus is from ", paste(split_str[6:length(split_str)], collapse = " "), ".")
    return(descr)
}
# subtype descriptions
describe_subtype_var <- function(var) {
    split_str <- unlist(strsplit(var, ".", fixed = TRUE))
    descr <- paste0("Binary indicator variable denoting the subtype of the virus. 1 denotes that the virus is subtype ", split_str[length(split_str)], ".")
    return(descr)
}
# AA descriptions
describe_aa_var <- function(var) {
    split_str <- unlist(strsplit(var, ".", fixed = TRUE))
    descr <- paste0("Binary indicator (1 denotes present) of residue containing ", paste(split_str[3:(length(split_str) - 1)], collapse = " "), " at HXB2-referenced site ", split_str[2], ".")
    return(descr)
}
# geometry descriptions
describe_geom_var <- function(var) {
    split_str <- unlist(strsplit(var, ".", fixed = TRUE))
    if (grepl("length", var)) {
        descr <- paste0("Length of ", split_str[length(split_str)])
    } else if (grepl("num", var)) {
        descr <- paste0("Number of ", ifelse(grepl("cysteine", var), "cysteines", "sequons"), " in ", split_str[length(split_str)])
    } else {
        # do nothing
    }
    return(descr)
}
# describe the number of complete sequences for each outcome
get_complete_data_description <- function(opts, ncomplete_ic50, ncomplete_ic80, ncomplete_ic50ic80, one_nab, est_fillin) {
    num_obs <- ""
    outcome_txt <- ""
    print_ic50 <- "ic50" %in% opts$outcomes | "sens1" %in% opts$outcomes | "sens2" %in% opts$outcomes | "iip" %in% opts$outcomes
    print_ic80 <- "ic80" %in% opts$outcomes
    print_ic5080 <- "iip" %in% opts$outcomes
    if (print_ic50) {
        num_obs <- paste0(num_obs, ncomplete_ic50)
        outcome_txt <- paste0(outcome_txt, ifelse(one_nab, "", est_fillin), "IC$_{50}$")
    }
    if (print_ic80 & print_ic5080) {
        num_obs <- paste0(num_obs, ifelse(!print_ic50, "", ", "), ncomplete_ic80)
        outcome_txt <- paste0(outcome_txt, ifelse(!print_ic50, "", ", "), ifelse(one_nab, "", est_fillin), "IC$_{80}$, ", " and both ", ifelse(one_nab, "", est_fillin), "IC$_{50}$ and ", ifelse(one_nab, "", est_fillin), "IC$_{80}$, respectively,")
    } else if (print_ic80) {
        num_obs <- paste0(num_obs, ifelse(!print_ic50, "", " and "), ncomplete_ic80)
        outcome_txt <- paste0(outcome_txt, ifelse(!print_ic50, "", " and "), ifelse(one_nab, "", est_fillin), "IC$_{80}$, respectively,")
    } else {
        # do nothing
    }
    paste0(num_obs, " sequences had complete data for ", outcome_txt, " and were used in the analysis")
}

# ------------------------------------------------------------------------------
# Functions for returning requested objects after completion of SLAPNAP
# ------------------------------------------------------------------------------
# get analysis dataset name (only a problem if you've reused a mount directory)
get_analysis_dataset_name <- function(all_nms, opts) {
    if (length(all_nms) > 1) {
        nms_with_requested_nabs <- all_nms[grepl(paste(opts$nab, collapse = "_"), all_nms)]
        nms_with_only_requested_nabs <- nms_with_requested_nabs[unlist(lapply(strsplit(nms_with_requested_nabs, "_", fixed = TRUE), function(x) length(x) == 2 + length(opts$nab)))]
        date_nab_only <- gsub(".csv", "", gsub("slapnap_", "", nms_with_requested_nabs))
        date_only <- gsub(paste0(paste(opts$nab, collapse = "_"), "_"), "", date_nab_only)
        current_date <- format(as.Date(Sys.getenv('current_date'), "%d%b%Y"), "%d%b%Y")
        closest_date <- which.min(as.Date(current_date, "%d%b%Y") - as.Date(date_only, "%d%b%Y"))
        nm <- nms_with_only_requested_nabs[closest_date]
    } else {
        nm <- all_nms
    }
    return(nm)
}
# get outcome names
get_outcome_names <- function(opts) {
    outcome_names <- c(
        switch("ic50" %in% opts$outcomes, "log10.pc.ic50", NULL),
        switch("ic80" %in% opts$outcomes, "log10.pc.ic80", NULL),
        switch("iip" %in% opts$outcomes, "iip", NULL),
        switch("sens1" %in% opts$outcomes, "dichotomous.1", NULL),
        switch("sens2" %in% opts$outcomes, "dichotomous.2", NULL)
    )
    return(outcome_names)
}
# clean analysis dataset, by changing column names to match how we refer to them
# in the report
# @param input_data the dataset that we use internally
# @param opts the options
# @return a cleaned dataset, with outcome columns matching the report
clean_analysis_dataset <- function(input_data, opts) {
    # drop nab_ columns
    minus_nab_cols <- input_data[, !grepl("nab_", names(input_data))]
    # rename ic50, ic80 cols
    minus_nab_nms <- names(minus_nab_cols)
    correct_ic50_ic80_nms <- gsub("log10.pc.", "", minus_nab_nms)
    # rename sens/estsens/multsens
    if (length(opts$nab) == 1) {
        corrected_nms <- gsub("dichotomous.1", "sens", correct_ic50_ic80_nms)
    } else {
        corrected_nms <- gsub("dichotomous.2", "multsens", gsub("dichotomous.1", "estsens", correct_ic50_ic80_nms))
    }
    corrected_cols <- setNames(minus_nab_cols, corrected_nms)
    # drop unneccessary cols
    analysis_dataset <- corrected_cols[, !(grepl("pc.", corrected_nms) | grepl("dichot", corrected_nms))]
    return(analysis_dataset)
}
# get full learner fit names
get_learner_fit_names <- function(all_fit_nms, opts) {
    fit_nms <- all_fit_nms[grepl("learner_", all_fit_nms)]
    fit_only_outcomes <- gsub(".rds", "", gsub("cv", "", gsub("learner_", "", fit_nms)))
    outcome_names <- get_outcome_names(opts)
    these_outcome_fit_nms <- fit_nms[!is.na(pmatch(fit_only_outcomes, outcome_names, duplicates.ok = TRUE))]
    return(these_outcome_fit_nms)
}

# get vimp object names
get_vimp_object_names <- function(all_fit_nms, opts) {
    vimp_nms <- all_fit_nms[grepl("_vimp", all_fit_nms, fixed = TRUE)]
    vimp_only_outcome <- gsub(".rds", "", gsub("_cv", "", gsub("_vimp", "", vimp_nms)))
    outcome_names <- get_outcome_names(opts)
    these_outcome_vimp_nms <- vimp_nms[!is.na(pmatch(vimp_only_outcome, outcome_names, duplicates.ok = TRUE))]
    return(these_outcome_vimp_nms)
}
