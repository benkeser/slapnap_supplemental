subset.summary.myCV.SuperLearner <- function (object, obsWeights = NULL, 
                                              method = NULL, opts, 
                                              subset_vec = NULL, ...) {
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
                these_folds <- folds[[ii]][folds[[ii]] %in% which(subset_vec == 1)]
                Risk.SL[ii] <- mean(obsWeights[these_folds] * (Y[these_folds] -
                                                                   SL.predict[these_folds])^2)
                if(is_cvsl){
                    Risk.dSL[ii] <- mean(obsWeights[these_folds] * (Y[these_folds] -
                                                                        discreteSL.predict[these_folds])^2)
                    Risk.library[, ii] <- apply(library.predict[these_folds,
                                                                , drop = FALSE], 2, function(x) mean(obsWeights[these_folds] *
                                                                                                         (Y[these_folds] - x)^2))
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
                these_folds <- folds[[ii]][folds[[ii]] %in% which(subset_vec == 1)]
                sl_auc <- tryCatch(cvAUC::ci.cvAUC(predictions = SL.predict[these_folds],
                                          labels = Y[these_folds], folds = NULL),
                                   error = function(e) list(cvAUC = NA, se = NA, ci = c(NA, NA), confidence = 0.95))
                Risk.SL[ii] <- sl_auc$cvAUC
                se.SL[ii] <- sl_auc$se
                if(is_cvsl){
                    dsl_auc <- tryCatch(cvAUC::ci.cvAUC(predictions = discreteSL.predict[these_folds],
                                               labels = Y[these_folds], folds = NULL),
                                        error = function(e) list(cvAUC = NA, se = NA, ci = c(NA, NA), confidence = 0.95))
                    Risk.dSL[ii] <- dsl_auc$cvAUC
                    se.dSL[ii] <- dsl_auc$se
                    library_auc <- apply(library.predict[these_folds, , drop = FALSE], 2, function(x){
                        tmp <- tryCatch(cvAUC::ci.cvAUC(predictions = x, labels = Y[these_folds], folds = NULL),
                                        error = function(e) list(cvAUC = NA, se = NA, ci = c(NA, NA), confidence = 0.95))
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
                                                      library.names), Ave = c(mean(Risk.SL, na.rm = TRUE), mean(Risk.dSL, na.rm = TRUE),
                                                                              apply(Risk.library, 1, mean, na.rm = TRUE)), se = se, Min = c(min(Risk.SL, na.rm = TRUE),
                                                                                                                              min(Risk.dSL, na.rm = TRUE), apply(Risk.library, 1, min, na.rm = TRUE)), Max = c(max(Risk.SL, na.rm = TRUE),
                                                                                                                                                                                   max(Risk.dSL, na.rm = TRUE), apply(Risk.library, 1, max, na.rm = TRUE)))
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

plot_predicted_prob_boxplots_subset <- function(cv_fit, topRank = 1, opts,
                                                subset_vec = NULL, shape_lab = NULL){
    allAlgos <- subset.summary.myCV.SuperLearner(cv_fit, method = "method.AUC", 
                                                 opts = opts, 
                                                 subset_vec = subset_vec)$Table %>% 
        mutate(Algorithm = relabel_library(as.character(Algorithm), opts))
    colnames(cv_fit[["library.predict"]]) <- relabel_library(colnames(cv_fit[["library.predict"]]), opts)
    # if CV.super learner then top two rows will be super learner and DSL
    if(length(opts$learners) > 1){
        # if super learner show ROC for SL, DSL, and top topRank algorithms
        allCandidates <- allAlgos[-(1:2), ] %>% arrange(-Ave)
        sortedAlgos <- rbind(allAlgos[1:2,], allCandidates[1:topRank,])
        
        predict <- cv_fit[["library.predict"]] %>% as.data.frame() %>%
            bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("cv selector"))) %>%
            bind_cols(cv_fit[["SL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("super learner"))) %>%
            bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y"))) %>%
            bind_cols(tibble::tibble(subset = subset_vec)) %>% 
            gather("algo", "pred", -Y, -subset) %>%
            filter(algo %in% sortedAlgos$Algorithm[1:(2+topRank)])
        cv_folds <- rep(NA, length(cv_fit$Y))
        for(v in seq_along(cv_fit$folds)){
            cv_folds[cv_fit$folds[[v]]] <- v
        }
        predict$cv_folds <- rep(cv_folds, 3)
    }else{
        if("SuperLearner" %in% class(cv_fit)){
            # if here, then just one algo, so just show curve for one algo
            predict <- cv_fit[["Z"]][,1] %>% as.data.frame() %>%
                bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y")))
            predict$algo <- relabel_library(cv_fit$libraryNames[1], opts)
            colnames(predict)[1] <- "pred"
            cv_folds <- rep(NA, length(cv_fit$Y))
            for(v in seq_along(cv_fit$validRows)){
                cv_folds[cv_fit$validRows[[v]]] <- v
            }
            predict$cv_folds <- cv_folds
        }else{
            # if here, then just want to show discrete super learner curve
            predict <- bind_cols(cv_fit[["discreteSL.predict"]] %>% as.data.frame() %>% `colnames<-`(c("cv selector"))) %>%
                bind_cols(cv_fit[["Y"]] %>% as.data.frame() %>% `colnames<-`(c("Y")))
            predict$algo <- "cv selector"
            colnames(predict)[1] <- "pred"
            cv_folds <- rep(NA, length(cv_fit$Y))
            for(v in seq_along(cv_fit$folds)){
                cv_folds[cv_fit$folds[[v]]] <- v
            }
            predict$cv_folds <- cv_folds
        }
    }
    # need to double check this code
    cv_fold_palette <- RColorBrewer::brewer.pal(max(c(cv_folds,3)), "Set1")
    predict %>%
        mutate(Sensitivity = if_else(Y==1, "Sensitive", "Resistant")) %>%
        filter(!is.na(Sensitivity)) %>%
        ggplot(aes(x = Sensitivity, y = pred, color = factor(cv_folds),
                   shape = subset)) +
        facet_grid(. ~ algo) +
        geom_boxplot(outlier.shape = NA) +
        labs(color = "CV fold", shape = shape_lab) +
        geom_point(position=position_jitterdodge(), aes(colour = factor(cv_folds)))+
        ylab(paste0("Predicted Probability of Sensitivity")) + xlab("") +
        theme_bw() + coord_cartesian(ylim=c(0,1)) +
        scale_color_manual(values = cv_fold_palette)
}

