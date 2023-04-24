#' Creates folds for cross validation.
#' 
#' From package caret.
#' 
#' @param y The labels.
#' @param k The number of folds.
#' @param list Logical; should the result be a list or not.
#' @param returnTrain Logical; return training sets or not (test sets).
creat_folds <- function (y, k = 10L, list = TRUE, returnTrain = FALSE) 
{
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i] %/% k
      
      if (min_reps > 0L) {
        spares <- numInClass[i] %% k
        seqVector <- rep(1:k, min_reps)
        
        if (spares > 0L) 
          seqVector <- c(seqVector, sample(1:k, spares))
        
        foldVector[which(y == names(numInClass)[i])] <- 
          sample(seqVector)
      }
      else {
        foldVector[which(y == names(numInClass)[i])] <- 
          sample(1:k, size = numInClass[i])
      }
    }
  }
  else 
    foldVector <- seq_along(y)
  
  if (list) {
    out <- split(seq_along(y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq_along(out))), sep = "")
    
    if (returnTrain) 
      out <- lapply(out, function(data, y) y[-data], y = seq_along(y))
  }
  else 
    out <- foldVector
  
  out
}


#' Calculates cross-validation errors.
#' 
#' Optimizes the regularization cost for svm.
#' 
#' @param train Training set.
#' @param test Test set.
#' @param costs A vector of costs.
#' @param ... Additional arguments for svm.
cv_svm <- function (train, test, costs = c(10E-2, 10E-1, 1, 5, 50), ...)
{
  len  <- length(costs)
  tabs <- errs <- vector("list", len)
  
  for (i in 1:len) {
    m <- e1071::svm(y. ~ ., data = train, cost = costs[i], ...)
    tab_i <- tabs[[i]] <- table(pred = predict(m, test), true = test[["y."]])
    tot <- tab_i[1, 1] + tab_i[1, 2] + tab_i[2, 1] + tab_i[2, 2]
    errs[[i]] <- (tab_i[1, 1] + tab_i[1, 2])/tot
  }

  list(tab = tabs, err = errs)
}


#' Percolator
#'
#' @param df A data frame of \code{psmC.txt}.
#' @param prob_cos Probability cot-offs (as a function of pep_len).
#' @param fct_score The factor in converting probability p-values to scores. The
#'   value is always 10.
#' @param k The k-folds for cross validation.
#' @param cross_valid Logical; to perform cross validations or not.
#' @param costs The costs for cross validations.
#' @param def_cost The default cost.
#' @param svm_tol Tolerance in FDR.
#' @param svm_iters The number of iterations.
#' @inheritParams matchMS
perco_svm <- function (prob_cos = NULL, out_path = NULL, df = NULL, 
                       target_fdr = .01, fdr_type = "protein", 
                       min_len = 7L, max_len = 40L, max_pepscores_co = 50, 
                       min_pepscores_co = 0, enzyme = "trypsin_p", 
                       fdr_group = "base", nes_fdr_group = "base", 
                       fct_score = 10, k = 10, cross_valid = FALSE, 
                       costs = c(.1, .3, 1, 3, 10), def_cost = 1L, 
                       svm_kernel = "radial",
                       svm_feats  = c("pep_score", "pep_ret_range", 
                                      "pep_delta", "pep_n_ms2", 
                                      "pep_expect", # "pep_len", 
                                      "pep_exp_mz", "pep_exp_mr", 
                                      "pep_tot_int", # "pep_mod_group", 
                                      "pep_n_matches2", 
                                      "pep_ms2_deltas_mean"), 
                       svm_iters = 10L, svm_tol = 1E-4, ...)
{
  if (!all(costs > 0))
    costs <- c(.1, .3, 1, 3, 10)
  
  if (def_cost <= 0)
    def_cost <- 1L
  
  if (svm_iters <= 0)
    svm_iters <- 10L
  
  if (svm_tol <= 0)
    svm_tol <- 1E-4
  
  fileC <- file.path(out_path, "psmC.txt")
  fileP <- file.path(out_path, "temp", "prob_cos.rds")
  
  # --- preparation 
  if (is.null(df)) {
    if (is.null(out_path)) 
      stop("Argument \"out_path\" cannot be NULL.")
    
    if (!file.exists(fileC))
      stop("File not found: ", fileC)
    
    df <- readr::read_tsv(fileC)
  }
  
  if (is.null(prob_cos)) {
    if (!file.exists(fileP))
      stop("File not found: ", fileP)
    
    prob_cos <- qs::qread(fileP)
  }
  prob_cos0 <- prob_cos

  if (FALSE) {
    if (nrow(df) <= 500L) {
      message("No SVM post-processing with fewer than 500 observations.")
      return(df)
    }
  }

  if (!"pep_delta" %in% names(df)) 
    df$pep_delta <- df$pep_exp_mr - df$pep_calc_mr
  
  cnms <- names(df)
  
  if (!"pep_issig" %in% cnms) {
    warning("No SVM post-processing without data column \"pep_issig\".")
    return(df)
  }
  
  if (!all(c("raw_file", "pep_scan_num") %in% cnms)) {
    warning("Require columns \"raw_file\" and \"pep_scan_num\".")
    return(df)
  }
  
  if (!"pep_score" %in% cnms) {
    warning("Column \"pep_score\" not found.")
    return(df)
  }
  
  # information already used in getting initial `prob_cos`
  if ("pep_len" %in% svm_feats) {
    warning("Excluded feature \"pep_len\" (information already used).")
    svm_feats <- svm_feats[svm_feats != "pep_len"]
  }
  
  if (FALSE) {
    if ("pep_z_expect" %in% svm_feats && !"pep_z_expect" %in% cnms) {
      df <- df %>% 
        dplyr::left_join(calc_z_pepfdr(out_path = out_path), by = "pep_exp_z") %>% 
        dplyr::mutate(pep_z_expect = 10^((pep_z_prob_co - pep_score)/10) * target_fdr)
    }
  }

  rm(list = c("cnms"))
  
  
  # --- initialization
  # metric for selecting high-quality training PSMs
  if (!"pep_expect" %in% svm_feats)
    svm_feats <- c("pep_expect", svm_feats)
  
  if (!"pep_score" %in% svm_feats)
    svm_feats <- c("pep_score", svm_feats)
  
  if (!"pep_expect" %in% names(df))
    df <- dplyr::mutate(df, pep_expect = 10^((pep_score_co - pep_score)/10) * target_fdr)

  if (!"pep_prob" %in% names(df))
    df <- dplyr::mutate(df, pep_prob = 10^(-pep_score/fct_score))
  
  # note: `df` being altered
  if ("pep_exp_z" %in% names(df)) 
    df[["pep_exp_z"]] <- as.integer(factor(df[["pep_exp_z"]]))
  
  td <- prep_pepfdr_td(df, 
                       out_path = out_path, 
                       enzyme = enzyme, 
                       nes_fdr_group = nes_fdr_group, 
                       fdr_group = fdr_group)
  td <- keep_pepfdr_best(td, cols = c("pep_scan_num", "raw_file"))
  td[["y."]] <- as.factor(td[["pep_issig"]])
  
  oks <- svm_feats %in% names(df)
  
  if (!all(oks)) {
    warning("SVM features not found: ", svm_feats[!oks])
    svm_feats <- svm_feats[oks]
  }
  
  oks <- unlist(lapply(df[svm_feats], is.numeric))
  
  if (!all(oks)) {
    warning("Non-numeric features excluded: ", svm_feats[!oks])
    svm_feats <- svm_feats[oks]
  }
  
  for (pf in svm_feats)
    td[[paste0(pf, ".")]] <- td[[pf]]
  
  rm(list = c("pf", "oks"))
  
  svm_feats <- paste0(svm_feats, ".")
  
  if ("pep_expect." %in% svm_feats)
    td[["pep_expect."]] <- -log10(td[["pep_expect."]])
  
  if (FALSE) {
    if ("pep_z_expect." %in% svm_feats)
      td[["pep_z_expect."]] <- -log10(td[["pep_z_expect."]])
  }

  nas <- lapply(svm_feats, function (x) is.na(td[[x]]))
  nas <- Reduce(`|`, nas)
  td0 <- td[nas, ]
  td1 <- td[!nas, ]
  
  rows <- td1[["pep_isdecoy"]]
  ta <- td1[!rows, c("y.", svm_feats), drop = FALSE]
  de <- td1[ rows, c("y.", svm_feats), drop = FALSE]
  rm(list = c("rows"))
  
  ## (2) train and test for each of target and decoy sets
  oks <- ta[["pep_expect."]] >= median(ta[["pep_expect."]], na.rm = TRUE)
  ttrain <- ta[oks, ]
  ttest  <- ta[!oks, ]
  
  rows <- sample(c(TRUE, FALSE), nrow(de), replace = TRUE)
  dtrain <- de[rows, ]
  dtest  <- de[rows, ]
  train  <- dplyr::bind_rows(ttrain, dtrain)
  rm(list = c("dtrain", "rows", "oks"))
  
  if (cross_valid) {
    mses  <- vector("numeric", length(costs))
    
    folds <- creat_folds(ta[["y."]], k = k)
    tests <- trains <- vector("list", k)
    
    for (i in seq_len(k)) {
      tests[[i]]  <- ta[folds[[i]], ]
      trains[[i]] <- ta[-folds[[i]], ]
    }
    
    n_cores <- min(mzion:::detect_cores(16L), k)
    cl  <- parallel::makeCluster(getOption("cl.cores", n_cores))
    cvs <- parallel::clusterMap(cl, cv_svm, trains, tests, 
                                MoreArgs = list(costs = costs), 
                                SIMPLIFY = FALSE, USE.NAMES = FALSE)
    parallel::stopCluster(cl)
    errs <- lapply(cvs, `[[`, "err")
    
    for (i in seq_along(costs))
      mses[[i]] <- mean(unlist(lapply(errs, `[[`, i)))
    
    best_co <- costs[which.min(mses)]
  }
  else 
    best_co <- def_cost
  
  message("Regularization cost: ", best_co)
  
  fit <- tryCatch(
    e1071::svm(y. ~ ., data = train, kernel = svm_kernel, 
               cost = best_co, ...), 
    error = function (e) NULL)
  
  if (is.null(fit))
    return(prob_cos)

  pred <- tryCatch(
    as.logical(predict(fit, td1[, svm_feats])), 
    error = function (e) NULL)
  
  if (is.null(pred))
    return(prob_cos)

  ###
  # tdr_bare - target-decoy rate (TDR) based on one feature (pep_len)
  # tdr_svm  - TDR based on multiple features; initially tdr_svm < tdr_bare ->
  #            ^target_fdr -> ^tdr_svm ... tdr_svm == tdr_bare
  ###
  
  tdr_bare <- sum(td1[td1$pep_isdecoy, "pep_issig"])/sum(td1[!td1$pep_isdecoy, "pep_issig"])
  td1[, "pep_issig"] <- as.logical(pred)
  tdr_svm <- sum(td1[td1$pep_isdecoy, "pep_issig"])/sum(td1[!td1$pep_isdecoy, "pep_issig"])
  delta <- tdr_svm - tdr_bare
  rm(list = c("tdr_svm"))
  
  # if ((delta > 0) || (abs(delta) <= svm_tol)) return(prob_cos)

  
  # --- iteration 
  fdr0 <- target_fdr
  fdr1 <- target_fdr * 5
  fdrm <- (fdr0 + fdr1)/2
  step <- fdr1 - fdr0
  
  while((svm_iters > 0L) && (abs(delta) > svm_tol)) {
    prob_cos <- calc_pepfdr(target_fdr = fdrm, 
                            fdr_type = fdr_type, 
                            min_len = min_len, 
                            max_len = max_len, 
                            max_pepscores_co = max_pepscores_co, 
                            min_pepscores_co = min_pepscores_co, 
                            enzyme = enzyme, 
                            fdr_group = fdr_group, 
                            nes_fdr_group = nes_fdr_group, 
                            out_path = out_path)
    df <- post_pepfdr(prob_cos, out_path) # also updates pepfdr.rds
    
    if (!"pep_delta" %in% names(df)) 
      df$pep_delta <- df$pep_exp_mr - df$pep_calc_mr
    
    if (!"pep_expect" %in% names(df))
      df <- dplyr::mutate(df, pep_expect = 10^((pep_score_co - pep_score)/10) * fdrm)
    
    if (!"pep_prob" %in% names(df))
      df <- dplyr::mutate(df, pep_prob = 10^(-pep_score/fct_score))
    
    td <- prep_pepfdr_td(df, 
                         out_path = out_path, 
                         enzyme = enzyme, 
                         nes_fdr_group = nes_fdr_group, 
                         fdr_group = fdr_group)
    td <- keep_pepfdr_best(td, cols = c("pep_scan_num", "raw_file"))
    td[["y."]] <- as.factor(td[["pep_issig"]])
    
    svm_feats <- gsub("\\.", "", svm_feats)
    
    for (pf in svm_feats)
      td[[paste0(pf, ".")]] <- td[[pf]]
    
    rm(list = c("pf"))
    svm_feats <- paste0(svm_feats, ".")
    
    if ("pep_expect." %in% svm_feats)
      td[["pep_expect."]] <- -log10(td[["pep_expect."]])

    nas <- lapply(svm_feats, function (x) is.na(td[[x]]))
    nas <- Reduce(`|`, nas)
    
    if (length(nas)) {
      td0 <- td[nas, ]
      td1 <- td[!nas, ]
    }
    else {
      td0 <- NULL
      td1 <- td
    }
    
    rows <- td1[["pep_isdecoy"]]
    ta <- td1[!rows, c("y.", svm_feats), drop = FALSE]
    de <- td1[ rows, c("y.", svm_feats), drop = FALSE]
    rm(list = c("rows"))
    
    ## (2) train and test for each of target and decoy sets
    oks <- ta[["pep_expect."]] >= median(ta[["pep_expect."]], na.rm = TRUE)
    ttrain <- ta[oks, ]
    ttest  <- ta[!oks, ]
    
    rows <- sample(c(TRUE, FALSE), nrow(de), replace = TRUE)
    dtrain <- de[rows, ]
    dtest  <- de[rows, ]
    train  <- dplyr::bind_rows(ttrain, dtrain)
    rm(list = c("dtrain", "rows", "oks"))
    
    fit <- tryCatch(e1071::svm(y. ~ ., data = train, kernel = svm_kernel, 
                               cost = best_co, ...), 
                    error = function (e) NULL)
    
    if (is.null(fit))
      return(if (all(prob_cos[, 2] >= prob_cos0[, 2])) prob_cos else prob_cos0)
    
    pred <- tryCatch(
      as.logical(predict(fit, td1[, svm_feats])), 
      error = function (e) NULL)
    
    if (is.null(pred))
      return(if (all(prob_cos[, 2] >= prob_cos0[, 2])) prob_cos else prob_cos0)

    td1[, "pep_issig"] <- as.logical(pred)
    tdr_svm <- sum(td1[td1$pep_isdecoy, "pep_issig"])/sum(td1[!td1$pep_isdecoy, "pep_issig"])
    delta <- tdr_svm - tdr_bare
    
    if (delta < 0) {
      # next grid
      fdr0 <- fdrm
      fdrm <- (fdr0 + fdr1)/2
    }
    else {
      # left-half grid
      fdr1 <- fdrm
      fdrm <- (fdr0 + fdr1)/2
      step <- fdrm - fdr0
    }
    
    message("Range of adjusted FDR: ", fdr0, " : ", fdr1)
    svm_iters <- svm_iters - 1L
  }
  
  if (all(prob_cos[, 2] >= prob_cos0[, 2])) prob_cos else prob_cos0
}


