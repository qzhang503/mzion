#' A wrapper of \link[stats]{dist} with the handling of partial argument
#' matches.
#' 
#' Not yet used.
#' 
#' @param ... Arguments for \link[stats]{dist}
my_dist <- function (...) 
{
  dots <- as.list(substitute(...()))
  
  dummies <- c("p")
  
  msgs <- c(
    "`p` in `dist()` not used."
  )
  
  # stopifnot(length(dummies) == length(msgs))
  
  purrr::walk2(dummies, msgs, ~ {
    if (.x %in% names(dots)) {
      warning(.y, call. = FALSE)
      dots[[.x]] <<- NULL
    } 
  })
  
  do.call(stats::dist, dots)
}


#' Cosine similarity.
#' 
#' Not yet used.
#' 
#' @param M An input Matrix.
#' 
#' Vectors are in rows: mtcars[1, ].
#' 
#' @examples 
#' library(mzion)
#' 
#' sim <- mzion:::cos_sim(as.matrix(mtcars[1:3, ]))
#' 
#' # distances (against normalized vectors)
#' as.dist(1 - sim)
cos_sim <- function (M) 
{
  if (!is.matrix(M))
    stop("The input data is not a matrix.")

  L <- sqrt(rowSums(M * M)) # vector lengths; no `mean` subtraction
  Mn <- M / L # normalized M 
  Mn %*% t(Mn) # dot products; vectors in the rows of Mn
}


#' Match MS at no enzymatic specificity.
#' 
#' @param ... Variable arguments.
#' @export
matchMS_NES <- function (...)
{
  args <- as.list(substitute(...()))
  nms  <- names(args)
  oks  <- nms != ""
  args <- args[oks]
  nms  <- nms[oks]
  
  if ("max_miss" %in% nms)
    warning("Argument \"max_miss\" has no effect at no enzymatic specificity.")
  
  if ("enzyme" %in% nms && tolower(args[["enzyme"]]) != "noenzyme") {
    warning("Changed setting to `enzyme = Noenzyme`.")
    args["enzyme"] <- NULL
  }
  
  do.call(matchMS, c(list(enzyme = "noenzyme", 
                          # bypass old matchMS_noenzyme
                          bypass_noenzyme = TRUE, 
                          # call add_protacc at enzyme = "noenzyme"
                          direct_prot_acc = TRUE), 
                     args))
}


#' Re-match MS at no enzymatic specificity.
#' 
#' Bypasses ion matches.
#' 
#' @param ... Variable arguments.
#' @export
rematchMS_NES <- function (...)
{
  args <- as.list(substitute(...()))
  nms  <- names(args)
  oks  <- nms != ""
  args <- args[oks]
  nms  <- nms[oks]
  
  if ("max_miss" %in% nms)
    warning("Argument \"max_miss\" has no effect at no enzymatic specificity.")
  
  if ("enzyme" %in% nms && tolower(args[["enzyme"]]) != "noenzyme") {
    warning("Changed setting to `enzyme = Noenzyme`.")
    args["enzyme"] <- NULL
  }
  
  do.call(matchMS, c(list(enzyme = "noenzyme", 
                          bypass_noenzyme = TRUE, 
                          bypass_pepmasses = TRUE,
                          bypass_bin_ms1 = TRUE,
                          bypass_mgf = TRUE,
                          bypass_ms2match = TRUE), 
                     args))
}


