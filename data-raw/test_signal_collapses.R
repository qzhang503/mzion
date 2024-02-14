foo_collapse <- function (xs, ys, zs, ups, unv, lwr = 115L, step = 1e-5)
{
  if (FALSE) {
    xout <- matrix(nrow = 5, ncol = 11)
    
    for (i in 1:5) {
      oks <- ups[[i]]
      xout[i, oks] <- xs[[i]]
    }
    
    # mat <- matrix(rep_len(list(), 5*11),nrow = 5, ncol = 11) 
  }
  
  
  ##
  ps <- mzion:::find_gates(unv)
  lenp <- length(ps)
  
  ps2 <- lapply(ps, `[[`, 2)
  ps2 <- .Internal(unlist(ps2, recursive = FALSE, use.names = FALSE))
  
  for (i in 1:lenp) {
    c12 <- ps[[i]]
    c2 <- c12[2]
    
    # possible to have values in both columns but simply overwrite: 1 <- 2
    # c2 can be 0L, e.g., at three adjacent indexes: c(4, 5), c(6, 0)
    if (!c2)
      next
    
    cols2 <- lapply(ups, function (x) .Internal(which(x == c2)))
    rows <- .Internal(which(lengths(cols2) > 0L))
    cols2 <- .Internal(unlist(cols2[rows], recursive = FALSE, use.names = FALSE))
    
    ui <- ups[rows]
    xi <- xs[rows]
    yi <- ys[rows]
    
    for (j in seq_along(rows)) {
      uj <- ui[[j]]
      cj <- cols2[[j]]
      cj0 <- cj - 1L
      
      # adjacent in bin indexes
      if (uj[[cj]] - uj[[cj0]] == 1L) {
        xj <- xi[[j]]
        xj[[cj0]] <- xj[[cj]]
        xj[[cj]] <- NA_real_ # may store cj for the removals of NA entries
        xs[rows][[j]] <- xj
        
        yj <- yi[[j]]
        yj[[cj0]] <- yj[[cj]]
        yj[[cj]] <- NA_real_
        ys[rows][[j]] <- yj
      }
      else {
        uj[[cj]] <- uj[[cj]] - 1L
        ups[rows][[j]] <- uj
      }
    }
  }
  
  oks <- lapply(xs, function (x) !is.na(x))
  
  for (i in seq_along(xs)) {
    oki <- oks[[i]]
    ups[[i]] <- ups[[i]][oki]
    xs[[i]] <- xs[[i]][oki]
    ys[[i]] <- ys[[i]][oki]
  }
  
  list(xs = xs, ys = ys, zs = zs, ups = ups)
}




foo <- function ()
{
  unv <- c(1, 3, 4, 7:9, 12, 15:16, 19:20)
  
  ups <- list(
    c(1, 2),
    c(2:6, 8:9, 11),
    c(2:4, 7, 9), 
    c(1:2, 5:6, 10),
    c(4:5, 8)
  )
  
  xs <- list(
    c(100, 130.008),
    c(130.008, 130.1, 150.0, 150.001, 150.002, 199.009, 199.1, 220.01),
    c(130.008, 130.1, 150.0, 170, 199.1), 
    c(100, 130.008, 150.001, 150.002, 220.009),
    c(150.0, 150.001, 199.009)
  )
  
  ys <- list(
    c(5, 1),
    c(1, -1, 1, -4, -2, 1, -3, 1),
    c(1, -1, 1, 1, -3), 
    c(5, 1, -4, -2, 1),
    c(1, -4, 1)
  )
  
  zs <- list(
    c(2, 3),
    c(3, 3, 4, 4, 2, 1, 3, 1),
    c(3, 3, 2, 2, 3), 
    c(2, 3, 2, 2, 1),
    c(2, 2, 1)
  )
  
  ans <- foo_collapse(xs, ys, zs, ups, unv)
  
}