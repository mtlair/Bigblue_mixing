# ==============================================================================
# morris_shim.R - Offline fallback providing morris() / tell() when the CRAN
# package 'sensitivity' cannot be installed (no network access to CRAN).
#
# Code extracted verbatim from the 'sensitivity' package v1.31.0
# (Gilles Pujol, Bertrand Iooss, Alexandre Janon et al., licensed GPL-2),
# files R/base.R (tell generic), R/morris.R, R/morris_oat.R, R/morris_sfd.R.
# Only the Morris OAT screening path used by this project is included.
# ==============================================================================

tell <- function(x, y = NULL, ...) UseMethod("tell")


ind.rep <- function(i, p) {
# indices of the points of the ith trajectory in the DoE
  (1 : (p + 1)) + (i - 1) * (p + 1)
}


morris <- function(model = NULL, factors, r, design, binf = 0, bsup = 1, scale = TRUE, ...) {
  
  # argument checking: factor number and names
  if (is.character(factors)) {
    X.labels <- factors
    p <- length(X.labels)
  } else if (is.numeric(factors)) {
    p <- factors
    X.labels <- paste("X", 1 : p, sep="")
  } else {
    stop("invalid argument \'factors\', waiting for a scalar (number) or a character string vector (names)")
  }
  
  # argument checking: number of repetitions
  if (length(r) == 1) {
    r.max <- r
  } else {
    r.max <- r[2]
    r <- r[1]
  }
  
  # argument checking: design parameters
  if (! "type" %in% names(design)) {
    design$type <- "oat"
    warning("argument \'design$type\' not found, set at \'oat\'")
  }
  if (design$type == "oat") {
    # one-at-a-time design
    if (! "levels" %in% names(design)) {
      stop("argument \'design$levels\' not found")
    }
    nl <- design$levels
    if (length(nl) == 1) nl <- rep(nl, p)
    if ("grid.jump" %in% names(design)) {
      jump <- design$grid.jump
      if (any(round(jump, 0) != jump)) stop("grid.jump must be integer")
      if (length(jump) == 1) jump <- rep(jump, p)
    } else {
      jump <- rep(1, p)
      warning("argument \'design$grid.jump\' not found, set at 1")
    }
  } else if (design$type == "simplex") {
    # simplex-based design
    if (! "scale.factor" %in% names(design)) {
      stop("argument \'design$scale.factor\' not found")
    }
    h <- design$scale.factor
  } else {
    stop("invalid argument design$type, waiting for \"oat\" or \"simplex\"")
  }
  
  # argument checking: domain boundaries
  if (length(binf) == 1) binf <- rep(binf, p)
  if (length(bsup) == 1) bsup <- rep(bsup, p)
  
  # generation of the initial design
  if (design$type == "oat") {
    X <- random.oat(p, r.max, binf, bsup, nl, jump)
  } else if (design$type == "simplex") {
    X <- random.simplexes(p, r.max, binf, bsup, h)
  }
  
  # duplicated repetitions are removed
  X.unique <- array(t(X), dim = c(p, p + 1, r.max))
  X.unique <- unique(X.unique, MARGIN = 3)
  X <- matrix(X.unique, ncol = p, byrow = TRUE)
  colnames(X) <- X.labels
  r.unique <- nrow(X) / (p + 1)
  if (r.unique < r.max) {
    warning(paste("keeping", r.unique, "repetitions out of", r.max))
  }
  r.max <- r.unique
  
  # optimization of the design
  if (r < r.max) {
    ind <- morris.maximin(X, r)
    X <- X[sapply(ind, function(i) ind.rep(i, p)),]
  }
  
  # object of class "morris"
  x <- list(model = model, factors = factors, r = r, design = design,
            binf = binf, bsup = bsup, scale = scale, X = X, call =
              match.call())
  class(x) <- "morris"
  
  # computing the response if the model is given
  if (!is.null(x$model)) {
    response(x, other_types_allowed = TRUE, ...)
    tell(x)
  }
  
  return(x)
}
tell.morris <- function(x, y = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else if (is.null(x$y)) {
    stop("y not found")
  }
  
  X <- x$X
  y <- x$y
  
  if (x$scale) {
    #X <- scale(X)
    #y <- as.numeric(scale(y))
    Binf <- matrix(x$binf, nrow = nrow(X), ncol = length(x$binf), byrow = TRUE)
    Bsup <- matrix(x$bsup, nrow = nrow(X), ncol = length(x$bsup), byrow = TRUE)
    X <- (X - Binf) / (Bsup - Binf) 
  }
  
  if (x$design$type == "oat") {
    x$ee <- ee.oat(X, y)
  } else if (x$design$type == "simplex") {
    x$ee <- ee.simplex(X, y)
  }
  
  assign(id, x, parent.frame())
}


random.oat <- function(p, r, binf = rep(0, p), bsup = rep(0, p), nl, design.step) {
  # orientation matrix B
  B <- matrix(-1, nrow = p + 1, ncol = p)
  B[lower.tri(B)] <- 1
  # grid step
  delta <- design.step / (nl - 1)
  X <- matrix(nrow = r * (p + 1), ncol = p)
  for (j in 1 : r) {
    # directions matrix D
    D <- diag(sample(c(-1, 1), size = p, replace = TRUE), nrow = p)
    # permutation matrix P
    perm <- sample(p)
    P <- matrix(0, nrow = p, ncol = p)
    for (i in 1 : p) {
      P[i, perm[i]] <- 1
    }
    # starting point
    x.base <- matrix(nrow = p + 1, ncol = p)
    for (i in 1 : p) {
      x.base[,i] <- ((sample(nl[i] - design.step[i], size = 1) - 1) / (nl[i] - 1))
    }
    X[ind.rep(j,p),] <- 0.5 * (B %*% P %*% D + 1) %*% 
      diag(delta, nrow = p) + x.base
  }
  for (i in 1 : p) {
    X[,i] <- X[,i] * (bsup[i] - binf[i]) + binf[i]
  }
  return(X)
}

ee.oat <- function(X, y) {
  # compute the elementary effects for a OAT design
  p <- ncol(X)
  r <- nrow(X) / (p + 1)
  
#  if(is(y,"numeric")){
  if(inherits(y, "numeric")){
    one_i_vector <- function(i){
      j <- ind.rep(i, p)
      j1 <- j[1 : p]
      j2 <- j[2 : (p + 1)]
      # return((y[j2] - y[j1]) / rowSums(X[j2,] - X[j1,]))
      return(solve(X[j2,] - X[j1,], y[j2] - y[j1]))
    }
    ee <- vapply(1:r, one_i_vector, FUN.VALUE = numeric(p))
    ee <- t(ee)
    # "ee" is now a (r times p)-matrix.
#  } else if(is(y,"matrix")){
    ee <- vapply(1:r, one_i_vector, FUN.VALUE = numeric(p))
    ee <- t(ee)
    colnames(ee) <- colnames(X)
  } else {
    stop("morris_shim: only numeric y supported")
  }
  return(ee)
}
hausdorff.distance <- function(x, set1, set2) {
# Hausdorff distance function
# x: matrix of points.
# set1: indices of points (in x) of the first group.
# set2: indices of points (in x) of the second group.
# returns: the Haussdorf distance between the two sets of points.
  n1 <- length(set1)
  n2 <- length(set2)
  d <- matrix(nrow = n1, ncol = n2)
  for (i1 in 1 : n1) {
    for (i2 in 1 : n2) {
      d[i1,i2] <- sqrt(sum((x[set1[i1],] - x[set2[i2],])^2))
    }
  }
  return(max(mean(apply(d, 1, min)), mean(apply(d, 2, min))))
}
kennard.stone <- function(dist.matrix, n) {
# Kennard & Stone algorithm (1969).
# dist.matrix: distance matrix (N * N) (cf help(dist)).
# n: number of points to keep (n < N).
# returns: the indices of the n chosen points.
  out <- numeric(n)
  out[1] <- 1
  for (i in 2 : n) {
    tmp <- dist.matrix[out, -out, drop = FALSE]
    # Remark: drop = FALSE since 'out' is of length 1 at the first
    # iteration, cf help(Extract) for the meaning of 'drop'
    out[i] <- (1 : nrow(dist.matrix))[-out][which.max(apply(tmp, 2, min))]
  }
  return(out)
}


morris.maximin <- function(x, r) {
# Select r repetitions (out of the R ones of the "morris" design x)
# that are "space-filling".
# returns: the indices (in 1:R) of the r selected repetitions.
  p <- ncol(x)
  R <- nrow(x) / (p + 1)
  d <- matrix(0, nrow = R, ncol = R)
  if (requireNamespace("pracma", quietly = TRUE)) {
    for (i in 1 : (R - 1)) {
      for (j in (i + 1) : R) {
        d[i,j] <- d[j,i] <- hausdorff.distance2(x, ind.rep(i, p), ind.rep(j, p))
      }
    }
  } else {
    for (i in 1 : (R - 1)) {
      for (j in (i + 1) : R) {
        d[i,j] <- d[j,i] <- hausdorff.distance(x, ind.rep(i, p), ind.rep(j, p))
      }
    }
  }
  kennard.stone(d, r)
}
