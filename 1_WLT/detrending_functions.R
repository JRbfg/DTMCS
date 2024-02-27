# function dctn
dctn <- function(y) {
  y <- matrix(y)
  siz <- pracma::size(y)
  n <- siz[1]
  vec1 <- seq(from = 1, to = n, by = 2)
  vec2 <- seq(from = 2 * floor(n/2), to = 2, by = -2)
  vec <- c(vec1, vec2)
  y <- y[vec, ]
  # y <- reshape(y,n,[])
  y <- y * sqrt(2 * n)
  y <- ifft(y)
  y <- Real(y)
  y[1] <- y[1]/sqrt(2)
  y <- c(y)
  return(y)
}

# function idctn
idctn <- function(y) {
  y <- matrix(y)
  siz <- pracma::size(y)
  n <- siz[1]
  # siz = size(y); n = siz(1); y = reshape(y,n,[]); y =
  # bsxfun(@times,y,w{dim});
  y[1] <- y[1]/sqrt(2)
  y <- c(y)
  y <- ifft(y)
  y <- Real(y * sqrt(2 * n))
  I <- c(1:n) * 0.5 + 0.5
  vec1 <- seq(from = 2, to = length(I), by = 2)
  vec2 <- seq(from = 1, to = length(I) - 1, by = 2)
  I[vec1] <- n - I[vec2] + 1
  y <- y[I]
  y <- c(y)
  return(y)
}

# function searchclosest
searchclosest <- function(x, v) {
  i <- data.frame()
  lo <- 1
  hi <- length(x)
  
  # Phase 1: Binary Search
  while (hi - lo > 1) {
    mid <- round((hi + lo)/2)
    diff <- x[mid] - v
    if (diff == 0) {
      i <- mid
      cv <- v
      closest <- list(i, cv)
      return(closest)
    } else if (diff < 0) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  if (abs(x[lo] - v) < abs(x[hi] - v)) {
    i <- lo
    cv <- x[lo]
    closest <- list(i, cv)
    return(closest)
  } else {
    i <- hi
    cv <- x[hi]
    closest <- list(i, cv)
    return(closest)
  }
}

# function SMOOTHN function for filtering signal
SMOOTHN <- function(y, s){
  k <- 0
  sizy <- pracma::size(y)
  noe <- prod(sizy)  # number of elements
  if(noe < 2){
    z <- y
  }
  
  # Smoothness parameter and weights
  W <- matrix(1, sizy[1], sizy[2])
  # if (k==2){ if (is.null(s) || is_scalar_atomic(s)){ # smoothn(y,s) s
  # <- s # smoothness parameter else { # smoothn(y,W) W <- s} # weight
  # array} elseif (k==3){ # smoothn(y,W,s) W = varargin{2} # weight array
  # s = varargin{3} # smoothness parameter end if ~isequal(size(W),sizy)
  # error('MATLAB:smoothn:SizeMismatch',...  'Arrays for data and weights
  # must have same size.') elseif ~isempty(s) && (~isscalar(s) || s<0)
  # error('MATLAB:smoothn:IncorrectSmoothingParameter',...  'The
  # smoothing parameter must be a scalar >=0') end
  
  MaxIter <- 100
  TolZ <- 0.001
  isinitial <- FALSE
  
  IsFinite <- is.finite(y)
  nof <- table(IsFinite)
  nof <- data.frame(nof)
  nof <- filter(nof, nof$IsFinite == "TRUE")
  nof <- nof$Freq
  W <- W * IsFinite
  
  if(any(W < 0)){
    print("MATLAB:smoothn:NegativeWeights", "Weights must all be >=0")
  }else{
    W <- W/max(W[, ])
  }
  
  # Weighted or missing data?
  isweighted <- any(W[, ] < 1)
  
  # Robust smoothing?
  isrobust <- FALSE  #any(strcmpi(varargin,'robust'));
  
  # Automatic smoothing?
  isauto = is.null(s)
  
  ## Creation of the Lambda tensor
  #---
  # Lambda contains the eingenvalues of the difference matrix used in
  # this penalized least squares process.
  d <- ndims(matrix(y))
  Lambda <- matrix(0, sizy[1], sizy[2])
  for(i in 1:d){
    
    siz0 <- matrix(1, 1, d)
    siz0[i] <- sizy[i]
    help <- cos(pi * (c(1:sizy[i]) - 1)/sizy[i])
    Lambda <- Lambda + help
  }
  
  Lambda <- -2 * (d - Lambda)
  if(isauto == F){
    Gamma <- 1/(1 + s * Lambda^2)
  }
  
  ## Upper and lower bound for the smoothness parameter The average
  ## leverage (h) is by definition in [0 1]. Weak smoothing occurs if h is
  ## close to 1, while over-smoothing appears when h is near 0. Upper and
  ## lower bounds for h are given to avoid under- or over-smoothing. See
  ## equation relating h to the smoothness parameter (Equation #12 in the
  ## referenced CSDA paper).
  N <- sum(sizy != 1)  # tensor rank of the y-array
  hMin <- 1e-06
  hMax <- 0.99
  sMinBnd <- (((1 + sqrt(1 + 8 * hMax^(2/N)))/4/hMax^(2/N))^2 - 1)/16
  sMaxBnd <- (((1 + sqrt(1 + 8 * hMin^(2/N)))/4/hMin^(2/N))^2 - 1)/16
  
  
  ## Initialize before iterating
  #---
  Wtot <- W
  #--- Initial conditions for z
  if(isweighted == TRUE){
    #--- With weighted/missing data
    # An initial guess is provided to ensure faster convergence. For that
    # purpose, a nearest neighbor interpolation followed by a coarse
    # smoothing are performed.
    #---
    if(isinitial == TRUE){
      # an initial guess (z0) has been provided
      z <- z0
    }else{
      z <- InitialGuess(y, IsFinite)
    }
  }else{
    z <- matrix(0, sizy[1], sizy[2])
  }
  #---
  z0 <- z
  y[!IsFinite] <- 0  # arbitrary values for missing y-data
  #---
  tol <- 1
  RobustIterativeProcess <- TRUE
  RobustStep <- 1
  nit <- 0
  #--- Error on p. Smoothness parameter s = 10^p
  errp <- 0.1
  # opt <- optimset('TolX',errp)
  #--- Relaxation factor RF: to speedup convergence
  RF <- 1 + 0.75 * isweighted
  
  ## Main iterative process
  #---
  while(RobustIterativeProcess == T){
    #--- 'amount' of weights (see the function GCVscore)
    aow <- sum(Wtot[, ])/noe  # 0 < aow <= 1
    #---
    while(tol > TolZ && nit < MaxIter){
      nit <- nit + 1
      DCTy <- dctn(Wtot * (y - z) + z)
      # if (isauto==T && log2(nit)%%1==0){
      #---
      # The generalized cross-validation (GCV) method is used.  We seek the
      # smoothing parameter s that minimizes the GCV score i.e. s =
      # Argmin(GCVscore).  Because this process is time-consuming, it is
      # performed from time to time (when nit is a power of 2)
      #---
      # fminbnd(log10(sMinBnd),log10(sMaxBnd),opt) end
      z <- RF * idctn(Gamma * DCTy) + (1 - RF) * z
      
      # if no weighted/missing data => tol=0 (no iteration)
      tol <- isweighted * Norm((z0 - z)/Norm(z))
      
      z0 <- z
    }  # re-initialization
    exitflag <- nit < MaxIter
    
    if(isrobust == T){
      #-- Robust Smoothing: iteratively re-weighted process
      #--- average leverage
      h <- sqrt(1 + 16 * s)
      h <- sqrt(1 + h)/sqrt(2)/h
      h <- h^N
      #--- take robust weights into account
      Wtot <- W * RobustWeights(y - z, IsFinite, h)
      #--- re-initialize for another iterative weighted process
      isweighted <- TRUE
      tol <- 1
      nit <- 0
      #---
      RobustStep <- RobustStep + 1
      RobustIterativeProcess <- RobustStep < 4
    }else{
      RobustIterativeProcess <- FALSE
    }
  }  # stop the whole process
  
  
  
  gcv <- function(p) {
    # Search the smoothing parameter s that minimizes the GCV score
    #---
    s <- 10^p
    Gamma <- 1/(1 + s * Lambda^2)
    #--- RSS = Residual sum-of-squares
    if (aow > 0.9) {
      # aow = 1 means that all of the data are equally weighted very much
      # faster: does not require any inverse DCT
      RSS <- Norm(DCTy * (Gamma - 1))^2
    } else {
      # take account of the weights to calculate RSS:
      yhat <- idctn(Gamma * DCTy)
      RSS <- Norm(sqrt(Wtot[IsFinite]) * (y[IsFinite] - yhat[IsFinite]))^2
    }
    #---
    TrH <- sum(Gamma[])
    GCVscore <- RSS/nof/(1 - TrH/noe)^2
    return(GCVscore)
  }
  return(z)
}


# function peakdet function to detect all peaks in a vector
peakdet <- function(v, delta) {
  # v <- fun delta <- param
  x <- c((1:length(v)))
  
  if (length(v) != length(x)) {
    print("Input vectors v and x must have same length")
  }
  
  if ((length(delta) > 1)) {
    print("Input argument DELTA must be a scalar")
  }
  
  if (delta <= 0) {
    print("Input argument DELTA must be positive")
  }
  
  mn <- Inf
  mx <- -Inf
  mnpos <- NaN
  mxpos <- NaN
  
  maxtab <- data.frame(mxpos = rep(NA, length(v)), mx = rep(NA, length(v)))
  mintab <- data.frame(mnpos = rep(NA, length(v)), mn = rep(NA, length(v)))
  
  lookformax <- 1
  
  for (i in 1:length(v)) {
    
    this <- v[i]
    if (this > mx) {
      mx <- this
      mxpos <- x[i]
    }
    if (this < mn) {
      mn <- this
      mnpos <- x[i]
    }
    
    if (lookformax == 1) {
      if (this < mx - delta) {
        # temp <- data.frame(mxpos,mx) maxtab <- rbind(maxtab,temp)
        maxtab[i, 1] <- mxpos
        maxtab[i, 2] <- mx
        mn <- this
        mnpos <- x[i]
        lookformax <- 0
      }
    } else if (this > mn + delta) {
      # temp2 <- data.frame(mnpos,mn) mintab <- rbind(mintab,temp2)
      mintab[i, 1] <- mnpos
      mintab[i, 2] <- mn
      mx <- this
      mxpos <- x[i]
      lookformax <- 1
    }
  }
  
  maxtab <- filter(maxtab, !is.na(maxtab$mxpos))
  mintab <- filter(mintab, !is.na(mintab$mnpos))
  
  tabs <- list(maxtab, mintab)
  return(tabs)
}