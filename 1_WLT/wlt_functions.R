WLT <- function(BEP, mother_wavelet, wlt_parameter, delta_freq, 
    signif_level, ID) {
    
  ### functions of Bedforms ATM
  
  # function significance
  significance <- function(motherWavelet, wltParameter, signalVariance, nSignalData, 
                           sampleFreq, deltaFreq, signifLevel, analysisMethod, lag1, wltScale, 
                           lowerScale, upperScale) {
    
    stringMotherWlt <- motherWavelet
    if (stringMotherWlt == "MORLET") {
      
      degreOfFrdm <- 2
      fourierFactor <- (4 * pi)/(wltParameter + sqrt(2 + wltParameter^2))  # Scale-->Fourier [Sec.3h]
      empFactors = c(2, -1, -1, -1)
      if (wltParameter == 6) {
        empFactors[2:4] = c(0.776, 2.32, 0.6)
      }
    }
    
    if (stringMotherWlt == "MEXICANHAT") {
      degreOfFrdm <- 1
      fourierFactor <- 2 * pi/sqrt(wltParameter + 0.5)
      empFactors = c(1, -1, -1, -1)
      if (wltParameter == 2) {
        empFactors[2:4] = c(3.541, 1.43, 1.4)
      }
      if (wltParameter == 6) {
        empFactors[2:4] = c(1.966, 1.37, 0.97)
      }
    }
    
    # Calculator
    nWltScale <- length(wltScale) - 1
    # (variance=1 for the normalized SST)
    signalVariance <- signalVariance
    dj <- log(wltScale[2]/wltScale[1])/log(2)
    period <- wltScale * fourierFactor
    freqIndex <- sampleFreq/period  # normalized frequency
    minDOF <- empFactors[1]  # Degrees of freedom with no smoothing
    recFactor <- empFactors[2]  # reconstruction factor
    gammaFactor <- empFactors[3]  # time-decorrelation factor
    dj0 <- empFactors[4]  # scale-decorrelation factor
    redNoiseSpectr <- (1 - lag1^2)/(1 - 2 * lag1 * cos(freqIndex * 2 * 
                                                         pi) + lag1^2)  # from Eq(16)
    redNoiseSpectr <- signalVariance * redNoiseSpectr  # include time-series variance
    signalSignif <- redNoiseSpectr
    
    if (analysisMethod == 0) {
      
      degreOfFrdm <- minDOF
      X <- chi2inv(signifLevel, degreOfFrdm)
      chiSq <- X/degreOfFrdm
      signalSignif <- redNoiseSpectr * chiSq
      return(signalSignif)
    }
    
    if (analysisMethod == 1) {
      
      degreOfFrdm <- nSignalData - wltScale  # rearanged from the original code
      if (length(degreOfFrdm) == 1) {
        degreOfFrdm <- replicate(nWltScale + 1, 0) + degreOfFrdm
      }
      trunc <- which(degreOfFrdm < 1)
      
      if (length(trunc) > 0) {
        degreOfFrdm[trunc] = replicate(length(trunc), 1)
      }
      degreOfFrdm <- minDOF * sqrt(1 + (degreOfFrdm * sampleFreq/gammaFactor/wltScale)^2)  #Eq (23)
      
      trunc <- which(degreOfFrdm < minDOF)
      
      if (length(trunc) > 0) {
        degreOfFrdm[trunc] <- minDOF * replicate(length(trunc), 1)
      }
      
      iterm <- c(1:(nWltScale + 1))
      funk_singnif <- function(iterm) {
        Y <- chi2inv(signifLevel, degreOfFrdm[iterm])
        chiSq <- Y/degreOfFrdm[iterm]
        signalSignif <- redNoiseSpectr[iterm] * chiSq
        return(signalSignif)
      }
      signalSignif <- sapply(iterm, funk_singnif)
      
      return(signalSignif)
    }
    
    if (analysisMethod == 2) {
      s1 <- lowerScale
      s2 <- upperScale
      wltScale <- data.frame(wltScale)
      avgx <- as.numeric(rownames(wltScale))
      avgx <- data.frame(avgx, wltScale)
      avgx <- filter(avgx, wltScale >= s1 && wltScale <= s2)
      avrg <- avgx[, 1]
      navg <- length(avrg)
      if (navg == 0) {
        print("error scales")
      }
      
      Savg <- 1/sum(1/wltScale[avrg, ])
      Smid <- exp((log(s1) + log(s2))/2)
      redNoiseSpectr <- Savg * sum(redNoiseSpectr[avrg]/wltScale[avrg, 
      ])
      degreOfFrdm <- (minDOF * navg * Savg/Smid) * sqrt(1 + (navg * dj/dj0)^2)
      Z <- chi2inv(signifLevel, degreOfFrdm)
      chiSq <- Z/degreOfFrdm
      signalSignif <- (dj * sampleFreq/recFactor/Savg) * redNoiseSpectr * 
        chiSq
      
      return(signalSignif)
    }
  }
  
  
  # function waveletTransform
  waveletTransform <- function(motherWavelet, wltParameter, signalData, sampleFreq, 
                               deltaFreq) {
    
    nSignalData <- length(signalData)
    s0 <- 2 * sampleFreq
    largestScale <- floor(log2(nSignalData * sampleFreq/s0))/deltaFreq
    wltScale <- s0 * 2^(deltaFreq * (0:largestScale))
    df1 <- deltaFreq * c(0:largestScale)
    wltScale <- s0 * 2^df1
    
    # Normalizing with respect to the mean
    padedSignal <- signalData - mean(signalData)
    padedSignal <- data.frame(padedSignal)
    
    # Padding with zeros
    nSignalBase2 <- nextpow2(nSignalData)
    zeros <- c(1:(2^(nSignalBase2) - nSignalData)) * 0
    zeros <- data.frame(zeros)
    colnames(zeros) <- colnames(padedSignal)
    padedSignal <- rbind(padedSignal, zeros)  # padedSignal, zeros(1,2^(nSignalBase2)-nSignalData)]  #Signal padded with zeros
    nPaded <- nrow(padedSignal)  # Length of the signal padded with zeros 
    
    # Estimate angular frequency - equation (5) Compo's angular frequency
    fourierFreq <- c(1:(nPaded/2))
    fourierFreq <- fourierFreq * ((2 * pi)/(nPaded * sampleFreq))
    cc <- c(floor((nPaded - 1)/2):1)
    fourierFreq = c(0, fourierFreq, -fourierFreq[cc])
    
    # Fourier transform of the paded signal
    ftPaded_Signal = fft(as.numeric(padedSignal$padedSignal))
    
    wlt0 <- list()
    condition1 <- motherWavelet == "MORLET"
    condition2 <- motherWavelet == "MEXICANHAT"
    
    ff <- data.frame(fourierFreq)
    ff <- mutate(ff, V2 = ifelse(ff$fourierFreq > 0, 1, 0))
    ff <- ff[, 2]
    
    # Estimate the Parseval product system.time( for (iScale in
    # 1:(largestScale+1)){
    funk_0 <- function(y) {
      daughterWLT <- list()
      
      if (condition1) {
        centralFreq <- wltParameter
        nFourierFreq <- length(fourierFreq)
        normalizFactor <- sqrt(wltScale[y] * fourierFreq[2]) * (pi^(-0.25)) * 
          sqrt(nFourierFreq)
        expArgument <- -(wltScale[y] * fourierFreq - centralFreq)^2/2 * 
          ff
        daughterWavelet <- normalizFactor * exp(expArgument)
        daughterMorlet <- daughterWavelet * (ff)
        fourierFMorlet <- (4 * pi)/(centralFreq + sqrt(2 + centralFreq^2))
        coneMorlet <- fourierFMorlet/sqrt(2)
        
        daughterWavelet <- daughterMorlet
        fourierFactor <- fourierFMorlet
        coneOfInfl <- coneMorlet
        daughterWLT[[y]] <- daughterWavelet
        rere <- list(a = daughterWavelet, b = fourierFactor, c = coneOfInfl)
        return(rere)
      } else if (condition2) {
        normalizFactor <- sqrt(2 * pi/sampleFreq)/sqrt(gamma(wltParameter + 
                                                               0.5)) * wltScale[y]^0.5
        expArgument <- -(wltScale[y] * fourierFreq)^2/2
        daughterDOG <- -normalizFactor * ((0 + (0 + (0+1i)))^wltParameter) * 
          ((wltScale[y] * fourierFreq)^wltParameter) * exp(expArgument)
        fourierFDOG <- 2 * pi/sqrt(wltParameter + 0.5)
        coneDOG <- fourierFDOG/sqrt(2)
        
        daughterWavelet <- daughterDOG
        fourierFactor <- fourierFDOG
        coneOfInfl <- coneDOG
        daughterWLT[[y]] <- daughterWavelet
      } else {
        print("Wrong MotherWavelet")
      }
    }
    
    scale <- c(1:(largestScale + 1))
    re_1 <- lapply(scale, funk_0)
    re_2 <- transpose(re_1)
    daughter <- as.list(re_2[[1]])
    fourierFactor <- as.numeric(re_2[[2]][largestScale + 1])
    coneOfInfl <- as.numeric(re_2[[3]][largestScale + 1])
    
    funk <- function(x) {
      re <- ifft(ftPaded_Signal * x)
      return(re)
    }
    wlt0 <- sapply(daughter, funk)
    wlt0 <- t(wlt0)
    
    # timeEnergy = ifft(dautherWavelet); Take out the zeros of the paded
    # wavelet transform
    wltTransform <- wlt0[, 1:nSignalData]
    
    # Power of the wavelet
    wltPower <- (abs(wltTransform))^2
    
    # Fourier Period of each scale
    fourierPeriod <- fourierFactor * wltScale
    
    # Global Cone of influence
    vec <- c(1e-05, (1:((nSignalData + 1)/2 - 1)), floor((nSignalData/2 - 
                                                            1)):1, 1e-05)
    coneOfInfluence <- coneOfInfl * sampleFreq * vec
    wlts <- list(wltPower, fourierPeriod, coneOfInfluence)
    return(wlts)
  }
  
  ### function statisticsWlt
  statisticsWlt <- function(motherWavelet, wltParameter, signalData, sampleFreq, 
                            deltaFreq, signifLevel, wltPower, lowerScale, upperScale) {
    
    ## General data Estimate the lag1 coefficient for the Red Noise Power
    ## Spectrum for 10 Lag values (-10 to 10) according to matlab reference
    ## it is good for getting the lag of the white noise
    stringMotherWlt <- motherWavelet
    # Autocorellation berechnen
    lag <- 1
    sig <- data.frame()
    j <- c(1:(length(signalData) - 1))
    
    funk_lag <- function(j) {
      temp <- signalData[j] * signalData[j + 1]
      return(temp)
    }
    
    sig <- sapply(j, funk_lag)
    sig <- as.data.frame(sig)
    sum <- sum(sig[, 1])
    norm <- sum/sum(signalData^2)
    lags <- 1
    df <- data.frame(norm, lags)
    acf <- df
    indexLag1 <- filter(acf, lags == 1)
    lag1 <- indexLag1$norm
    nSignalData <- length(signalData)
    signalVariance <- std(signalData)^2
    s0 <- 2 * sampleFreq
    largestScale <- floor(log2(nSignalData * sampleFreq/s0))/deltaFreq
    wltScale <- s0 * 2^(deltaFreq * (0:largestScale))
    # Get the variance of the wavelet function Significance of the power
    # wavelet Significance levels: (variance=1 for the normalized SST)
    # [signalSignif, redNoiseSpectr] =
    # significance(1.0,dt,scale,0,lag1,SIGLVL,-1,mother);
    signalSignif <- significance(motherWavelet, wltParameter, 1, nSignalData, 
                                 sampleFreq, deltaFreq, signifLevel, 0, lag1, wltScale, lowerScale, 
                                 upperScale)
    
    m <- matrix(1, 1, nSignalData)
    wltSignif <- (signalSignif) %*% m
    wltSignif <- as.matrix(wltPower)/wltSignif
    
    # Global spectrum and significance
    globalWs <- signalVariance^1 * (colSums(t(as.matrix(wltPower)))/nSignalData)
    # globalFftps = globalFftps/max(ftSignal); dof = nSignalData-wltScale;
    # it is inside significance function (RONALD)***
    globalSignif <- significance(motherWavelet, wltParameter, signalVariance, 
                                 nSignalData, sampleFreq, deltaFreq, signifLevel, 1, lag1, wltScale, 
                                 lowerScale, upperScale)
    # Scale Average Significance
    wltScale <- data.frame(wltScale)
    
    avgx <- as.numeric(rownames(wltScale))
    avgx <- data.frame(avgx, wltScale)
    avgx <- filter(avgx, wltScale >= lowerScale && wltScale < upperScale)
    avg <- avgx[, 1]
    if (stringMotherWlt == "MORLET") {
      Cdelta <- 2
      if (wltParameter == 6) {
        Cdelta <- 0.776
      }
    }
    
    if (stringMotherWlt == "MEXICANHAT") {
      Cdelta <- 1
      if (wltParameter == 2) {
        Cdelta <- 3.541
      }
      if (wltParameter == 6) {
        Cdelta <- 1.966
      }
    }
    scaleAvg <- as.matrix(wltScale[, 1]) %*% replicate(nSignalData, 1)
    scaleAvg <- as.matrix(wltPower)/scaleAvg
    scaleAvg <- signalVariance * sampleFreq * sampleFreq/Cdelta * colSums(scaleAvg[avg, 
    ])
    wltScale <- wltScale$wltScale
    
    scaleAveSignif <- significance(motherWavelet, wltParameter, signalVariance, 
                                   nSignalData, sampleFreq, deltaFreq, signifLevel, 2, lag1, wltScale, 
                                   lowerScale, upperScale)
    statistics <- list(wltSignif, globalSignif, globalWs, lag1, scaleAveSignif, 
                       scaleAvg)
    return(statistics)
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
  
  # function SMOOTHN function for filtering signal
  SMOOTHN <- function(y, s){
    k <- 0
    sizy <- size(y)
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
  
  # function Robust weights
  RobustWeights <- function(r, I, h) {
    # weights for robust smoothing.
    MAD <- median(abs(r[I] - median(r[I])))  # median absolute deviation
    u <- abs(r/(1.4826 * MAD)/sqrt(1 - h))  # studentized residuals
    c <- 4.685
    W <- (1 - (u/c)^2)^2 * ((u/c) < 1)  # bisquare weights
    # c = 2.385; W = 1./(1+(u/c).^2); # Cauchy weights c = 2.795; W = u<c;
    # # Talworth weights
    W[is.nan(W)] <- 0
    return(W)
  }
  
  # function dctn
  dctn <- function(y) {
    y <- matrix(y)
    siz <- size(y)
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
    siz <- size(y)
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
  
  
  
  
  # read BEP and make required adaptions
  signal_data <- BEP[ ,2]
  nsignal_data <- length(signal_data)
  signal_dataX <- BEP[, 1]
  sample_freq <- BEP[2, 1] - BEP[1, 1]
    
  # exceute Wavelet-Transformation and get return-parameters
  wlts <- waveletTransform(mother_wavelet, wlt_parameter, signal_data, 
      sample_freq, delta_freq)
  wlt_power <- data.frame(wlts[1])
  fourier_period <- data.frame(wlts[2])
  cone_of_influence <- data.frame(wlts[3])
    
  # execute statisticsWLT and get return-parameters
  lower_scale <- 0.001
  upper_scale <- 10000
  statistics <- statisticsWlt(mother_wavelet, wlt_parameter, signal_data, 
      sample_freq, delta_freq, signif_level, wlt_power, lower_scale, 
      upper_scale)
  wlt_signif <- data.frame(statistics[1])
  global_signif <- data.frame(statistics[2])
  global_ws <- data.frame(statistics[3])
  lag1 <- data.frame(statistics[4])
  scale_ave_signif <- data.frame(statistics[5])
  scale_avg <- data.frame(statistics[6])
    
  # calculating the value and coordinates of the peaks in the global
  # wavelet significance of the wavelet transform The peaks that are
  # greater than the global significance
    
  length_ws <- nrow(global_ws)
  xpos_fm <- 0
  print_f_peak <- 0
  peaks_cross <- data.frame()
  x_peakspos <- data.frame()
  y_peakspos <- data.frame()
    
  for (rr in 1:length_ws) {
      if (global_ws[rr, 1] > global_signif[rr, 1]) {
          if (global_ws[rr, 1] > xpos_fm) {
              print_f_peak <- 1
              xpos_fm <- global_ws[rr, 1]
              yposFM <- log2(fourier_period[rr, 1])
          }
      } else if (print_f_peak == 1) {
          peaks_cross <- rbind(peaks_cross, round(2^yposFM))
          x_peakspos <- rbind(x_peakspos, xpos_fm)
          y_peakspos <- rbind(y_peakspos, yposFM)
          print_f_peak <- 0
      }
  }
    
  peaks_all <- data.frame()
  for(i in 2:(nrow(global_ws)-1)){
    if(global_ws[i,1]>global_ws[i+1,1]&global_ws[i,1]>global_ws[i-1,1]){
      df <- data.frame(x=round(fourier_period[i,1]), y=global_ws[i,1])
      peaks_all <- rbind(peaks_all, df)
    }
  }
    
  # return identified peaks
  if(nrow(peaks_cross)!=0){
    peaks_sig <- data.frame(x=peaks_cross[ ,1], y=x_peakspos[ ,1], layer=rownames(peaks_cross))
    return(list(peaks_sig, peaks_all))
  }else{
    print(paste0("WARNING: no significant peaks determined in BEP ", ID))
    return(list(peaks_all))
  }
}

