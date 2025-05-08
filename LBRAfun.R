#fun definition 
{
  
  trun_model <- function(len, Linf, vbk, Lc, La, t0, mL = NULL, nboot = 100){
    len <- len[!is.na(len)]
    lengths <- len[len >= Lc]
    results <- data.frame(meanlen = NA, n = NA, Z = NA, SE = NA, stringsAsFactors =  FALSE)
    cnt = 0
    
    mean.boot <- function(x,i){
      f <- function(x){
        if(is.null(mL) == TRUE){
          mL <- exp(sum(log(x))/length(x))
          mL <- rnorm(1, mL, sd(lengths)/sqrt(length(lengths)))
        }
        u <- function(Z){
          d1 <- ((Linf-La)/(Linf-Lc))^(Z/vbk)
          d2 <- (Z*(Lc-mL)+vbk*(Linf-mL))/(Z*(La-mL)+vbk*(Linf-mL))
          diff <- (d1-d2)^2
          diff
        }
        return(optimize(u, c(0.01,10), tol = 1e-08)$minimum)
      }
      f(lengths[i])
    }
    dd <- boot(lengths, mean.boot, R = nboot)
    mL <- mean(lengths)
    mLsd <- sd(lengths)
    cnt <- cnt+1
    results <- list(mL = mL, mLsd = mLsd, lengths = length(lengths), meanZ = dd[[1]], Zsd = apply(dd$t,2,sd))
    
    return(results)
  }
  
  # calculate equilibrium abundance
  calc_abund <- function(ages, S_a, M, FM, R0){
    N_a <- R0
    
    Fmat <- NA
    for(i in 1:length(S_a)){
      Fmat[i] <- S_a[i]*FM
    }
    
    dt <- diff(ages)
    
    for (i in 2:length(ages)){
      if (i < length(ages))
        N_a[i] <- N_a[i-1] * exp(-(M+Fmat[i-1])*dt[i-1])
      if (i == length(ages))
        N_a[i] <- N_a[i-1] * exp(-(M+Fmat[i-1])*dt[i-1])/(1-exp(-(M+Fmat[i-1])*dt[i-1]))
    }
    return(N_a)
  }
  # calculate per recruit and SPR
  # ref - FALSE outputs SPR, ref = value between 0 and 1, used with uniroot to find F at which SPR = ref
  calc_ref <- function(ages, Mat_a, W_a, M, FM, S_a, ref = FALSE){
    dt <- diff(ages)[1]
    
    Na0 <- calc_abund(ages = ages, M = M, FM = 0, S_a = S_a, R0 = 1)
    NaF <- calc_abund(ages = ages, M = M, FM = FM, S_a = S_a, R0 = 1)
    
    SB0 <- sum(Na0*Mat_a*W_a)*dt
    SBF <- sum(NaF*Mat_a*W_a)*dt
    
    ratio <- SBF/SB0
    if(ref == FALSE) return(ratio)
    if(ref != FALSE){
      diff <- ref - ratio
      return(diff)
    }
  }
  
  # calculate YPR
  calc_msy <- function(FM, ages, M, R0, W_a, S_a){
    dt <- diff(ages)[1]
    Nage <- calc_abund(ages = ages, M = M, FM = FM, R0 = R0, S_a)
    YPR <- sum(Nage*W_a*(1-exp(-M-FM*S_a))*(FM*S_a)/(M+FM*S_a))*dt
    return(YPR)
  }
  
  # calculate spawning biomass at equilibrium
  calc_SSB <- function(ages, S_a, Mat_a, W_a, M, FM, R0){
    dt <- diff(ages)[1]
    Nage <- calc_abund(ages = ages, S_a = S_a, M = M, FM = FM, R0 = R0)
    SSB <- sum(Nage*Mat_a*W_a)*dt
    return(SSB)
  }
  # calculate per recruit and SPR
  # ref - FALSE outputs SPR, ref = value between 0 and 1, used with uniroot to find F at which SPR = ref
  calc_ref <- function(ages, Mat_a, W_a, M, FM, S_a, ref = FALSE){
    dt <- diff(ages)[1]
    
    Na0 <- calc_abund(ages = ages, M = M, FM = 0, S_a = S_a, R0 = 1)
    NaF <- calc_abund(ages = ages, M = M, FM = FM, S_a = S_a, R0 = 1)
    
    SB0 <- sum(Na0*Mat_a*W_a)*dt
    SBF <- sum(NaF*Mat_a*W_a)*dt
    
    ratio <- SBF/SB0
    if(ref == FALSE) return(ratio)
    if(ref != FALSE){
      diff <- ref - ratio
      return(diff)
    }
  }
  
  # calculate YPR
  calc_msy <- function(FM, ages, M, R0, W_a, S_a){
    dt <- diff(ages)[1]
    Nage <- calc_abund(ages = ages, M = M, FM = FM, R0 = R0, S_a)
    YPR <- sum(Nage*W_a*(1-exp(-M-FM*S_a))*(FM*S_a)/(M+FM*S_a))*dt
    return(YPR)
  }
}
LBRA_function <- function(
    len,
    LF_data,
    Linf,
    vbk,
    t0,
    Amax,
    LWa,
    LWb,
    Lmat,
    Lc,
    La,
    nboot = 200,
    sim_run = 100,
    start_age = 0,
    seed = 1234,
    CVlen = 0.07,
    nseason = 1 
)
{
  set.seed(seed)
  Nyears <- dim(LF_data)[2] # number of years of data
  
  # outputs to save for each simulation
  SPR <- matrix(nrow = sim_run, ncol = Nyears)
  YPR <- matrix(nrow = sim_run, ncol = Nyears)
  Fmsy <- matrix(nrow = sim_run, ncol = Nyears)
  FFmsy <- matrix(nrow = sim_run, ncol = Nyears)
  SSB_eq <- matrix(nrow = sim_run, ncol = Nyears)
  Bmsy <- matrix(nrow = sim_run, ncol = Nyears)
  BBmsy <- matrix(nrow = sim_run, ncol = Nyears)
  SPRmsy <- matrix(nrow = sim_run, ncol = Nyears)
  SSB <- matrix(nrow = sim_run, ncol = Nyears)
  B <- matrix(nrow = sim_run, ncol = Nyears)
  AvgN <- matrix(nrow = sim_run, ncol = Nyears)
  AvgB <- matrix(nrow = sim_run, ncol = Nyears)
  FM_list <- matrix(nrow = sim_run, ncol = Nyears)
  Z_list <- matrix(nrow = sim_run, ncol = Nyears)
  Zsd_list <- matrix(nrow = sim_run, ncol = Nyears)
  mL_list <- matrix(nrow = sim_run, ncol = Nyears)
  mLsd_list <- matrix(nrow = sim_run, ncol = Nyears)
  M_list <- vector(length = sim_run)
  beta_list <- vector(length = sim_run)
  
  # risk counters
  Fishedcount <- 0 # biomass
  Fishingcount <- 0 # fishing mortality
  SPRcount <- 0 # SPR
  
  
  for(i in 1:sim_run){
    # mean length, number of lengths above Lc, Z, sd
    res <- lapply(1:Nyears, function(x) trun_model(len = len[[x]], Linf = Linf, vbk = vbk, Lc = Lc, La = La[x], nboot = nboot))
    # warning - this model has trouble running if Linf < La
    
    ## probabilistic mortality rates ##
    # mean length - normal distribution
    mL <- sapply(1:Nyears, function(x) res[[x]]$mL)
    mLsd <- sd(mL)
    mL_probs <- sapply(1:Nyears, function(x) rnorm(1, mL[[x]], mLsd))
    # total mortality - reapply new mean lengths into model
    res <- lapply(1:Nyears, function(x) trun_model(len = len[[x]], Linf = Linf, vbk = vbk, mL = mL_probs[x], Lc = Lc, La = La[x], nboot = nboot))
    Z <- sapply(1:Nyears, function(x) res[[x]]$meanZ)
    
    # natural mortality
    # natural mortality
    M_init <- -log(0.05)/Amax # 5% survivorship 
    a_lambda <- -log(0.001)/M_init # theoretical age - based on 0.1%
    a_delta <- a_lambda-Amax # age interval
    one_beta <- rexp(1, -log(0.001)/a_delta) # exponential distribution
    M_final <- -log(0.05)/(one_beta+Amax)
    
    # fishing mortality
    FM <- sapply(1:Nyears, function(x) Z[[x]]-M_final)
    FM[FM<0] <- 0 # FM trap
    
    
    ## life history ##
    # length bins
    binwidth <- as.numeric(rownames(LF_data)[2])-as.numeric(rownames(LF_data)[1])
    mids <- as.numeric(rownames(LF_data)) # length bins
    highs <- mids+(binwidth/2) # length bins high
    lows <- mids-(binwidth/2) # length bins low
    
    # length and weight at age
    Ages <- seq(from = start_age, to = Amax, by = (1/nseason)) # all ages including age at 0 if start_ages = 0 and seasonality
    L_a <- Linf*(1-exp(-vbk*(Ages-t0))) # length at age
    W_a <- LWa*L_a^LWb # weight at age
    Age_length <- length(Ages)
    
    # probability length bin given age
    lbprobs <- function(mnl,sdl) return(pnorm(highs, mnl, sdl)-pnorm(lows, mnl, sdl)) # setting up the probability between high and low length bins
    vlprobs <- Vectorize(lbprobs, vectorize.args = c("mnl", "sdl"))
    plba <- t(vlprobs(L_a, L_a*CVlen)) # probability of being in length bin given length at age and CV of length
    plba <- plba/rowSums(plba)
    plba[is.na(plba)] <- 0 # need this especially if start_age = 0
    
    # knife edge selectivity
    S_l <- rep(0, length(mids))
    S_l[mids >= Lc] <- 1
    S_a <- colSums(t(plba)*S_l)
    S_a <- S_a/max(S_a)
    
    # knife edge maturity
    Mat_l <- rep(0, length(mids))
    Mat_l[mids >= Lmat] <- 1
    Mat_a <- colSums(t(plba)*Mat_l)
    Mat_a <- Mat_a/max(Mat_a)
    
    
    ## abundance calc ##
    N_a <- matrix(nrow = Nyears, ncol = Age_length) # fished equilibrium calculation
    # fishing mortality
    F_a <- matrix(nrow = Nyears, ncol = Age_length)
    for(y in 1:Nyears){
      for(a in 1:Age_length){
        F_a[y,a] <- FM[y]*S_a[a]
      }
    }
    
    # need a better way to calculate recruitment, can replace here
    R0 <- 1 
    R_t <- vector(length = Nyears)
    rmaxBH <- 1000
    betaBH <- 1 # constant recruitment
    
    # year 1
    for(a in 1:Age_length){
      if(a == 1) N_a[1,a] <- R0
      if(a > 1 & a < Age_length) N_a[1,a] <- N_a[1,a-1]*exp(-M_final-F_a[y,a])
      if(a == Age_length) N_a[1,a] <- (N_a[1,a-1]*exp(-M_final-F_a[y,a]))/(1-exp(-M_final-F_a[y,a]))
    }
    SSB[i,1] <- sum(N_a[1,]*W_a*Mat_a)
    
    # numbers at age
    for(y in 2:Nyears){
      R_t[y] <- rmaxBH*SSB[i,y-1]/(betaBH*SSB[i,y-1])
      N_a[y,1] <- R_t[y]
      for(a in 2:Age_length){
        if(a < Age_length) N_a[y,a] <- N_a[y-1,a-1]*exp(-(M_final+F_a[a-1]))
        if(a == Age_length) N_a[y,a] <- N_a[y-1,a-1]*exp(-(M_final+F_a[a-1])/(1-exp(-(M_final+F_a[a-1]))))
      }
      SSB[i,y] <- sum(N_a[y,]*W_a*Mat_a)
    }
    
    
    ## Calculations ##
    # spawning biomass
    # SSB[i,] <- sapply(1:Nyears, function(x) sum(N_a[x,]*Mat_a*W_a))
    # biomass
    B[i,] <- sapply(1:Nyears, function(x) sum(N_a[x,]*W_a))
    
    # average numbers
    for(y in 1:Nyears){
      N_temp <- N_a[y,]*(1-exp(-Z[y]))/Z[y]
      B_temp <- N_temp*W_a
      AvgN[i,y] <- sum(N_temp)
      AvgB[i,y] <- sum(B_temp) 
    }
    
    
    ## Sustainability Analysis 
    SPR[i,] <- sapply(1:Nyears, function(x) calc_ref(ages = Ages, Mat_a = Mat_a, W_a = W_a, M = M_final, S_a = S_a, FM = FM[x]))
    YPR[i,] <- sapply(1:Nyears, function(x) calc_msy(FM = FM[x], R0 = R0, ages = Ages, M = M_final, W_a = W_a, S_a = S_a))
    Fmsy[i,] <- sapply(1:Nyears, function(x) optimize(calc_msy, ages = Ages, M = M_final, R0 = R0, W_a = W_a, S_a = S_a, lower = 0, upper = 10, maximum = TRUE)$maximum)
    FFmsy[i,] <- sapply(1:Nyears, function(x) FM[x]/Fmsy[i,x])
    Bmsy[i,] <- sapply(1:Nyears, function(x) calc_SSB(ages = Ages, S_a = S_a, Mat_a = Mat_a, W_a = W_a, M = M_final, FM = Fmsy[i,x], R0 = R0))
    SSB_eq[i,] <- sapply(1:Nyears, function(x) calc_SSB(ages = Ages, S_a = S_a, Mat_a = Mat_a, W_a = W_a, M = M_final, FM = FM[x], R0 = R0))
    BBmsy[i,] <- sapply(1:Nyears, function(x) SSB_eq[i,x]/Bmsy[i,x])
    SPRmsy[i,] <- sapply(1:Nyears, function(x) calc_ref(ages = Ages, Mat_a = Mat_a, W_a = W_a, M = M_final, S_a = S_a, FM = Fmsy[i,x]))
    
    ## save others
    FM_list[i,] <- FM
    M_list[i] <- M_final
    Z_list[i,] <- Z
    # Zsd_list[i,] <- Zsd
    mL_list[i,] <- mL_probs
    # mLsd_list[i,] <- mLsd
    beta_list[i] <- one_beta
    
    # Counters
    if(FFmsy[i,Nyears] > 1) Fishingcount <- Fishingcount+1/sim_run
    if(BBmsy[i,Nyears] < 1) Fishedcount <- Fishedcount+1/sim_run
    if(SPR[i,Nyears] < 0.4) SPRcount <- SPRcount+1/sim_run
  } # end of sim_run loop
  
  return(list(FM = FM_list,
              M = M_list,
              Z = Z_list,
              # Zsd = Zsd_list,
              mL = mL_list,
              # mLsd = mLsd_list,
              onebeta = beta_list,
              SPR = SPR,
              YPR = YPR,
              Fmsy = Fmsy,
              FFmsy = FFmsy,
              Bmsy = Bmsy,
              BBmsy = BBmsy,
              SSB = SSB,
              B = B,
              SPRmsy = SPRmsy,
              AvgN = AvgN,
              AvgB = AvgB,
              Fishingcount = Fishingcount,
              Fishedcount = Fishedcount,
              SPRcount = SPRcount
  ))
} # end of function
