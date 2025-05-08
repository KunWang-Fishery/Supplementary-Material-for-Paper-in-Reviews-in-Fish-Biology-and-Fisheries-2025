# The LBB function was rewrite based on Froese et al.(2018).
LBBest<-function (lfq,  
          LinfUser = NA, LcutUser = NA, LcUser = NA, LstartUser = NA, 
          MKUser = NA, mmUser = FALSE, GausSel = FALSE, MergeLF = FALSE, 
          n.chains = 3, n.cluster = 3, plot = FALSE, mfrow = NA) 
{
  set.seed(123)
  if (!"Lm50" %in% names(lfq)) {
    stop("Lm50 missing! Please provide the length at maturity as element 'Lm50' in the argument 'lfq'.")
  }
  res<-lfq
  years <- as.numeric(format(lfq$dates, "%Y"))
 if (class(lfq$catch) == "numeric") {
    nrowC <- length(lfq$catch)
    ncolC <- 1
  }
  StartYear <- min(years)
  EndYear <- max(years)
  AllYear <- rep(years, each = nrowC)
  AllLength <- rep(lfq$midLengths, ncolC)
  AllFreq <- as.numeric(lfq$catch)
  Years <- sort(unique(AllYear))
  nYears <- length(years)
  Stock <- ifelse(!is.null(res$stock), res$stock, "unknown stock")
  Comment <- ifelse(!is.null(res$comment), res$comment, "no comment")
  Species <- ifelse(!is.null(res$species), res$species, "unknown species")
  idx <- which(AllFreq != 0)
  AllLength <- AllLength[idx]
  AllYear <- AllYear[idx]
  AllFreq <- AllFreq[idx]
  Ldat <- data.frame(Stock = rep(Stock, nYears), Year = rep(NA, 
                                                            nYears), Linf = rep(NA, nYears), Linf.lcl = rep(NA, nYears), 
                     Linf.ucl = rep(NA, nYears), Lc = rep(NA, nYears), Lc.lcl = rep(NA, 
                                                                                    nYears), Lc.ucl = rep(NA, nYears), Lmean = rep(NA, 
                                                                                                                                   nYears), r.alpha = rep(NA, nYears), r.alpha.lcl = rep(NA, 
                                                                                                                                                                                         nYears), r.alpha.ucl = rep(NA, nYears), r.GLmean = rep(NA, 
                                                                                                                                                                                                                                                nYears), r.SD = rep(NA, nYears), MK = rep(NA, nYears), 
                     MK.lcl = rep(NA, nYears), MK.ucl = rep(NA, nYears), FK = rep(NA, 
                                                                                  nYears), FK.lcl = rep(NA, nYears), FK.ucl = rep(NA, 
                                                                                                                                  nYears), ZK = rep(NA, nYears), ZK.lcl = rep(NA, nYears), 
                     ZK.ucl = rep(NA, nYears), FM = rep(NA, nYears), FM.lcl = rep(NA, 
                                                                                  nYears), FM.ucl = rep(NA, nYears), r.Lopt = rep(NA, 
                                                                                                                                  nYears), BB0 = rep(NA, nYears), BB0.lcl = rep(NA, 
                                                                                                                                                                                nYears), BB0.ucl = rep(NA, nYears), YR = rep(NA, 
                                                                                                                                                                                                                             nYears), YR.lcl = rep(NA, nYears), YR.ucl = rep(NA, 
                                                                                                                                                                                                                                                                             nYears), perc.mat = rep(NA, nYears), L95 = rep(NA, 
                                                                                                                                                                                                                                                                                                                            nYears))
  df <- data.frame(AllYear, AllLength, AllFreq)
  names(df) <- c("Year", "Length", "Freq")
  AG  <- function(dat) { # where dat contains dat$Year, dat$Length in cm, dat$CatchNo
    
    # aggregate normalized annual LFs by weighing with square root of sample size
    # get sum of frequencies per year
    sum.Ny  <- aggregate(Freq~Year,dat,sum)$Freq  
    # get the sqrt of the sum of frequencies for every year
    sqrt.Ny <- sqrt(sum.Ny) 
    # get highest frequency in each year
    max.Ny <- aggregate(Freq~Year,dat,max)$Freq
    # get Number of Length bins in each year
    binsN <- aggregate(Freq~Year,dat,length)$Freq    
    # create vectors for sqrt.Ni and sum.Ni to weigh LF data
    sqrt.Ni = rep(sqrt.Ny,binsN)
    sum.Ni = rep(sum.Ny,binsN)
    #Do weighing
    # Divide all years by sum.Ni and multiply by sqrt.Ni
    LF.w = dat$Freq/sum.Ni*sqrt.Ni  
    # Aggregate
    LF = aggregate(LF.w, by=list(dat$Length),FUN=sum)
    # Add correct column names
    colnames(LF) <- c("Length","Freq")         
    return(LF)
  } #end of aggregate function
  BH <- function(AllLength,Linf,MK,FK,GausSel,selpar1,selpar2) {
    if(GausSel==F) {
      r.Lc     <- selpar1
      r.alpha  <- selpar2 
      Lx       <- AllLength[AllLength >= Linf*(r.Lc-4.59/r.alpha)][1]
    } else if(GausSel==T) {
      r.GLmean <- selpar1
      r.SD     <- selpar2
      Lx       <- AllLength[AllLength >= Linf*(r.GLmean-3*r.SD)][1]
    }
    class.width  <- median(diff(sort(unique(AllLength))))
    FM <- FK/MK
    
    r            <- vector() # auxilliary reduction factor
    G            <- vector() # product of reduction factors
    SL.bh        <- vector() # selection at length
    YR1.2        <- vector() # relative yield per recruit per length class
    CPUER1.2     <- vector() # relative CPUE per recruit per length class
    B1.2         <- vector() # relative unexploited biomass per recruit by length class
    L.bh         <- seq(from=Lx, to=Linf, by=class.width) # lengths to be considered
    r.L.bh       <-  L.bh / Linf # standardized lengths
    
    # calculate selection, Y'/R and CPUE'/R for every length class
    for(o in 1 : length(r.L.bh)) { 
      if(GausSel==F) {
        if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc))), # mean selection in length class
                                                  1/(1+exp(-r.alpha*(r.L.bh[o+1]-r.Lc)))))
        } else SL.bh[o] <- 1/(1+exp(-r.alpha*(r.L.bh[o]-r.Lc)))
      } else if(GausSel==T) { # gill net selection 
        if(o<length(r.L.bh)) { SL.bh[o] <- mean(c(exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2))), # mean selection in length class
                                                  exp(-((r.L.bh[o+1]-r.GLmean)^2/(2*r.SD^2)))))
        } else SL.bh[o] <- exp(-((r.L.bh[o]-r.GLmean)^2/(2*r.SD^2)))
      } # end of calculation of selectivity loop
      
      if(o<length(r.L.bh)) {
        r[o]       <- (1-r.L.bh[o+1])^(FK*SL.bh[o])/(1-r.L.bh[o])^(FK*SL.bh[o]) 
        G[o]       <- prod(r[1:o]) }
      if(o==1) {
        YR1.2[o] <-(FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                       (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                    (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) -
          (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                  (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                             (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o] 
      } else if(o==length(r.L.bh)) {
        YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                        (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                     (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] 
      } else {
        YR1.2[o] <- (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/
                                                                                        (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                                     (1-r.L.bh[o])^3/(1+3/(MK+FK*SL.bh[o])))) * G[o-1] -
          (FM*SL.bh[o]/(1+FM*SL.bh[o])*(1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/
                                                                                  (MK+FK*SL.bh[o]))+3*(1-r.L.bh[o+1])^2/(1+2/(MK+FK*SL.bh[o]))-
                                                             (1-r.L.bh[o+1])^3/(1+3/(MK+FK*SL.bh[o]))))*G[o]              
      } # end of loop to calculate yield per length class
      
      CPUER1.2[o] <- YR1.2[o] / FM # CPUE/R = Y/R divided by F/M
      
      if(o<length(r.L.bh)) {
        B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                        (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)) -
                      (1-r.L.bh[o+1])^MK*(1-3*(1-r.L.bh[o+1])/(1+1/MK)+3*(1-r.L.bh[o+1])^2/
                                            (1+2/MK)-(1-r.L.bh[o+1])^3/(1+3/MK)))*SL.bh[o]
      } else {
        B1.2[o] <- ((1-r.L.bh[o])^MK*(1-3*(1-r.L.bh[o])/(1+1/MK)+3*(1-r.L.bh[o])^2/
                                        (1+2/MK)-(1-r.L.bh[o])^3/(1+3/MK)))*SL.bh[o]
      }
    } # end of B&H loop through length classes
    BB0   <- sum(CPUER1.2)/sum(B1.2)
    YR    <- sum(YR1.2)
    if(BB0 < 0.25) YR <- YR * BB0 / 0.25 # reduce YR if recruitment and thus productivity is reduced
    return(list(BB0,YR))
    
  } # end of BH function
  LF.all <- AG(dat = df)
  LF.all$Freq = LF.all$Freq/max(LF.all$Freq)
  LF.all <- LF.all[which(LF.all$Freq > 0)[1]:length(LF.all$Length), 
  ]
  LF.all <- LF.all[1:which(LF.all$Length == max(LF.all$Length[LF.all$Freq > 
                                                                0])), ]
  n.LF.all <- length(LF.all$Length)
  Lmax <- LF.all$Length[n.LF.all]
  Lmax.med <- median(as.numeric(by(rep(res$midLengths, ncolC)[as.numeric(res$catch) > 
                                                                0], rep(res$dates, each = nrowC)[as.numeric(res$catch) > 
                                                                                                   0], max)))/10
  L10 <- LF.all$Length[which(LF.all$Freq > 0.1)[1]]
  L90 <- LF.all$Length[which(LF.all$Freq > 0.9)[1]]
  Lc.st <- ifelse(is.na(LcUser), (L10 + L90)/2, LcUser)
  alpha.st <- max(0.01,abs(-log(1/LF.all$Freq[which(LF.all$Freq > 0.1)[1]])/(L10 - 
                                                                  Lc.st)))
  Linf.st <- max(LF.all$Length)
  Lmean.st <- sum(LF.all$Length[LF.all$Length >= Lc.st] * LF.all$Freq[LF.all$Length >= 
                                                                        Lc.st])/sum(LF.all$Freq[LF.all$Length >= Lc.st])
  MK.st <- ifelse(is.na(MKUser), 1.5, MKUser)
  ZK.st <- (Linf.st - Lmean.st)/(Lmean.st - Lc.st)
  FK.st <- ifelse((ZK.st - MK.st) > 0, ZK.st - MK.st, 0.3)
  if (!is.na(LstartUser)) {
    Lstart <- LstartUser
  }
  else {
    Lstart <- (alpha.st * Lc.st - log(1/0.95 - 1))/alpha.st
    Lstart.i <- which(LF.all >= Lstart)[1]
    Lmax.i <- length(LF.all$Length)
    peak.i <- which.max(LF.all$Freq)
    if (Lstart.i < (peak.i + 1)) 
      Lstart <- LF.all$Length[peak.i + 1]
    if ((Lmax.i - Lstart.i) < 4) 
      Lstart <- LF.all$Length[Lstart.i - 1]
  }
  L.L <- LF.all$Length[LF.all$Length >= Lstart & LF.all$Length < 
                         Linf.st]
  L.Freq <- LF.all$Freq[LF.all$Length >= L.L[1] & LF.all$Length < 
                          Linf.st]
  if (length(L.L) < 4) {
    plot(x = LF.all$Length, y = LF.all$Freq, bty = "l", main = Stock)
    lines(x = c(Lstart, Lstart), y = c(0, 0.9 * max(LF.all$Freq)), 
          lty = "dashed")
    text(x = Lstart, y = max(LF.all$Freq), "Lstart")
    lines(x = c(Linf.st, Linf.st), y = c(0, 0.9 * max(LF.all$Freq)), 
          lty = "dashed")
    text(x = Linf.st, y = max(LF.all$Freq), "Lmax")
    stop("Too few fully selected data points: set Lstart.user\n")
  }
  sum.L.Freq <- sum(L.Freq)
  L.Freq <- L.Freq/sum.L.Freq
  if (is.na(LinfUser)) {
    Linf.mod <- nls(L.Freq ~ ((Linf - L.L)/(Linf - Lstart))^ZK/sum(((Linf - 
                                                                       L.L)/(Linf - Lstart))^ZK), start = list(ZK = ZK.st, 
                                                                                                               Linf = Linf.st), lower = c(0.5 * ZK.st, 0.999 * Linf.st), 
                    upper = c(1.5 * ZK.st, 1.2 * Linf.st), algorithm = "port")
    ZK.nls <- as.numeric(coef(Linf.mod)[1])
    ZK.nls.sd <- as.numeric(coef(summary(Linf.mod))[, 2][1])
    ZK.nls.lcl <- ZK.nls - 1.96 * ZK.nls.sd
    ZK.nls.ucl <- ZK.nls + 1.96 * ZK.nls.sd
    Linf.nls <- as.numeric(coef(Linf.mod)[2])
    Linf.nls.sd <- as.numeric(coef(summary(Linf.mod))[, 2][2])
    Linf.lcl <- Linf.nls - 1.96 * Linf.nls.sd
    Linf.ucl <- Linf.nls + 1.96 * Linf.nls.sd
  }
  else {
    Linf.nls <- LinfUser
    Linf.nls.sd <- 0.001 * LinfUser
    ZK.mod <- nls(L.Freq ~ exp(ZK * (log(1 - L.L/Linf.nls) - 
                                       log(1 - L.L[1]/Linf.nls)))/sum(exp(ZK * (log(1 - 
                                                                                      L.L/Linf.nls) - log(1 - L.L[1]/Linf.nls)))), start = list(ZK = ZK.st), 
                  lower = c(0.7 * ZK.st), upper = c(1.3 * ZK.st), algorithm = "port")
    ZK.nls <- max(as.numeric(coef(ZK.mod)[1]),MKUser)
    ZK.nls.sd <- as.numeric(coef(summary(ZK.mod))[, 2][1])
    ZK.nls.lcl <- max((ZK.nls - 1.96 * ZK.nls.sd),0)
    ZK.nls.ucl <- ZK.nls + 1.96 * ZK.nls.sd
  }
  AllFreq <- AllFreq[AllLength <= Linf.nls]
  AllYear <- AllYear[AllLength <= Linf.nls]
  AllLength <- AllLength[AllLength <= Linf.nls]
  cat("Running Jags model to fit SL and N distributions for", 
      Species, "\nin", Years, "....\n")
  i = 0
  for (Year in Years) {
    i = i + 1
    if (i > 1 & MergeLF & substr(Stock, start = nchar(Stock) - 
                                 2, stop = nchar(Stock)) != "Sim") {
      AG.yr <- c(Years[i - 1], Year)
    }
    else AG.yr <- Year
    df <- data.frame(AllYear[AllYear %in% AG.yr], AllLength[AllYear %in% 
                                                              AG.yr], AllFreq[AllYear %in% AG.yr])
    names(df) <- c("Year", "Length", "Freq")
    LF.y <- AG(dat = df)
    LF.y$Freq <- LF.y$Freq/sum(LF.y$Freq)
    LF.y <- LF.y[which(LF.y$Freq > 0)[1]:length(LF.y$Length), 
    ]
    LF.y <- LF.y[1:which.max(LF.y$Length[LF.y$Freq > 0]), 
    ]
    L.y <- LF.y$Length
    r.Freq.y <- LF.y$Freq
    r.Freq.y[r.Freq.y == 0] <- min(r.Freq.y[r.Freq.y > 0], 
                                   na.rm = T)/100
    Ldat$Year[i] <- Year
    n.L <- length(L.y)
    Linf.pr <- Linf.nls
    Linf.sd.pr <- ifelse(Linf.nls.sd/Linf.nls < 0.01, Linf.nls.sd, 
                         0.01 * Linf.nls)
    MK.pr <- MK.st
    MK.sd.pr <- ifelse(is.na(MKUser) == TRUE, 0.15, 0.075)
    if (GausSel == FALSE) {
      Lc.pr <- ifelse(is.na(LcUser) == TRUE, 1.02 * Lc.st, 
                      LcUser)
      Lc.sd.pr <- ifelse(is.na(LcUser) == TRUE, 0.1 * Lc.pr, 
                         0.005 * Lc.pr)
      r.max.Freq <- max(r.Freq.y, na.rm = T)
      r.alpha.pr <- max(0.01,abs(-log(r.max.Freq/r.Freq.y[which(r.Freq.y > 
                                                     (0.1 * r.max.Freq))[1]])/(L10/Linf.nls - Lc.st/Linf.nls)))
      r.alpha.sd.pr <- 0.025 * r.alpha.pr
      FK.pr <- ifelse((ZK.nls - MK.st) > 0, ZK.nls - MK.st, 
                      0.3)
      jags.data <- list("r.Freq.y", "L.y", "n.L", "Linf.pr", 
                        "Linf.sd.pr", "Lc.pr", "Lc.sd.pr", "r.alpha.pr", 
                        "r.alpha.sd.pr", "MK.pr", "MK.sd.pr", "FK.pr")
      jags.params <- c("r.alpha.d", "Lc.d", "SL", "xN", 
                       "FK.d", "MK.d", "Linf.d")
      sink("SLNMod.jags")
      cat("\n             model {\n             r.alpha.d_tau  <- pow(r.alpha.sd.pr, -2) \n             r.alpha.d      ~ dnorm(r.alpha.pr,r.alpha.d_tau) \n\n              Lc.d_tau  <- pow(Lc.sd.pr,-2)\n              Lc.d      ~ dnorm(Lc.pr,Lc.d_tau) #       \n\n              MK.d_tau  <-pow(MK.sd.pr, -2) # strong prior on M/K\n              MK.d      ~ dnorm(MK.pr, MK.d_tau)\n\n              Linf.tau  <- pow(Linf.sd.pr,-2) \n              Linf.d    ~ dnorm(Linf.pr,Linf.tau)\n\n              FK.d       ~ dlnorm(log(FK.pr),4) # wide prior range for F/K\n\n              SL[1]       ~ dlogis(0,1000)\n              Freq.pred[1]<-0\n              xN[1]       <-1\n\n              for(j in 2:n.L) {\n               SL[j]<- 1/(1+exp(-r.alpha.d*(L.y[j]/Linf.d-Lc.d/Linf.d))) # selection at length L[j]\n\n               xN[j] <- xN[j-1]*((Linf.d-L.y[j])/(Linf.d-L.y[j-1]))^(MK.d+FK.d*SL[j])\n\n                           Freq.pred[j]<-xN[j]*SL[j]\n\n               # normalize frequencies by dividing by sum of frequencies; multiply with 10 to avoid small numbers and with 1000 for effective sample size\n               r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10*1000\n             }\t\n\n             #><> LIKELIHOOD FUNCTION\n             #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)\n             r.Freq.y[2:n.L] ~ ddirch(r.Freq.pred[2:n.L])  \n\n             } # END OF MODEL\n               ", 
          fill = TRUE)
      sink()
      MODEL = "SLNMod.jags"
      if (n.cluster > 1 & n.cluster == n.chains) {
        jagsfitSLN <- R2jags::jags.parallel(data = jags.data, 
                                            working.directory = NULL, inits = NULL, parameters.to.save = jags.params, 
                                            model.file = paste(MODEL), n.burnin = 300, 
                                            n.thin = 10, seed=100,n.iter = 600, n.chains = n.chains, 
                                            n.cluster = n.cluster)
      }
      else {
        jagsfitSLN <- R2jags::jags(data = jags.data, 
                                   working.directory = NULL, inits = NULL, parameters.to.save = jags.params, 
                                   model.file = paste(MODEL), n.burnin = 300, 
                                   n.thin = 10,n.iter = 600, n.chains = n.chains)
      }
      Ldat$Lc[i] <- median(jagsfitSLN$BUGSoutput$sims.list$Lc.d)
      Ldat$Lc.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d, 
                                 0.025)
      Ldat$Lc.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Lc.d, 
                                 0.975)
      Ldat$Lmean[i] <- sum(L.y[L.y >= Ldat$Lc[i]] * r.Freq.y[L.y >= 
                                                               Ldat$Lc[i]])/sum(r.Freq.y[L.y >= Ldat$Lc[i]])
      Ldat$r.alpha[i] <- median(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d)
      Ldat$r.alpha.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d, 
                                      0.025)
      Ldat$r.alpha.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$r.alpha.d, 
                                      0.975)
      Ldat$MK[i] <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
      Ldat$MK.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d, 
                                 0.025)
      Ldat$MK.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d, 
                                 0.975)
      Ldat$FK[i] <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
      Ldat$FK.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d, 
                                 0.025)
      Ldat$FK.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d, 
                                 0.975)
      FMi <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
      Ldat$FM[i] <- median(FMi)
      Ldat$FM.lcl[i] <- quantile(FMi, 0.025)
      Ldat$FM.ucl[i] <- quantile(FMi, 0.975)
      ZKi <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
      Ldat$ZK[i] <- median(ZKi)
      Ldat$ZK.lcl[i] <- quantile(ZKi, 0.025)
      Ldat$ZK.ucl[i] <- quantile(ZKi, 0.975)
      Ldat$r.Lopt[i] <- 3/(3 + Ldat$MK[i])
      Ldat$Linf[i] <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
      Ldat$Linf.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d, 
                                   0.025)
      Ldat$Linf.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d, 
                                   0.975)
    }
    if (GausSel == TRUE) {
      GLmean.st <- L.y[which.max(r.Freq.y)]
      Lc.pr <- L.y[which(r.Freq.y >= (0.5 * max(r.Freq.y)))][1]
      SD.st <- max(GLmean.st - Lc.pr, 0.25 * GLmean.st)
      cat("Running Jags model to fit SL and N distributions\n")
      n.L <- length(L.y)
      jags.data <- list("n.L", "GLmean.st", "L.y", "SD.st", 
                        "ZK.nls", "r.Freq.y", "Linf.pr", "Linf.sd.pr", 
                        "MK.pr")
      jags.params <- c("GLmean.d", "SD.d", "SL", "xN", 
                       "FK.d", "MK.d", "Linf.d")
      sink("SLNMod.jags")
      cat("\n              model {\n              GLmean.tau <- pow(0.1*GLmean.st,-2) \n              GLmean.d   ~ dnorm(GLmean.st,GLmean.tau)\n\n              SD.tau    <- pow(0.2*SD.st,-2)\n              SD.d      ~ dnorm(SD.st,SD.tau)\n\n              MK.d_tau  <-pow(0.15,-2)\n              MK.d      ~ dnorm(MK.pr,MK.d_tau)\n\n              Linf.tau  <- pow(Linf.sd.pr,-2)\n              Linf.d    ~ dnorm(Linf.pr,Linf.tau)\n\n              FK        <- (ZK.nls-1.5) # ZK overestimated in gillnet selection, used as upper range\n              FK.d      ~ dunif(0,FK)  \n\n              SL[1]~ dlogis(0,1000)\n              Freq.pred[1]<-0\n              xN[1]<-1\n\n              for(j in 2:n.L) {\n                SL[j]<- exp(-((L.y[j]-GLmean.d)^2/(2*SD.d^2)))\n\n                xN[j]<-xN[j-1]*exp((MK.d+FK.d*SL[j])*(log(1-L.y[j]/Linf.d)-log(1-L.y[j-1]/Linf.d)))\n\n                Freq.pred[j]<-xN[j]*SL[j]\n\n                #><> add effective sample size (try 100 typical for LF data)\n                r.Freq.pred[j]<- Freq.pred[j]/sum(Freq.pred)*10000\n              }\t\n\n              #><> LIKELIHOOD FUNCTION\n              #><> Fit observed to predicted LF data using a Dirichlet distribution (more robust in JAGS)\n              r.Freq.y[2:n.L]~ddirch(r.Freq.pred[2:n.L])  \n\n           } # END OF MODEL\n              ", 
          fill = TRUE)
      sink()
      MODEL = "SLNMod.jags"
      if (n.cluster > 1 & n.cluster == n.chains) {
        jagsfitSLN <- R2jags::jags.parallel(data = jags.data, 
                                            working.directory = NULL, inits = NULL, parameters.to.save = jags.params, 
                                            model.file = paste(MODEL), n.burnin = 300, 
                                            n.thin = 10, n.iter = 1000, n.chains = n.chains, 
                                            n.cluster = n.cluster)
      }
      else {
        jagsfitSLN <- R2jags::jags(data = jags.data, 
                                   working.directory = NULL, inits = NULL, parameters.to.save = jags.params, 
                                   model.file = paste(MODEL), n.burnin = 300, 
                                   n.thin = 10, n.iter = 1000, n.chains = n.chains)
      }
      Ldat$GLmean[i] <- median(jagsfitSLN$BUGSoutput$sims.list$GLmean.d)
      Ldat$GLmean.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d, 
                                     0.025)
      Ldat$GLmean.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$GLmean.d, 
                                     0.975)
      Ldat$SD[i] <- median(jagsfitSLN$BUGSoutput$sims.list$SD.d)
      Ldat$SD.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d, 
                                 0.025)
      Ldat$SD.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$SD.d, 
                                 0.975)
      Ldat$MK[i] <- median(jagsfitSLN$BUGSoutput$sims.list$MK.d)
      Ldat$MK.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d, 
                                 0.025)
      Ldat$MK.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$MK.d, 
                                 0.975)
      Ldat$FK[i] <- median(jagsfitSLN$BUGSoutput$sims.list$FK.d)
      Ldat$FK.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d, 
                                 0.025)
      Ldat$FK.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$FK.d, 
                                 0.975)
      FMi <- jagsfitSLN$BUGSoutput$sims.list$FK.d/jagsfitSLN$BUGSoutput$sims.list$MK.d
      Ldat$FM[i] <- median(FMi)
      Ldat$FM.lcl[i] <- quantile(FMi, 0.025)
      Ldat$FM.ucl[i] <- quantile(FMi, 0.975)
      ZKi <- jagsfitSLN$BUGSoutput$sims.list$MK.d + jagsfitSLN$BUGSoutput$sims.list$FK.d
      Ldat$ZK[i] <- median(ZKi)
      Ldat$ZK.lcl[i] <- quantile(ZKi, 0.025)
      Ldat$ZK.ucl[i] <- quantile(ZKi, 0.975)
      Ldat$r.Lopt[i] <- 3/(3 + Ldat$MK[i])
      Ldat$Linf[i] <- median((jagsfitSLN$BUGSoutput$sims.list$Linf.d))
      Ldat$Linf.lcl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d, 
                                   0.025)
      Ldat$Linf.ucl[i] <- quantile(jagsfitSLN$BUGSoutput$sims.list$Linf.d, 
                                   0.975)
    }
    BH.list <- BH(AllLength = unique(AllLength[AllYear == 
                                                 Year]), Linf = Ldat$Linf[i], MK = Ldat$MK[i], FK = Ldat$FK[i], 
                  GausSel = GausSel, selpar1 = ifelse(GausSel == T, 
                                                      Ldat$GLmean[i]/Ldat$Linf[i], Ldat$Lc[i]/Ldat$Linf[i]), 
                  selpar2 = ifelse(GausSel == T, Ldat$SD[i]/Ldat$Linf[i], 
                                   Ldat$r.alpha[i]))
    Ldat$BB0[i] <- as.numeric(BH.list[1])
    Ldat$YR[i] <- as.numeric(BH.list[2])
    rel.lcl <- sqrt(((Ldat$FM[i] - Ldat$FM.lcl[i])/Ldat$FM[i])^2 + 
                      ((Ldat$MK[i] - Ldat$MK.lcl[i])/Ldat$MK[i])^2 + ((Ldat$FK[i] - 
                                                                         Ldat$FK.lcl[i])/Ldat$FK[i])^2 + ((Ldat$Linf[i] - 
                                                                                                             Ldat$Linf.lcl[i])/Ldat$Linf[i])^2)
    rel.ucl <- sqrt(((Ldat$FM.ucl[i] - Ldat$FM[i])/Ldat$FM[i])^2 + 
                      ((Ldat$MK.ucl[i] - Ldat$MK[i])/Ldat$MK[i])^2 + ((Ldat$FK.ucl[i] - 
                                                                         Ldat$FK[i])/Ldat$FK[i])^2 + ((Ldat$Linf.ucl[i] - 
                                                                                                         Ldat$Linf[i])/Ldat$Linf[i])^2)
    Ldat$BB0.lcl[i] <- Ldat$BB0[i] - Ldat$BB0[i] * rel.lcl
    Ldat$BB0.ucl[i] <- Ldat$BB0[i] + Ldat$BB0[i] * rel.ucl
    Ldat$YR.lcl[i] <- Ldat$YR[i] - Ldat$YR[i] * rel.lcl
    Ldat$YR.ucl[i] <- Ldat$YR[i] + Ldat$YR[i] * rel.ucl
    Ldat$L95[i] <- Hmisc::wtd.quantile(x = L.y, weights = r.Freq.y, 
                                       probs = c(0.95))
    Ldat$perc.mat[i] <- ifelse(is.na(res$Lm50) == F, sum(r.Freq.y[L.y > 
                                                                    res$Lm50])/sum(r.Freq.y), NA)
    if (plot) {
      if (which(Years == Year) == 1) {
        if (is.na(mfrow)) {
          ncols <- ifelse(length(Years) <= 4, 2, 3)
          if (length(Years) == 1) 
            ncols <- 1
          opar <- par(mfrow = c(ceiling(length(Years)/3), 
                                ncols))
        }
        else {
          opar <- par(mfrow = mfrow)
        }
      }
      r.L.y <- L.y[L.y < Ldat$Linf[i]]/Ldat$Linf[i]
      r.Freq.y <- r.Freq.y[L.y < Ldat$Linf[i]]
      plotLBB.year(r.L.y = r.L.y, r.Freq.y = r.Freq.y, 
                   r.Lopt = Ldat$r.Lopt[i], SL1 = ifelse(GausSel == 
                                                           T, Ldat$GLmean[i], Ldat$Lc[i]), SL2 = ifelse(GausSel == 
                                                                                                          T, Ldat$SD[i], Ldat$r.alpha[i]), MK = Ldat$MK[i], 
                   FK = Ldat$FK[i], Linf = Ldat$Linf[i], Year = Year, 
                   GausSel = GausSel)
      if (which(Years == Year) == length(Years)) 
        par(opar)
    }
  }
  Linf.med <- median(Ldat$Linf)
  Linf.lcl <- median(Ldat$Linf.lcl)
  Linf.ucl <- median(Ldat$Linf.ucl)
  if (GausSel == F) {
    Lc.med <- median(Ldat$Lc)
    r.alpha.med <- median(Ldat$r.alpha)
  }
  else {
    GLmean.med <- median(Ldat$GLmean)
    SD.med <- median(Ldat$SD)
  }
  MK.med <- median(Ldat$MK)
  MK.lcl <- median(Ldat$MK.lcl)
  MK.ucl <- median(Ldat$MK.ucl)
  FK.med <- median(Ldat$FK)
  FK.lcl <- median(Ldat$FK.lcl)
  FK.ucl <- median(Ldat$FK.ucl)
  FM.med <- median(Ldat$FM)
  FM.lcl <- median(Ldat$FM.lcl)
  FM.ucl <- median(Ldat$FM.ucl)
  ZK.med <- median(Ldat$ZK)
  ZK.lcl <- median(Ldat$ZK.lcl)
  ZK.ucl <- median(Ldat$ZK.ucl)
  r.Lopt.med <- median(Ldat$r.Lopt)
  Lopt.med <- r.Lopt.med * Linf.med
  Lc_opt.med <- Linf.med * (2 + 3 * FM.med)/((1 + FM.med) * 
                                               (3 + MK.med))
  BB0.med <- median(Ldat$BB0)
  BB0.lcl <- median(Ldat$BB0.lcl)
  BB0.ucl <- median(Ldat$BB0.ucl)
  YR.med <- median(Ldat$YR)
  YR.lcl <- median(Ldat$YR.lcl)
  YR.ucl <- median(Ldat$YR.ucl)
  BFM1B0.list <- BH(AllLength = unique(AllLength), Linf = Linf.med, 
                    MK = MK.med, FK = MK.med, GausSel = GausSel, selpar1 = ifelse(GausSel == 
                                                                                    T, r.Lopt.med, 5/(2 * (3 + MK.med))), selpar2 = ifelse(GausSel == 
                                                                                                                                             T, SD.med/Linf.med, r.alpha.med))
  BFM1B0 <- as.numeric(BFM1B0.list[1])
  YRFM1 <- as.numeric(BFM1B0.list[2])
  cat("\n----------------------------------------------------------------------\n")
  cat("Results for", Species, ", stock", Stock, ",", StartYear, 
      "-", EndYear, ifelse(GausSel == T, ", Gaussian selection", 
                           ""), "\n")
  cat("-----------------------------------------------------------------------\n")
  cat("Linf prior =", Linf.pr, ", SD =", Linf.sd.pr, "(cm)", 
      ifelse(is.na(LinfUser) == TRUE, "", "(user-defined)"), 
      "\n")
  cat("Z/K prior  =", ZK.nls, ", SD =", ZK.nls.sd, ", M/K prior  =", 
      MK.pr, ", SD =", MK.sd.pr, ifelse(is.na(MKUser) == TRUE, 
                                        "", "(user-defined)"), "\n")
  if (GausSel == F) {
    cat("F/K prior  =", FK.pr, "(wide range with tau=4 in log-normal distribution)\n")
    cat("Lc prior   =", Lc.pr, ", SD =", Lc.sd.pr, "(cm)", 
        ifelse(is.na(LcUser) == TRUE, "", "(user-defined)"), 
        ", alpha prior=", r.alpha.pr, ", SD =", 0.1 * r.alpha.pr, 
        "\n\n")
  }
  cat("General reference points [median across years]: \n")
  cat("Linf               =", format(Linf.med, digits = 3), 
      paste("(", format(Linf.lcl, digits = 3), "-", format(Linf.ucl, 
                                                           digits = 3), ifelse(mmUser == F, ") cm", ") mm"), 
            sep = ""), "\n")
  cat("Lopt               =", format(Lopt.med, digits = 2), 
      paste(ifelse(mmUser == F, "cm,", "mm,"), "Lopt/Linf ="), 
      format(r.Lopt.med, digits = 2), "\n")
  cat("Lc_opt             =", format(Lc_opt.med, digits = 2), 
      paste(ifelse(mmUser == F, "cm,", "mm,"), "Lc_opt/Linf ="), 
      format(Lc_opt.med/Linf.med, digits = 2), "\n")
  cat("M/K                =", format(MK.med, digits = 3), paste("(", 
                                                                format(MK.lcl, digits = 3), "-", format(MK.ucl, digits = 3), 
                                                                ")", sep = ""), "\n")
  cat("F/K                =", format(FK.med, digits = 3), paste("(", 
                                                                format(FK.lcl, digits = 3), "-", format(FK.ucl, digits = 3), 
                                                                ")", sep = ""), "\n")
  cat("Z/K                =", format(ZK.med, digits = 3), paste("(", 
                                                                format(ZK.lcl, digits = 3), "-", format(ZK.ucl, digits = 3), 
                                                                ")", sep = ""), "\n")
  cat("F/M                =", format(FM.med, digits = 3), paste("(", 
                                                                format(FM.lcl, digits = 3), "-", format(FM.ucl, digits = 3), 
                                                                ")", sep = ""), "\n")
  cat(ifelse(GausSel == F, "B/B0 F=M Lc=Lc_opt =", "B/B0 F=M Lmean=Lopt="), 
      format(BFM1B0, digits = 3), "\n")
  cat("B/B0               =", format(BB0.med, digits = 3), 
      paste("(", format(BB0.lcl, digits = 3), "-", format(BB0.ucl, 
                                                          digits = 3), ")", sep = ""), "\n")
  cat(ifelse(GausSel == F, "Y/R' F=M Lc=Lc_opt =", "Y/R' F=M Lmean=Lopt="), 
      format(YRFM1, digits = 3), "\n")
  cat("Y/R'               =", format(YR.med, digits = 3), paste("(", 
                                                                format(YR.lcl, digits = 3), "-", format(YR.ucl, digits = 3), 
                                                                ")", sep = ""), "(linearly reduced if B/B0 < 0.25)\n\n")
  cat("Estimates for last year", EndYear, ":\n")
  last <- which(Ldat$Year == EndYear)
  if (GausSel == F) {
    cat("Lc         =", format(Ldat$Lc[last], digits = 3), 
        paste("(", format(Ldat$Lc.lcl[last], digits = 3), 
              "-", format(Ldat$Lc.ucl[last], digits = 3), ifelse(mmUser == 
                                                                   F, ") cm, Lc/Linf = ", ") mm, Lc/Linf = "), 
              format(Ldat$Lc[last]/Ldat$Linf[last], digits = 2), 
              " (", format(Ldat$Lc.lcl[last]/Ldat$Linf[last], 
                           digits = 3), "-", format(Ldat$Lc.ucl[last]/Ldat$Linf[last], 
                                                    digits = 3), ")", sep = ""), "\n")
    cat("alpha      =", format(Ldat$r.alpha[last], digits = 3), 
        "(", format(Ldat$r.alpha.lcl[last], digits = 3), 
        "-", format(Ldat$r.alpha.ucl[last], digits = 3), 
        ") \n")
    cat("Lmean/Lopt =", format(Ldat$Lmean[last]/(Ldat$r.Lopt[last] * 
                                                   Ldat$Linf[last]), digits = 2), ", Lc/Lc_opt =", format(Ldat$Lc[last]/Lc_opt.med, 
                                                                                                          digits = 2), ", L95th =", format(Ldat$L95[last], 
                                                                                                                                           digits = 3), ifelse(mmUser == F, "cm", "mm"), ", L95th/Linf =", 
        format(Ldat$L95[last]/Ldat$Linf[last], digits = 2), 
        ", \nLm50 =", format(res$Lm50, digits = 3), ifelse(mmUser == 
                                                             F, "cm", "mm"), ", Mature =", format(Ldat$perc.mat[last] * 
                                                                                                    100, digits = 2), "%\n")
  }
  else if (GausSel == T) {
    cat("GLmean/Linf=", format(Ldat$GLmean[last]/Ldat$Linf[last], 
                               digits = 2), ",SD/Linf =", format(Ldat$SD[last]/Ldat$Linf[last], 
                                                                 digits = 3), "\n")
    cat("GLmean     =", format(Ldat$GLmean[last], digits = 3), 
        ",SD =", format(Ldat$SD[last], digits = 3), "\n")
  }
  cat("F/K        =", format(Ldat$FK[last], digits = 2), "(", 
      format(Ldat$FK.lcl[last], digits = 3), "-", format(Ldat$FK.ucl[last], 
                                                         digits = 3), ")\n")
  cat("F/M        =", format(Ldat$FK[last]/Ldat$MK[last], digits = 2), 
      "(", format(Ldat$FM.lcl[last], digits = 3), "-", format(Ldat$FM.ucl[last], 
                                                              digits = 3), ")\n")
  cat("Z/K        =", format(Ldat$ZK[last], digits = 3), "(", 
      format(Ldat$ZK.lcl[last], digits = 3), "-", format(Ldat$ZK.ucl[last], 
                                                         digits = 3), ")\n")
  cat("Y/R'       =", format(Ldat$YR[last], digits = 2), "(", 
      format(Ldat$YR.lcl[last], digits = 3), "-", format(Ldat$YR.ucl[last], 
                                                         digits = 3), ") (linearly reduced if B/B0 < 0.25)\n")
  cat("B/B0       =", format(Ldat$BB0[last], digits = 2), "(", 
      format(Ldat$BB0.lcl[last], digits = 3), "-", format(Ldat$BB0.ucl[last], 
                                                          digits = 3), ")\n")
  cat("B/Bmsy     =", format(Ldat$BB0[last]/BFM1B0, digits = 2), 
      "(", format(Ldat$BB0.lcl[last]/BFM1B0, digits = 3), "-", 
      format(Ldat$BB0.ucl[last]/BFM1B0, digits = 3), ")\n")
  if (Comment != "" & !is.na(Comment)) 
    cat("Comment:", Comment, "\n")
  if (Ldat$MK[last] < 0 | Ldat$FK[i] < 0) 
    cat("Data unsuitable for LF analysis, negative mortality rates are impossible\n")
  if (Ldat$BB0[last] > 1.1) 
    cat("Data unsuitable for LF analysis, biomass exceeds carrying capacity\n")
  flush.console()
  priors <- t(data.frame(Linf = c(Linf.pr, Linf.sd.pr, LinfUser), 
                         ZK = c(ZK.nls, ZK.nls.sd, NA), MK = c(MK.pr, MK.sd.pr, 
                                                               MKUser), FK = c(FK.pr, NA, NA), Lc = c(Lc.pr, Lc.sd.pr, 
                                                                                                      LcUser), alpha = c(r.alpha.pr, 0.1 * r.alpha.pr, 
                                                                                                                         NA)))
  colnames(priors) <- c("prior", "SD", "user")
  medianRefLev <- t(data.frame(Linf = c(Linf.med, Linf.lcl, 
                                        Linf.ucl, mmUser), Lopt = c(Lopt.med, NA, NA, mmUser), 
                               Lc_opt = c(Lc_opt.med, NA, NA, mmUser), ZK = c(ZK.med, 
                                                                              ZK.lcl, ZK.ucl, NA), MK = c(MK.med, MK.lcl, MK.ucl, 
                                                                                                          NA), FK = c(FK.med, FK.lcl, FK.ucl, NA), FM = c(FM.med, 
                                                                                                                                                          FM.lcl, FM.ucl, NA), BB0 = c(BB0.med, BB0.lcl, BB0.ucl, 
                                                                                                                                                                                       NA), YR = c(YR.med, YR.lcl, YR.ucl, NA)))
  colnames(medianRefLev) <- c("median", "lo", "up", "user")
  last <- which(Ldat$Year == EndYear)
  lastRefLev <- t(data.frame(ZK = c(Ldat$ZK[last], Ldat$ZK.lcl[last], 
                                    Ldat$ZK.ucl[last], NA), MK = c(Ldat$MK[last], Ldat$MK.lcl[last], 
                                                                   Ldat$MK.ucl[last], NA), FK = c(Ldat$FK[last], Ldat$FK.lcl[last], 
                                                                                                  Ldat$FK.ucl[last], NA), FM = c(Ldat$FM[last], Ldat$FM.lcl[last], 
                                                                                                                                 Ldat$FM.ucl[last], NA), BB0 = c(Ldat$BB0[last], Ldat$BB0.lcl[last], 
                                                                                                                                                                 Ldat$BB0.ucl[last], NA), YR = c(Ldat$YR[last], Ldat$YR.lcl[last], 
                                                                                                                                                                                                 Ldat$YR.ucl[last], NA), BBmsy = c(Ldat$BB0[last]/BFM1B0, 
                                                                                                                                                                                                                                   Ldat$BB0.lcl[last]/BFM1B0, Ldat$BB0.ucl[last]/BFM1B0, 
                                                                                                                                                                                                                                   NA)))
  colnames(lastRefLev) <- c("median", "lo", "up", "mm")
  if (!GausSel) {
    lastRefLev <- rbind(lastRefLev, t(data.frame(Lc = c(Ldat$Lc[last], 
                                                        Ldat$Lc.lcl[last], Ldat$Lc.ucl[last], mmUser), LcLinf = c(Ldat$Lc[last]/Ldat$Linf[last], 
                                                                                                                  Ldat$Lc.lcl[last]/Ldat$Linf[last], Ldat$Lc.ucl[last]/Ldat$Linf[last], 
                                                                                                                  NA), alpha = c(Ldat$r.alpha[last], Ldat$r.alpha.lcl[last], 
                                                                                                                                 Ldat$r.alpha.ucl[last], NA), LmeanLopt = c(Ldat$Lmean[last]/(Ldat$r.Lopt[last] * 
                                                                                                                                                                                                Ldat$Linf[last]), NA, NA, NA), LcLc_opt = c(Ldat$Lc[last]/Lc_opt.med, 
                                                                                                                                                                                                                                            NA, NA, NA), L95th = c(Ldat$L95[last], NA, NA, mmUser), 
                                                 L95thLinf = c(Ldat$L95[last]/Ldat$Linf[last], NA, 
                                                               NA, NA), Lm50 = c(res$Lm50, NA, NA, mmUser), 
                                                 Mature = c(Ldat$perc.mat[last] * 100, NA, NA, NA))))
  }
  else if (GausSel) {
    lastRefLev <- rbind(lastRefLev, t(data.frame(GLmean = c(Ldat$GLmean[last], 
                                                            Ldat$SD[last], NA, NA), GLmeanLinf = c(Ldat$GLmean[last]/Ldat$Linf[last], 
                                                                                                   Ldat$SD[last]/Ldat$Linf[last], NA, NA))))
  }
  ret <- c(res, list(GausSel = GausSel, priors = priors, refLev = Ldat, 
                     medianRefLev = medianRefLev, lastRefLev = lastRefLev, 
                     LFall = LF.all))
  return(ret)
}
