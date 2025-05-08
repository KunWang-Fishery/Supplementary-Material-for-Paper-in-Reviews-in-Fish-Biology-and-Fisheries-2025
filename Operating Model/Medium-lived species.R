rm(list=ls())
{library(LIME)
  library(ggplot2)
  library(dplyr)
}

##Medium-lived Species ====
{lh <- create_lh_list(vbk=0.21, 
                      linf=64, 
                      t0=-0.01,
                      lwa=2.45e-2,
                      lwb=2.79, 
                      M50=24,
                      M95=24,
                      maturity_input="length",
                      M=0.32,
                      S50=25, ## starting value
                      S95=25+0.0001, ## starting value
                      selex_input="length",
                      selex_type=c("logistic"),
                      CVlen=0.05,
                      binwidth=1,
                      h=0.7,
                      theta=10,
                      R0=1,
                      Fequil = 0.5,
                      SigmaC = 0,
                      SigmaI = 0,
                      SigmaR = 0,
                      SigmaF = 0,
                      nseasons=1,
                      nfleets=1)

true <- generate_data(modpath=NULL,
                      itervec=1, 
                      Fdynamics="None",# Oneway,Endogenous,None
                      Rdynamics="BevertonHolt",# Constant,BevertonHolt,
                      lh=lh,
                      Nyears=30,
                      Nyears_comp=30,
                      comp_sample=2000,
                      init_depl=1,
                      seed=1,
                      fleet_proportions=1)
I0<-true$I_ft[1]## use abundance index to reflect Biomass

}
## Equilibrium====
lh <- create_lh_list(vbk=0.21, 
                     linf=64, 
                     t0=-0.01,
                     lwa=2.45e-2,
                     lwb=2.79, 
                     M50=24,
                     M95=24,
                     maturity_input="length",
                     M=0.32,
                     S50=25, ## starting value
                     S95=25+0.0001, ## starting value
                     selex_input="length",
                     selex_type=c("logistic"),
                     CVlen=0.05,
                     binwidth=1,
                     h=0.7,
                     theta=10,
                     R0=1,
                     Fequil = 0.5,
                     SigmaC = 0,
                     SigmaI = 0,
                     SigmaR = 0,
                     SigmaF = 0,
                     nseasons=1,
                     nfleets=1)

true <- generate_data(modpath=NULL,
                      itervec=1, 
                      Fdynamics="Constant",# Oneway,Endogenous,None
                      Rdynamics="Constant",# Constant,BevertonHolt,
                      lh=lh,
                      Nyears=30,
                      Nyears_comp=30,
                      comp_sample=2000,
                      init_depl=0.25,
                      seed=1,
                      fleet_proportions=1)

## Length comp data input options
{LF_matrix <- true$LF[,,1]
  LF_data <- true$LF[20,,1]
  LF_array <- true$LF ## array with rows = years, columns = upper length bins, 3rd dimension = fleetst
  LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) 
  ## convert matrix, array, or list to data frame
  LF_df <- LFreq_df(LF=LF_list)
  
  # Filter only the last year
  last_year <- max(as.numeric(as.character(LF_df$Year)))
  LF_df_last <- LF_df %>% filter(Year == last_year)
  
  # Calculate lag-1 autocorrelation
  autocorr_last <- cor(LF_df_last$Length[1:(nrow(LF_df_last)-1)], LF_df_last$Length[2:nrow(LF_df_last)])
  
  # Calculate Effective Sample Size (ESS)
  N_last <- nrow(LF_df_last)
  ESS_last <- round(N_last / (1 + 2 * autocorr_last))
  
  # set sample size as ESS
  set.seed(1)
  LF_df <- LF_df %>%
    group_by(Year) %>%
    sample_n(size = ESS_last) %>%
    ungroup()
  ## length data only
  data_LF <- list("years"=1:true$Nyears, "LF"=LF_df)
  ## create model inputs with life history information and data
  ## outputs length data as array
  inputs_LC <- create_inputs(lh=lh, input_data=data_LF)
  #frequency table data for LBB
  freq_table <- LF_df %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    filter(Year == max(Year)) %>%
    mutate(bin = cut(Length,
                     breaks = seq(0, ceiling(max(Length)) + 0.5, by = 1),
                     include.lowest = TRUE, right = FALSE)) %>%
    count(bin, name = "count") %>%
    mutate(lengthClass = as.numeric(sub("\\[(.*),.*", "\\1", bin)) + 0.5) %>%
    right_join(
      data.frame(lengthClass = seq(0.5, ceiling(max(LF_df$Length)) + 0.5, by = 1)),
      by = "lengthClass"
    ) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    dplyr::select(lengthClass, count) %>%    
    arrange(lengthClass)
}
LF_eq<-ggplot(LF_df, aes(x = Length)) +
  geom_histogram(binwidth = 1, fill = "#a8e6a3", color = "#a8e6a3", boundary = 0) +
  facet_wrap(~ Year) +
  labs(title = "Length Frequency Distribution by Year", x = "Length (cm)", y = "Frequency") +
  theme_minimal()
LF_eq
## Non-equilibrium====
lh <- create_lh_list(vbk=0.21, 
                     linf=64, 
                     t0=-0.01,
                     lwa=2.45e-2,
                     lwb=2.79, 
                     M50=24,
                     M95=24,
                     maturity_input="length",
                     M=0.32,
                     S50=25, ## starting value
                     S95=25+0.0001, ## starting value
                     selex_input="length",
                     selex_type=c("logistic"),
                     CVlen=0.05,
                     binwidth=1,
                     h=0.7,
                     theta=10,
                     R0=1,
                     Fequil = 0.5,
                     SigmaC = 0,
                     SigmaI = 0,
                     SigmaR = 0.05,
                     SigmaF = 0.05,
                     nseasons=1,
                     nfleets=1)

# generate simulated data
true <- generate_data(modpath=NULL,
                      itervec=1, 
                      Fdynamics="Oneway",# Oneway,Endogenous,None
                      Rdynamics="BevertonHolt",# Constant,BevertonHolt,
                      lh=lh,
                      Nyears=30,
                      Nyears_comp=30,
                      comp_sample=2000,
                      init_depl=0.45,
                      seed=1,
                      fleet_proportions=1)
## Length comp data input options
{LF_matrix <- true$LF[,,1]
  LF_data <- true$LF[20,,1]
  LF_array <- true$LF ## array with rows = years, columns = upper length bins, 3rd dimension = fleetst
  LF_list <- lapply(1:lh$nfleets, function(x) true$LF[,,x]) 
  ## convert matrix, array, or list to data frame
  LF_df <- LFreq_df(LF=LF_list)
  
  # Filter only the last year
  last_year <- max(as.numeric(as.character(LF_df$Year)))
  LF_df_last <- LF_df %>% filter(Year == last_year)
  
  # Calculate lag-1 autocorrelation
  autocorr_last <- cor(LF_df_last$Length[1:(nrow(LF_df_last)-1)], LF_df_last$Length[2:nrow(LF_df_last)])
  
  # Calculate Effective Sample Size (ESS)
  N_last <- nrow(LF_df_last)
  ESS_last <- round(N_last / (1 + 2 * autocorr_last))
  
  # set sample size as ESS
  set.seed(1)
  LF_df <- LF_df %>%
    group_by(Year) %>%
    sample_n(size = ESS_last) %>%
    ungroup()
  ## length data only
  data_LF <- list("years"=1:true$Nyears, "LF"=LF_df)
  ## create model inputs with life history information and data
  ## outputs length data as array
  inputs_LC <- create_inputs(lh=lh, input_data=data_LF)
  #frequency table data for LBB
  freq_table <- LF_df %>%
    mutate(Year = as.numeric(as.character(Year))) %>%
    filter(Year == max(Year)) %>%
    mutate(bin = cut(Length,
                     breaks = seq(0, ceiling(max(Length)) + 0.5, by = 1),
                     include.lowest = TRUE, right = FALSE)) %>%
    count(bin, name = "count") %>%
    mutate(lengthClass = as.numeric(sub("\\[(.*),.*", "\\1", bin)) + 0.5) %>%
    right_join(
      data.frame(lengthClass = seq(0.5, ceiling(max(LF_df$Length)) + 0.5, by = 1)),
      by = "lengthClass"
    ) %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    dplyr::select(lengthClass, count) %>%    
    arrange(lengthClass)
}
LF_neq<-ggplot(LF_df, aes(x = Length)) +
  geom_histogram(binwidth = 1, fill = "#a8e6a3", color = "#a8e6a3", boundary = 0) +
  facet_wrap(~ Year) +
  labs(title = "Length Frequency Distribution by Year", x = "Length (cm)", y = "Frequency") +
  theme_minimal()
LF_neq
# comparison of different scenario
LF_eq
LF_neq
