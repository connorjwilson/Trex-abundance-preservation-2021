# Provides the functions needed by trex_monte_carlo.R to run the Monte Carlo simulations to generate the output distributions (including those shown in Fig. 2 and the standing and total biomasses) using the procedures described in Marshall et al. (2021), "Absolute abundance and preservation rate of 'Tyrannosaurus rex'. 
#--------------------------------------
#--- Returns the vector of body masses for the ages in 't'.
#--- For this paper 't' is the vector of T. rex ages from the age at sexual maturity to maximum age (28 years).
#--- Body mass computed using Erickson et al's (2004) logistic growth equation, see Fig. 2c, eq. S9.
#--- 'M_max' is the asymptotic body mass.
#--------------------------------------
logistic <- function(M_max,b = 0.55,c = 16.2,d = 5,t){
  (M_max/(1+exp(-b*(t-c))))+d
}

#--------------------------------------
#--- Returns lx values, the proportion of individuals surviving to each of the ages in 't'.
#--- For this paper 't' is the vector of T. rex ages from the age at sexual maturity to age at death (28 years).
#--- Uses a Gompertz equation (see Erickson et al 2006) with parameters 'a' and 'g', see Fig. 2d, eq. S19.
#--------------------------------------
surv <- function(t, a, g){
  exp((a/g)*(1-exp(g*t)))
}

#--------------------------------------
#--- Returns the proportion of individuals in each cohort beginning at sexual maturity, which we call the weighted survivorship, see Eq. S10.
#--- Needed to calculate ecological body mass, eqs. S11, S12.
#--- Arguments are the same as for function 'surv' above.
#--------------------------------------
w_surv <- function(t, a, g){surv(t, a, g)/sum(surv(t, a, g))}

#--------------------------------------
#--- Computes the set of Gompertz equations to assess the range of plausible l_x values by finding the Gompertz parameters 'a' and 'g' from the bootstrapped  T. rex ages from Erickson et al 2006 (see also our Table S3).
#--- 'n' is the number of Monte Carlo simulations.
#--------------------------------------
bootstrap_survivorship <- function(n, obsAges = c(2,6,8,9,11,14,14,15,16,17,18,18,18,18,18,19,21,21,21,22,22,22,22,22,22,23,23,24,24,28)){
  aEstimates <- vector()
  gEstimates <- vector()
  for(i in 1:n){
    while(TRUE){
      sampledAges <- sample(obsAges, size = length(obsAges), replace = T)
      propSurvSamp <- vector()
      for(j in 1:length(sampledAges)){
        propSurvSamp[j] <- 1-sum(sampledAges<sampledAges[j])/length(sampledAges)
      }
      
      # Fits Gompertz equation to sampled ages
      gompertzFit <- try(nls(propSurvSamp ~ exp( (a/g)*(1-exp(g*sampledAges))), start = list(a = 0.002, g = 0.2214)), silent = TRUE)
      if(class(gompertzFit) != "try-error"){
        aEstimates[i] <- summary(gompertzFit)$parameters[1,1]
        gEstimates[i] <- summary(gompertzFit)$parameters[2,1]
        break()
      }
    }
  }
  return(cbind(aEstimates, gEstimates))
}

#--------------------------------------
#--- Performs 'n' Monte Carlo simulations to calculate the output values (see Fig. 2).
#--- Most input values are given and can be changed in 'Trex_monte_carlo.R'.
#--- Several variables require a discretization of time, for example, the age at sexual maturity, the cohort body massess, and survivorship values.  The temporal resolution of the discretization is named 'agesResolution' and is set to 1 year. 
#--------------------------------------
trex_ever_lived2 <- function(n, AsymBMmean, AsymBMsdev, damuth_slope = -0.75, logA_mid, logA_sd, SMmean, SMsdev, TRmin, TRmax, GRmean, GRsdev,  agesResolution = 1, maxAge = 28){

  #--------------------------------------
  #--- Age of sexual maturity (years) rounded to the nearest year, Fig. 2j.
  #--------------------------------------
  sex_mat <- round(rnorm(n, SMmean, SMsdev), digits = 1)
  
  #--------------------------------------
  #--- Vector of ages from age at sexual maturity to 28 yrs.
  #--- Generates the list of 'n' vectors, which are used to calculate the ecological body mass (eq. S12) and generation time (eq. S24).
  #--------------------------------------
  max_age <- maxAge
  ages <- lapply(sex_mat, seq, to = max_age, by = agesResolution)
  
  #--------------------------------------
  #--- Asymptotic body mass (kg), see inset in Fig. 2c.
  #--------------------------------------
  asymp_mass <- rnorm(n, AsymBMmean, AsymBMsdev)
  
  #--------------------------------------
  #--- Bootstrapped Gompertz survivorship equations, Fig. 2d.
  #--------------------------------------
  AxG <- bootstrap_survivorship(n)
  
  #--------------------------------------
  #--- Ecological Body Mass computation (grams), Fig. 2b, eq. S12.
  #--------------------------------------
  lx <- list()
  mean_body_mass <- vector()
  
  for(i in 1:n){
    
    foo <- w_surv(t = ages[[i]], a = AxG[i,1], g = AxG[i,2])
    mean_body_mass[i] <- sum(logistic(asymp_mass[i], t = ages[[i]])*foo)*1000
  }
  
  #--------------------------------------
  #--- Computes lx values from bootstrapped Gompertz survivorship equations.
  #--------------------------------------
  for(i in 1:n){
      
    lx[[i]] <- surv(t = ages[[i]], a = AxG[i,1], g = AxG[i,2])
    
  }
  
  #--------------------------------------
  #--- Computes generation time (yrs), eqs. S23 and S24.
  #--------------------------------------
  b_i <- sapply(lx, function(x){1/sum(x)}) 
  
  i_bi_li <- Map('*', Map('*', ages, lx), b_i)
  
  gen_time <- sapply(i_bi_li, sum)

  #--------------------------------------
  #--- Calculates mean population density (ind/km^2) using Damuth's equation, Fig. 2a.
  #--------------------------------------
  pop_density <- (10**(rnorm(n, logA_mid, logA_sd) + damuth_slope * log(mean_body_mass, base = 10)))
  
  #--------------------------------------
  #--- Temporal range (millions of years), Fig. 2h.
  #--------------------------------------
  temporal_range <- runif(n, TRmin, TRmax)
  
  #--------------------------------------
  #--- Number of generations (years), Fig. 2l.
  #--------------------------------------
  num_generations <- (temporal_range * 1e+06)/gen_time
  
  #--------------------------------------
  #--- Geographic range (millions of km^2), Fig. 2f.
  #--------------------------------------
  geographic_range <- rnorm(n, GRmean, GRsdev)
  
  #--------------------------------------
  #--- Standing number of individuals, Fig. 2g.
  #--------------------------------------
  standing_number <- geographic_range * 1e+06 * pop_density
  
  #--------------------------------------
  #--- Total number of individuals, Fig. 2m, eq. S2.
  #--------------------------------------
  total_number <- standing_number * num_generations
  
  #--------------------------------------
  #--- Returns output for the variables given in Fig. 2.
  #--------------------------------------
  return(list(total_number = total_number, standing_number = standing_number, geographic_range = geographic_range, num_generations = num_generations, temporal_range = temporal_range, gen_time = gen_time, pop_density = pop_density, mean_body_mass = mean_body_mass))
}
