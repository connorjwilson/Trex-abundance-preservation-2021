# Provides input values for Monte Carlo simulations and generates all output distributions, including those shown in Fig. 2 and the standing and total biomasses using the procedures detailed in Marshall et al. (2021), "Absolute abundance and preservation rate of 'Tyrannosaurus rex'. Requires the file trex_functions.R, and places the results into trex_output.RData.  

# -----------------------------------
# --- Asymptotic body mass (kg) for growth curve, Fig. 2c.
# -----------------------------------
AsymBMmean <- 7090
AsymBMsdev <- 1012.755

# ------------------------------
# ---- Intercept for population density (ind/km^2) versus ecological body mass using Damuth's rule, eq. 1, Fig. 2a.
# ------------------------------
logA_mid <- 2.991
logA_sd <- 0.6060698

# -------------------------------
# ---- Age at sexual maturity (yrs), Fig. 2j.
# -------------------------------
SMmean <- 15.5 
SMsdev <- 0.7653

# ----------------------------------------
# ---- Temporal range (millions of years), Fig. 2h.
# ----------------------------------------
TRmin <- 1.2
TRmax <- 3.6

# --------------------------------------
# --- Geographic range (millions of km^2), Fig. 2f.
# --------------------------------------
GRmean <- 2.3
GRsdev <- 0.449

source("trex_functions.R")

n <- 1000000 # Number of Monte Carlo simulations to perform

# Run Monte Carlo simulations
timestamp()
list_trex2 <- trex_ever_lived2(n, maxAge = 28, AsymBMmean, AsymBMsdev, damuth_slope = -3/4, logA_mid, logA_sd, SMmean, SMsdev, TRmin, TRmax, GRmean, GRsdev, agesResolution = 1) # note  that it is possible to change the slope for the Damuth equation here.
timestamp()

trex_ranges <- t(sapply(list_trex2, quantile, c(0.025, 0.5,0.975), na.rm = T))

trex_ranges

standing_biomass <- list_trex2$standing_number * list_trex2$mean_body_mass
quantile(standing_biomass, c(0.025, 0.5,0.975))
total_biomass <- list_trex2$total_number * list_trex2$mean_body_mass
quantile(total_biomass, c(0.025, 0.5,0.975))