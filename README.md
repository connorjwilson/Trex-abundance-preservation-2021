The R code used for the Monte Carlo simulations and the output data from Marshall et al. (2021) "Absolute abundance and preservation rate of 'Tyrannosaurus rex'”.

Code written by Connor J. Wilson and Daniel V. de Latorre.


Files:

**trex_monte_carlo.R** – Assigns input values (mean and standard deviation for the normal distributions, total range for the uniform distribution) and calls functions from **trex_functions.R** to perform the Monte Carlo simulations that compute the distributions of output values.  

**trex_functions.R** – Defines the R functions needed to perform the analyses.

**trex_output.RData** – Contains the output from the million Monte Carlo simulations reported in Marshall et al. (2021).  


Running the code:

Open the **trex_monte_carlo.R**. The file **trex_functions.R** must be in the same working directory as **trex_monte_carlo.R**.  The output will be in **trex_output.RData**.


Experimenting with the data:

New input values can be changed in **trex_monte_carlo.R**, as can the number of Monte Carlo simulations, the slope for Damuth’s Law (damuth_slope) and the maximum age (maxAge).  The parameter values for the logistic growth curve and the ages assigned to the ‘T. rex’ specimens for the survivorship analyses can be changed in **trex_functions.R**.


The content of the code is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
