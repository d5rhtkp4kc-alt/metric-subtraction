<!-- PROJECT LOGO -->

<br />

<h3 align="center">

Metric Subtraction for Non-Euclidean Random Objects with Applications to Residual Analysis

</h3>

</div>

<!-- ABOUT THE PROJECT -->

## Abstract

Residual analysis for non-Euclidean regression models is challenging due to the lack of a natural subtraction operation and a corresponding notion of conditionally centered residuals. We introduce a metric subtraction framework that yields residuals in the original response space, with their conditional Fr\'echet mean equal to a fixed reference object. We establish sufficient conditions for this conditional centering property, together with a separate condition ensuring distance preservation, providing a construction principle for suitable metric subtraction operators. The proposed framework accommodates several important classes of object spaces, including spheres, symmetric positive definite matrices, simplices, and Wasserstein spaces. Building on this framework, we develop a randomized test for detecting remaining conditional location effects and establish its asymptotic size and consistency against fixed alternatives. We further introduce two complementary diagnostic tools, signed distance plots and distance-profile-based boxplots, to visualize systematic residual patterns and identify atypical observations. For uncertainty quantification, we construct residual-based prediction sets and establish their asymptotic validity. Simulation studies and applications to New York City bicycle usage distributions and U.S. electricity generation compositions  demonstrate the practical utility of the proposed framework for model assessment, residual diagnostics, and prediction in non-Euclidean regression. 


### Code

1.  basic_function.R --- basic R functions for numerical studies. 
2.  visual_plot.R --- R functions to generate plots in methodology sections.
3.  simulation.R --- R functions for simulation.
4.  real_data.R --- R functions for real data.

### Data

1.  US_energy.csv --- US energy generation composition dataset downloaded from <https://www.eia.gov/electricity/data/browser/>. 
2.  energy_AR_residual.RData --- Preproccessed US energy generation composition data after removing trend, periodicity, AR components.
3.  energy_AR_residual_no_out.RData --- Preproccessed US energy generation composition data after removing outliers, trend, periodicity, AR components.
4.  NY_2019_trip_dist.RData --- Preprocessed New York City Citi bike sharing system data with original dataset downloaded from <https://citibikenyc.com/system-data>. 

