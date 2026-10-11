# dG-ComBat for Multi-Site Data Harmonization

This is an implementation of the dG-ComBat algorithm for harmonizing multi-site neuroimage data. The method uses a distributed site-central workflow to estimate and remove site-related location, scale, and covariance effects while preserving biological covariates of interest.

# Codes

`dg_combat_function.R` contains the main site-side and central-side functions for the three-round dG-ComBat harmonization procedure, including data standardization, empirical Bayes correction, and covariance alignment.

`help_function.R` contains helper functions for estimating the shared graph structure, calculating hyperparameters, and iteratively updating the ComBat location and scale parameters.

`data_generation_function.R` contains functions for generating multiblock simulated data with latent covariance structure, feature-level batch effects, and inverse-gamma noise variances.

`single_run_demo.R` is a single-run example that generates a four-site simulated dataset, preserves the Alzheimer's disease (AD) effect through the model matrix, and applies the complete dG-ComBat harmonization workflow.

# Examples

We include a demonstration of dG-ComBat on simulated multi-site data. `single_run_demo.R` defines the simulation settings, generates the data, saves the generated dataset as `demo_data.rds`, and performs the three communication rounds between the individual sites and the central server. The final harmonized data are stored in `Y_new` in the R session.

