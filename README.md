# Multicarving
R Software for inference using multicarving as presented in https://arxiv.org/abs/2006.04613

The code for this project consists of 8 files:
 - In /inference/  
   - 3 files for methodology: hdi_adjustments.R; carving.R; sample_from_truncated.R;   
   - 1 file for error handling: tryCatch-W-E.R (not my own work);  
 - In /simulation_setups/  
   - 4 files for simulations: carving_simulation.R; carving_binomial_simulation.R; carving_group_simulation; carving_ci_simulation.R


hdi_adjustments.R: Functions adapted from the ‘hdi’-package in order to change multisplitting routines to multicarving routines
The user will mainly interact with the following four functions:
   - multi.carve: Function to execute the whole multicarving process, i.e. selecting a model and infering on each split as well as calculate multicarving p-values. Must at         least provide predictor matrix (x) and response vector (y).
   - carve100: Function to execute the whole process of pure post-selection inference, i.e. selecting a model and calculating p-values using all data for selection. Must at       least provide predictor matrix (x) and response vector (y).
   - multi.carve.group: Function to execute the whole multicarving process for groups, i.e. selecting a model and performing groupwise inference on each split as well as           calculate multicarving p-values. Must at least provide predictor matrix (x) and response vector (y) and a list of vectors for the groups that shall be tested.
   - multi.carve.ci.saturated: Function to determine multicarving confidence intervals using the saturated viewpoint, i.e. selecting a model and infering on each split as         well as calculate multicarving p-values and confindence intervals afterwards. Must at least provide predictor matrix (x) and response vector (y).

carving.R: Routines to calculate carving p-values. Does the necessary preprocessing as well as the evaluation of the MCMC samples. 

sample_from_truncated.R: Hit-and-run MCMC sampler for truncated, multivariate, white Gaussian.

tryCatch-W-E.R: Tries execution of the statement and stores potential errors and warnings. Provided to me by Claude Renaux.

carving_simulation.R: Executes simulations as described in Section 4.1., default is the Toeplitz design with inference in the selected model. For other scenarios, slight adjustments have to be done. Especially, undo the comments under # riboflavin if doing simulations for the riboflavin data. Do sigma.estimator = "modwise" to estimate σ on a per model base, use_lambda.min = TRUE (in args.model.selector) to do cross-validated lasso using λmin, selected = FALSE (in args.lasso.inference) for inference in the saturated model. Further changes are possible based on user’s needs.

carving_binomial_simulation.R: Executes simulations as described in Section 4.4. Changes are possible based on user’s needs.

carving__group_simulation.R: Executes simulations as described in Section 4.3, default is the dense alternative. Adjustments for the sparse alternative are described through commented code. Further changes are possible based on user’s needs.

carving_ci_simulation.R: Executes simulations as described in section 4.2. Changes are possible based on user’s needs.

All simulation files store the data needed for the analysis. Due to recent changes in the code, the obtained results will not be exactly the same as shown in the paper. However, as the methodology is not changed per se, the results should stay fairly similar.


