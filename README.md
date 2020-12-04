# Multicarving
R Software for inference using multicarving

The code for this thesis consists of 8 files:
    • 3 files for methodology: hdi_adjustments.R; optimal_inference.R; sample_from_truncated.R;
    • 4 files for simulations: all_multi_fix_fwer.R; binomial_all_multi_fix_fwer.R; group_example_mc.R; multi_ci.R
    • 1 file for error handling: tryCatch-W-E.R (not my own work)

hdi_adjustments.R: Functions adapted from the ‘hdi’-package in order to change multisplitting routines to multicarving routines

optimal_inference.R: Routines to calculate carving p-values. Does the necessary preprocessing as well as the evaluation of the MCMC samples

sample_from_truncated.R: Hit-and-run MCMC sampler for truncated, multivariate, white Gaussian 

all_multi_fix_fwer.R: Executes simulations as described in Section 4.1., default is the Toeplitz design with inference in the selected model. For other scenarios, slight adjustments have to be done. Especially, undo the comments under # riboflavin if doing simulations for the riboflavin data. Do usei <- TRUE to estimate σ on a per model base, use_lambda.min = TRUE to do cross-validated lasso using λmin, selected = FALSE for inference in the saturated model. Further changes are possible based on user’s needs.

binomial_all_multi_fix_fwer.R: Executes simulations as described in Section 4.4. Changes are possible based on user’s needs.

group_example_mc.R: Executes simulations as described in Section 4.3, default is the dense alternative. Adjustments for the sparse alternative are described through commented code. Further changes are possible based on user’s needs.

multi_ci.R: Executes simulations as described in section 4.2. Changes are possible based on user’s needs.

All simulation files store the data needed for the analysis. The save location has to be adjusted based on user’s folder structure. Due to slight changes in the code that decrease redundancy and increase efficiency, the obtained results will not be exactly the same as shown in the thesis. However, as the methodology is not changed per se, the results should stay fairly similar.

tryCatch-W-E.R: Tries execution of the statement and stores potential errors and warnings. Provided to me by Claude Renaux.
