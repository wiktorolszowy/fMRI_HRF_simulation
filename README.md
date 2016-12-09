Written by: Wiktor Olszowy, University of Cambridge
Contact: wo222@cam.ac.uk
Created: September-December 2016

# HRF_simulation_for_fMRI_experiments
The fMRI time series are simulated given a set of options the user can specify. Afterwards, the models are estimated and significance of the stimulus evoked response is assessed. There are 6 HRF models available: (1) canonical double-gamma, (2-3) extension by its derivatives, (4) nonlinear double gamma, (5) inverse logit and (6) finite impulse response.


(1) HRF_sim_est_inference.R
is the interface

(2) HRF_est_inference.R
estimates the model parameters and performs inference

(3) HRF_sim.R
simulates fMRI time series given a set of options: repetition time (TR), experiment paradigm/design, HRF model, CNR, AR1

(4) HRF_basic_functions.R
auxiliary functions for the above codes

(5) HRF_par_dist.R
makes figures showing the distributions of the parameter estimates

(6) HRF_FSL_gamma2_FIR.R
performs analysis in FSL for both the canonical model with the 1st derivative and for the finite impulse response model


You need the following packages in R: AnalyzeFMRI, oro.nifti, parallel, reshape2, ggplot2, GenSA, numDeriv, latticeExtra, gridExtra

The codes were tested in R 3.2.2 on Linux.
