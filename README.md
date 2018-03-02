| Written by: | Wiktor Olszowy, Department of Clinical Neurosciences, University of Cambridge     |
| ----------- | --------------------------------------------------------------------------------- |
| When:       | September 2016 - December 2016                                                    |
| Purpose:    | Simulating fMRI signal for different hemodynamic response functions (HRFs)        |
| Contact:    | wo222@cam.ac.uk                                                                   |

The fMRI time series are simulated given a set of options the user can specify. Afterwards, the models are estimated and significance of the stimulus evoked response is assessed. There are 6 HRF models available:

1. canonical double gamma
2. canonical double gamma with first derivative
3. canonical double gamma with first and second derivatives
4. nonlinear double gamma
5. inverse logit
6. finite impulse response

Software
==============

For the R computations I used the following packages (available from CRAN):

- AnalyzeFMRI
- GenSA
- ggplot2
- gridExtra
- latticeExtra
- numDeriv
- oro.nifti
- parallel
- reshape2

The codes were tested in R 3.2.2 on Linux.

Repository contents
==============

- `HRF_sim_est_inference.R`

  is the interface
- `HRF_est_inference.R`

  estimates the model parameters and performs inference
- `HRF_sim.R`

  simulates fMRI time series given a set of options: repetition time (TR), experimental paradigm/design, HRF model, CNR, AR1
- `HRF_basic_functions.R`

  auxiliary functions for the above codes
- `HRF_par_dist.R`

  makes figures showing the distributions of the parameter estimates
- `HRF_FSL_gamma2_FIR.R`
 
  performs analysis in FSL for both the canonical model with the 1st derivative and for the finite impulse response model
