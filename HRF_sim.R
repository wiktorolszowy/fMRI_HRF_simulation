

###############################################################################################
####   Simulation of fMRI time series: double-gamma, inverse logit, FIR. Parameters of the
####   models are defined in HRF_sim_est_inference.R
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     October-December 2016
###############################################################################################


HRF_sim = function (HRF, t, TR, CNR, AR1, stimulus_T, kernel_L) {
	N = length(t)
	#very low TR also used, to know how the true signal looks like exactly
	#e.g. for some TRs the simulated peak can be much lower than the peak from the model
	TR_fine = 0.1
	t_fine = seq(TR_fine/2,N*TR-TR_fine/2,by=TR_fine)
	if (HRF=="gamma2") {
		h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=gamma2_par)
		h_fine = gamma2(t=t_fine, TR=TR_fine, kernel_L=kernel_L*(TR/TR_fine), par=gamma2_par)
	} else if (HRF=="IL") {
		h = IL(t=t, TR=TR, kernel_L=kernel_L, par=IL_par)
		h_fine = IL(t=t_fine, TR=TR_fine, kernel_L=kernel_L*(TR/TR_fine), par=IL_par)
	} else if (HRF=="FIR") {
		#arbitrary choice
		h = c(FIR_par, rep(0,N-length(FIR_par)))
		#for FIR each time point after the stimulus (up to length(FIR_par)) corresponds to one df/par
		t_fine = t
		h_fine = h
		TR_fine = TR
	} else {
		return("the given HRF not implemented")
	}
	combined = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
	combined_fine = combine_h_stimulus(h=h_fine, stimulus_T=round(stimulus_T*(TR/TR_fine)), N=length(t_fine))
	#combined and combined_fine are not perfectly aligned to each other, there is some shift due to the time location of the stimulus
	#combined_fine is needed only to calculate the sd of the true signal, so the time shift is no problem (unless extremely short time series)!
	fMRI_sim = arima.sim(n=N, model=list(ar=AR1))
	fMRI_sim = (fMRI_sim/sd(fMRI_sim)) * (sd(combined_fine)/CNR)
	fMRI_sim = fMRI_sim + combined
	return(fMRI_sim)
}

