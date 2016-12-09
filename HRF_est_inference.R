

###############################################################################################
####   Estimation and inference of the HRF model given a simulated fMRI time series.
####   Double gamma model ('canonical'), the former with 1 and 2 first derivatives,
####   the nonlinear double gamma model, the inverse logit model (IL) and finite impulse
####   response (FIR) considered. The function checks whether response detected using
####   a given HRF model.
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     October-December 2016
###############################################################################################


HRF_est_inference = function (fMRI, HRF, t, TR, stimulus_T, kernel_L, par_max, max.time=10) {
	#making the fMRI time series, stimulus_T, kernel_L and TR global variables, needed for the cost function in SA
	fMRI <<- fMRI
	stimulus_T <<- stimulus_T
	kernel_L <<- kernel_L
	TR <<- TR
	N = length(fMRI)
	par = rep(NA,par_max)
	if (HRF=="gamma2") {
		#default SPM options
		h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=gamma2_par)
		h_combined = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		R = lm(fMRI~h_combined)
		fMRI_fit = R$fitted.values
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=2)
		par[1:length(R$coef)] = R$coef
		#pre-whitening the time-series
		#resid = ar(R$resid, order=1)$resid
		#the first element NA
		#resid[1] = 0
		#fMRI = fMRI_fit + resid
		#R = lm(fMRI~h_combined)
		#f = summary(R)$fstatistic
		#p-value from the F-test
		#p_val = as.numeric(1 - pf(f[1], f[2], f[3]))
	} else if (HRF=="gamma2_D") {
		#default SPM options
		h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=gamma2_par)
		h_combined = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		#temporal derivative included
		der1 = grad(gamma2_specified, x=t, kernel_L=kernel_L, stimulus_T=stimulus_T, N=N, method="simple")
		R = lm(fMRI~h_combined+der1)
		fMRI_fit = R$fitted.values
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=3)
		par[1:length(R$coef)] = R$coef
	} else if (HRF=="gamma2_DD") {
		#default SPM options
		h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=gamma2_par)
		h_combined = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		#temporal and dispersion derivatives included
		der1 = grad(gamma2_specified, x=t, kernel_L=kernel_L, stimulus_T=stimulus_T, N=N, method="simple")
		der2 = grad(der1_specified, x=t, kernel_L=kernel_L, stimulus_T=stimulus_T, N=N, method="simple")
		R = lm(fMRI~h_combined+der1+der2)
		fMRI_fit = R$fitted.values
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=4)
		par[1:length(R$coef)] = R$coef
	} else if (HRF=="gamma2_nonlinear") {
		#no derivatives included
		#max.time given in seconds
		G2NL_par_est = GenSA(lower=gamma2_par_lower, upper=gamma2_par_upper, fn=cost_gamma2, control=list(max.time=max.time))$par
		h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=G2NL_par_est)
		fMRI_fit = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=6)
		par[1:length(G2NL_par_est)] = G2NL_par_est
	} else if (HRF=="IL") {
		#Lindquist 2007 codes: LB = [0.05, 1, 0, 0.05, 5, 0, 10] UB = [10, 15, 5, 10, 15, 5, 30]
		#max.time in seconds
		par_GenSA = GenSA(lower=IL_par_lower, upper=IL_par_upper, fn=cost_IL, control=list(max.time=max.time))$par		
		A1 = par_GenSA[1]
		T1 = par_GenSA[2]
		T2 = par_GenSA[3]
		T3 = par_GenSA[4]
		D1 = par_GenSA[5]
		D2 = par_GenSA[6]
		D3 = par_GenSA[7]
		#fitted response should begin at zero (for time point 0)
		A2 =  A1*(L((-T3)/D3)-L((-T1)/D1)) / (L((-T3)/D3)+L((-T2)/D2))
		#fitted response ends at magnitude zero
		A3 = -(A1+A2)
		IL_par_est = c(A1, A2, A3, T1, T2, T3, D1, D2, D3)
		h = IL(t=t, TR=TR, kernel_L=kernel_L, par=IL_par_est)
		fMRI_fit = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=7)
		par[1:length(IL_par_est)] = IL_par_est
	} else if (HRF=="FIR") {
		#the number of time points -> the number of parameters
		tp_no = length(FIR_par)
		FIR_par_est_A = array(NA, dim=c(tp_no,length(stimulus_T)))
		for (i in 1:length(stimulus_T)) {
			FIR_par_est_A[1:tp_no,i] = fMRI[stimulus_T[i]:(stimulus_T[i]+(tp_no-1))]
		}
		FIR_par_est = rep(NA, tp_no)
		for (i in 1:tp_no) {
			FIR_par_est[i] = sum(FIR_par_est_A[i,], na.rm=T) / sum(!is.na(FIR_par_est_A[i,]))
		}
		h = FIR(N=N, kernel_L=kernel_L, par=FIR_par_est)
		fMRI_fit = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
		p_val = F_test(fMRI=fMRI, fMRI_fit=fMRI_fit, p2=tp_no)
		par[1:length(FIR_par_est)] = FIR_par_est
	} else if (HRF=="sFIR") {
		###################################################################
		##  NOT READY YET
		##  confusing how to choose the smoothing parameter
	} else {
		return("the given HRF not implemented (yet?)")
	}
	#the p-value and the fMRI fit returned
	return(list(p_val=p_val, fMRI_fit=fMRI_fit, par=par))
}

