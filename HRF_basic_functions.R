

###############################################################################################
####   Basic functions needed for HRF_sim_est_inference.R
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     October-December 2016
###############################################################################################


gamma2 = function (t, TR, kernel_L, par) {
	#parameters of the double gamma function:
	#alpha1, alpha2: control shape
	#beta1, beta2: control scale
	#c: determines the ratio of the undershoot to the response
	A = par[1]
	alpha1 = par[2]
	alpha2 = par[3]
	beta1 = par[4]
	beta2 = par[5]
	c = par[6]
	N = length(t)
	h = A*((t^(alpha1-1)*beta1^alpha1*exp((-1)*beta1*t))/gamma(alpha1) - c*
	       (t^(alpha2-1)*beta2^alpha2*exp((-1)*beta2*t))/gamma(alpha2))
	h[(kernel_L+1):N] = 0
	return(h)
}

combine_h_stimulus = function(h, stimulus_T, N) {
	combined = rep(0,N)
	for (i in 1:length(stimulus_T)) {
		combined[stimulus_T[i]:N] = combined[stimulus_T[i]:N] + h[1:(N-stimulus_T[i]+1)]
	}
	return (combined)
}

gamma2_specified = function (x, kernel_L=kernel_L, stimulus_T=stimulus_T, N=N) {
	h = gamma2(t=x, TR=TR, kernel_L=kernel_L, par=gamma2_par)
	combined = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
	return(combined)
}

der1_specified = function (x, kernel_L, stimulus_T=stimulus_T, N=N) {
	return(grad(gamma2_specified, x=x, kernel_L=kernel_L, stimulus_T=stimulus_T, N=N, method="simple"))
}

L = function (x) {
	return(1/(1+exp(-x)))
}

IL = function (t, TR, kernel_L, par) {
	A1 = par[1]
	A2 = par[2]
	A3 = par[3]
	T1 = par[4]
	T2 = par[5]
	T3 = par[6]
	D1 = par[7]
	D2 = par[8]
	D3 = par[9]
	N = length(t)
	h = A1*L((t-T1)/D1) + A2*L((t-T2)/D2) + A3*L((t-T3)/D3)
	h[(kernel_L+1):N] = 0
	return(h)
}

FIR = function (N, par, kernel_L) {
	h = par
	h = c(h, rep(0, N-length(h)))
	h[(kernel_L+1):N] = 0
	return(h)
}

#par = A, alpha1, alpha2, beta1, beta2, c
cost_gamma2 = function(par) {
	t = seq(TR/2,N*TR-TR/2,by=TR)
	h = gamma2(t=t, TR=TR, kernel_L=kernel_L, par=par)
	fMRI_fit = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
	return (mean((fMRI - fMRI_fit)^2))
}

#par = A1, T1, T2, T3, D1, D2, D3
cost_IL = function(par) {
	A1 = par[1]
	T1 = par[2]
	T2 = par[3]
	T3 = par[4]
	D1 = par[5]
	D2 = par[6]
	D3 = par[7]
	#fitted response should begin at zero (for time point 0)
	A2 =  A1*(L((-T3)/D3)-L((-T1)/D1)) / (L((-T3)/D3)+L((-T2)/D2))
	#fitted response ends at magnitude zero
	A3 = -(A1+A2)
	t = seq(TR/2,N*TR-TR/2,by=TR)
	par_in = c(A1, A2, A3, T1, T2, T3, D1, D2, D3)
	h = IL(t=t, TR=TR, kernel_L=kernel_L, par=par_in)
	fMRI_fit = combine_h_stimulus(h=h, stimulus_T=stimulus_T, N=N)
	return (mean((fMRI - fMRI_fit)^2))
}

F_test = function(fMRI, fMRI_fit, p2) {
	RSS1 = sum((fMRI - mean(fMRI))^2)
	RSS2 = sum((fMRI - fMRI_fit)^2)
	F = ((RSS1-RSS2)/(p2-1)) / (RSS2/(length(fMRI)-p2))
	p_val = 1-pf(F, df1=p2-1, df2=length(fMRI)-p2)
	return (p_val)
}

