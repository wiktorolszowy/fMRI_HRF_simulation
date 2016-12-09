

###############################################################################################
####   Figures showing distribution of the estimated HRF parameters. Run after
####   HRF_sim_est_inference.R. Red lines indicate the true values and green lines
####   indicate the search region boundaries.
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     October-December 2016
###############################################################################################


#ONLY ONE COMBINATION (WITH STRONG SIGNAL): LESS VARIABLE ESTIMATES
TR = "2"; CNR = paste(max(CNR_sample)); AR1 = "0.5"
pdf("par_dist_gamma2_NL.pdf")
	par(mfrow=c(4,2))
	hist(P_VAL_PARAMS["gamma2_nonlinear",TR,CNR,AR1,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	for (i in 1:length(gamma2_par)) {
		hist(P_VAL_PARAMS["gamma2_nonlinear",TR,CNR,AR1,i+1,], main=paste("gamma2 par", i))
		abline(v=gamma2_par[i], col="red")
		abline(v=c(gamma2_par_lower[i], gamma2_par_upper[i]), col="green")
	}
dev.off()
pdf("par_dist_IL.pdf")
	par(mfrow=c(4,2))
	hist(P_VAL_PARAMS["IL",TR,CNR,AR1,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	hist(P_VAL_PARAMS["IL",TR,CNR,AR1,2,], main=paste("IL par 1"))
	abline(v=IL_par[1], col="red")
	abline(v=c(IL_par_lower[1], IL_par_upper[1]), col="green")
	for (i in 5:(length(IL_par)+1)) {
		hist(P_VAL_PARAMS["IL",TR,CNR,AR1,i,], main=paste("IL par", i-1))
		#-1, because the very 1st entry a p-value
		abline(v=IL_par[i-1], col="red")
		#-2, because A2 and A3 are derived from rest, -1, because the very 1st entry a p-value
		abline(v=c(IL_par_lower[i-3], IL_par_upper[i-3]), col="green")
	}
dev.off()
pdf("par_dist_FIR.pdf")
	par(mfrow=c(4,3))
	hist(P_VAL_PARAMS["FIR",TR,CNR,AR1,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	for (i in 1:length(FIR_par)) {
		hist(P_VAL_PARAMS["FIR",TR,CNR,AR1,i+1,], main=paste("FIR tp", i))
		abline(v=FIR_par[i], col="red")
	}
dev.off()

#ALL COMBINATIONS: MORE VARIABLE ESTIMATES
pdf("par_dist_gamma2_NL_all.pdf")
	par(mfrow=c(4,2))
	hist(P_VAL_PARAMS["gamma2_nonlinear",,,,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	for (i in 1:length(gamma2_par)) {
		hist(P_VAL_PARAMS["gamma2_nonlinear",,,,i+1,], main=paste("gamma2 par", i))
		abline(v=gamma2_par[i], col="red")
		abline(v=c(gamma2_par_lower[i], gamma2_par_upper[i]), col="green")
	}
dev.off()
pdf("par_dist_IL_all.pdf")
	par(mfrow=c(4,2))
	hist(P_VAL_PARAMS["IL",,,,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	hist(P_VAL_PARAMS["IL",,,,2,], main=paste("IL par 1"))
	abline(v=IL_par[1], col="red")
	abline(v=c(IL_par_lower[1], IL_par_upper[1]), col="green")
	for (i in 5:(length(IL_par)+1)) {
		hist(P_VAL_PARAMS["IL",,,,i,], main=paste("IL par", i-1))
		#-1, because the very 1st entry a p-value
		abline(v=IL_par[i-1], col="red")
		#-2, because A2 and A3 are derived from rest, -1, because the very 1st entry a p-value
		abline(v=c(IL_par_lower[i-3], IL_par_upper[i-3]), col="green")
	}
dev.off()
pdf("par_dist_FIR_all.pdf")
	par(mfrow=c(4,3))
	hist(P_VAL_PARAMS["FIR",,,,1,], main="p values")
	abline(v=0.05, col="red"); abline(v=c(0,1), col="green")
	for (i in 1:length(FIR_par)) {
		hist(P_VAL_PARAMS["FIR",,,,i+1,], main=paste("FIR tp", i))
		abline(v=FIR_par[i], col="red")
	}
dev.off()

