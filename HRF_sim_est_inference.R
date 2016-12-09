

###############################################################################################
####   Simulation of fMRI time series using a 4-dim grid where different HRF models
####   (double-gamma, inverse logit, FIR), different values of TR, CNR and AR(1)
####   coefficient are considered. For each simulated fMRI time series the HRF model
####   is estimated and inference is conducted.
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     October-December 2016
###############################################################################################


#CRAN package with implemented simulated annealing
library(GenSA)
#numerical derivatives
library(numDeriv)
#parallelisation
library(parallel)
#3D bar plots
library(latticeExtra)
#par(mfrow) for lattice plots
library(gridExtra)

source("HRF_basic_functions.R")
source("HRF_sim.R")
source("HRF_est_inference.R")

#length of the time series
N=160
#stimulus times [seconds], up to time series length in sec when TR=2 multiplied with 10, so that different TRs work
#later in the simulaton only time points < max(t_i) ~ N considered
stimulus_T_sec=seq(16,304*10,by=32)
#length of kernel [seconds]
kernel_L_sec = 32
#in next weeks I will add randomized paradigms and proper convolution!!!

#canonical double gamma model
gamma2_par = c(1,6,16,1,1,1/6)
gamma2_par_lower = c(0, 1, 3, 0.1, 0.1, 0.01)
gamma2_par_upper = c(10, 12, 30, 2, 2, 2)
#parameters for IL chosen so that the evoked HRF as similar to canonical double gamma as possible
#IL_par has 9 values, lower and upper boundaries each has 7, because A2 and A3 are derived from rest
IL_par = c(0.24755081, 0.03119706, -0.27874786, 2.99383573, 16.85766860, 7.53035430, 0.72883708, 3.00000000, 1.91013312)
IL_par_lower = c(0, 0, 0, 0, 0.01, 0.01, 0.01)
IL_par_upper = c(10, 30, 30, 30, 5, 5, 5)
#parameters for FIR chosen arbitrarily
FIR_par = c(0,0,0,0,0,0.1,0.2,0.15,0.1,0.05,-0.05)

#maximum number of parameters, needed later to know how to save the parameter estimates
par_max = max(length(gamma2_par), length(IL_par), length(FIR_par))

HRF_sample = c("gamma2", "gamma2_D", "gamma2_DD", "gamma2_nonlinear", "IL", "FIR")
TR_sample = seq(1,4,by=0.5)
CNR_sample = seq(0.06,0.3,by=0.06)
AR1_sample = seq(0.1, 0.9, by=0.2)

#the simulation can take a long time, at the screen you will see the current iterations
#the number after mclappy -> number of iterations in the simulation
system.time({P_VAL_PARAMS = simplify2array(mclapply(1:14, function(i) {
	cat(i, "\n")
	P_VAL_PARAMS_i = array(-1, dim=c(length(HRF_sample), length(TR_sample), length(CNR_sample), length(AR1_sample), par_max+1))
	dimnames(P_VAL_PARAMS_i) = c(list(HRF_sample, TR_sample, CNR_sample, AR1_sample, NULL))
	#for HRF (now) 6 options considered
	for (HRF in HRF_sample) {
		#for all double gamma related HRF functions the same simulation used
		if (substr(HRF, 1, 6)=="gamma2") {
			HRF_short = "gamma2"
		} else {
			HRF_short = HRF
		}
		#for TR 7 options considered
		for (TR in TR_sample) {
			#t has to be defined each time again, since TR changes
			#assumed that the activation occurs exactly between two adjacent time points (no slice-timing correction!)
			t_i = seq(TR/2,N*TR-TR/2,by=TR)
			#only time points < max(t_i) ~ N considered
			stimulus_T = stimulus_T_sec[which(stimulus_T_sec<max(t_i))]
			#stimulus times, in [seconds/TR] = [scans]
			stimulus_T = round(stimulus_T/TR)
			#length of kernel, in [seconds/TR] = [scans]
			kernel_L = round(kernel_L_sec/TR)
			#for CNR 5 options considered
			for (CNR in CNR_sample) {
				#for AR1 5 options considered
				for (AR1 in AR1_sample) {
					fMRI = HRF_sim(HRF=HRF_short, t=t_i, TR=TR, CNR=CNR, AR1=AR1, stimulus_T=stimulus_T, kernel_L=kernel_L)
					#shift of 1 scan, uncomment when you want to consider shifts of the response!!!
					#fMRI[2:N] = fMRI[1:(N-1)]
					HRF_inference_R = HRF_est_inference(fMRI=fMRI, HRF=HRF, t=t_i, TR=TR, stimulus_T=stimulus_T, kernel_L=kernel_L, par_max=par_max, max.time=2)
					P_VAL_PARAMS_i[HRF, paste(TR), paste(CNR), paste(AR1), 1]             = HRF_inference_R$p_val
					P_VAL_PARAMS_i[HRF, paste(TR), paste(CNR), paste(AR1), 2:(par_max+1)] = HRF_inference_R$par
				}
			}
		}
	}
	return (P_VAL_PARAMS_i)
#the number of cores you want to use in the simulation
}, mc.cores=14))})

NUMBER_SIGN = array(-1, dim=c(length(HRF_sample), length(TR_sample), length(CNR_sample), length(AR1_sample)))
dimnames(NUMBER_SIGN) = c(list(HRF_sample, TR_sample, CNR_sample, AR1_sample))
for (HRF in HRF_sample) {
	for (TR in TR_sample) {
		for (CNR in CNR_sample) {
			for (AR1 in AR1_sample) {
				NUMBER_SIGN[HRF, paste(TR), paste(CNR), paste(AR1)] = sum(P_VAL_PARAMS[HRF, paste(TR), paste(CNR), paste(AR1), 1, ]<0.05)
			}
		}
	}
}

index = 0

#3D bar plots: values of 2 parameters varied: CNR, AR1
pdf("HRF_sim_results.pdf", width=16, height=10.5)
for (TR in TR_sample) {
	plots = list(NULL)
	for (HRF in HRF_sample) {
		index = index + 1
		plot = cloud(NUMBER_SIGN[HRF, paste(TR), , ],
			panel.3d.cloud=panel.3dbars, col.facet="grey", xlab="CNR", zlim=c(0,dim(P_VAL_PARAMS)[6]), ylab="AR1", zlab="#sign",
			xbase=0.4, ybase=0.4, main=paste("TR = ", TR, "s,       HRF = ", HRF, sep=""), scales=list(arrows=FALSE, col=1),
			par.settings = list(axis.line=list(col="transparent")), plot=F)
		plots[[(index-1)%%6+1]] = plot
	}
	#each page includes 6 figures, par(mfrow) doesn't work for cloud figures
	grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], ncol=3)
}
dev.off()

