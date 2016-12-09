

###############################################################################################
####   Different HRF models (double gamma and FIR) compared with each other in FSL using
####   FCP RS-fMRI scans and CRIC scans: RS and checkerboard.
####   Written by:  Wiktor Olszowy, University of Cambridge
####   Contact:     wo222@cam.ac.uk
####   Created:     September-October 2016
###############################################################################################


HRF_FSL_gamma2_FIR = function (data_source = "CRIC", model="FIR", paradigm="checkerboard", cores=16) {

	library(AnalyzeFMRI)
	library(oro.nifti)
	library(parallel)
	library(reshape2)

	path = paste("/home/wo222/HRF_modeling/scans_CRIC_", paradigm, "_", model, sep="")
	setwd(path)

	#function choosing the last 'n' characters from a string 'x'
	substrRight = function(x, n) {
		substr(x, nchar(x)-n+1, nchar(x))
	}

	#selecting only scans
	scans = sort(list.files()[which(substrRight(list.files(), 4)==".nii")])

	#for each subject only one scan taken, the most recent one
	#scans_single = scans[1]
	#for (i in (length(scans)-1):1) {
	#	if (substr(scans[i], 1, 5) != substr(scans[i+1], 1, 5)) {
	#		scans_single = c(scans_single, scans[i])
	#	} else {
	#		system(paste("rm", scans[i]))
	#	}
	#}
	#scans = scans_single

	#FSL 'pre-processing' + 'stats'
	res = mclapply(1:length(scans), function(i) {
		scan = scans[i]
		A = f.read.nifti.volume(scan)
		scan_short = substr(scan, 1, nchar(scan)-4)
		system(paste("cp /home/wo222/HRF_modeling/design_CRIC_", model, ".fsf ", path, "/design_CRIC_", scan_short, ".fsf" , sep=""))
		######### POSSIBLY HAS TO BE CHANGED FOR SOME CONSTELLATION OF THE OPTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		scan_def = "13676_20120222_0040_rotated.nii"
		system(paste("sed -i 's/'", scan_def, "'/'", scan, "'/g' ", "design_CRIC_", scan_short, ".fsf", sep=""))
		#default: FIR
		if (model=="gamma2") {
			option_def = paste("set fmri(convolve1) 6")
			option_new = paste("set fmri(convolve1) 3")
			system(paste("sed -i 's/", option_def, "/", option_new, "/g' design_CRIC_", scan_short, ".fsf", sep=""))
		}
		######### POSSIBLY HAS TO BE CHANGED FOR SOME CONSTELLATION OF THE OPTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		folder_def = paste("scans_CRIC_checkerboard_", model, sep="")
		folder_new = paste("scans_", data_source, "_", paradigm, "_", model, sep="")
		system(paste("sed -i 's/", folder_def, "/", folder_new, "/g' design_CRIC_", scan_short, ".fsf", sep=""))
		system(paste("feat design_CRIC_", scan_short, ".fsf", sep=""))
	}, mc.cores=cores)

	if (data_source == "CRIC") {
		FSL = array(0, dim=c(length(scans),64,64,32))
		masks = array(0, dim=c(length(scans),64,64,32))
	} else if (data_source == "FCP") {
		FSL = array(0, dim=c(length(scans),72,72,47))
		masks = array(0, dim=c(length(scans),72,72,47))
	}

	for (scan in scans) {
		#extracting zstat from FSL
		system(paste("gunzip ", path, "/", substr(scan, 1, nchar(scan)-4), ".feat/stats/zstat1.nii", sep=""))
		FSL[which(scans==scan),,,] = f.read.nifti.volume(paste(path, "/", substr(scan, 1, nchar(scan)-4), ".feat/stats/zstat1.nii", sep=""))[,,,1]
	}

	#the post-statistical analysis (cluster-based inference)
	post_stats = simplify2array(mclapply(1:length(scans), function(i) {
		scan = scans[i]
		setwd(paste(path, "/", substr(scan, 1, nchar(scan)-4), ".feat/stats/", sep=""))
		par = read.table("smoothness", header=F, sep=" ")
		system("rm cluster_zstat1.txt")
		system("rm cluster_index.nii")
		#3.09 is not default in FSL, but it is the default option in SPM and is suggested for FSL in Eklund et al. '16
		system(paste("cluster -i zstat1 -o cluster_index -t 3.09 -p 0.05 -d ", par$V2[1], " --volume=", par$V2[2], " > cluster_zstat1.txt", sep=""))
		system("gunzip cluster_index.nii.gz")
		mask = f.read.nifti.volume("cluster_index.nii")
		mask = mask > 0.5
		cluster = read.table("cluster_zstat1.txt", header=T, sep="\t")
		FSL_mask = mask[,,,1] * 1.0
		setwd(path)
		writeNIfTI(FSL_mask, filename=paste(substr(scan, 1, nchar(scan)-4), "_mask", sep=""))
		P_FSL_no = sum(cluster$Voxels)
		return(P_FSL_no)
	}, mc.cores=cores))

	setwd(path)
	for (i in 1:length(scans)) {
		scan = scans[i]
		system(paste("gunzip ", substr(scan, 1, nchar(scan)-4), "_mask.nii.gz", sep=""))
		masks[i,,,] = f.read.nifti.volume(paste(substr(scan, 1, nchar(scan)-4), "_mask.nii", sep=""))
	}
	masks_sum = apply(masks, c(2,3,4), "sum")

	return(list(scans=scans, post_stats=post_stats))

}

#barplots
library(ggplot2)
dd = data.frame(data_source=factor(c(rep("FCP: resting-state", 2), rep("CRIC: resting-state", 2), rep("CRIC: checkerboard", 2)),
	levels=c("FCP: resting-state", "CRIC: resting-state", "CRIC: checkerboard")),
	model=factor(rep(c("double gamma", "FIR"), 3)),
	#values manually copied from the output of the HRF_FSL_gamma2_FIR() function
	pos_prop=c(26.8, 21.2, 6.9, 2.7, 58.6, 35.7),
	vox_no_mean=c(79.3, 14.6, 87.0, 8.9, 366.2, 72.9))
pdf("HRF_models_means.pdf", width=5.6, height=5.6)
	ggplot(data=dd, aes(x=data_source, y=vox_no_mean, fill=model)) + geom_bar(stat="identity", position="dodge") +
		scale_fill_discrete(guide=guide_legend(title="HRF model \n")) +
		theme(legend.position=c(0.24, 0.84)) +
		scale_x_discrete(name="\n data source: paradigm") +
		scale_y_continuous(name="mean of the number of significant voxels \n")
dev.off()
pdf("HRF_models_props.pdf", width=5.6, height=5.6)
	ggplot(data=dd, aes(x=data_source, y=pos_prop, fill=model)) + geom_bar(stat="identity", position="dodge") +
		scale_x_discrete(name="\n data source: paradigm") +
		scale_y_continuous(name="proportion of significant scans [%]\n", limits=c(0,100)) +
		geom_hline(aes(yintercept=yint, colour="5% nominal FWER"), data=data.frame(yint=5, model="5% nominal FWER")) +
		theme(legend.position=c(0.24, 0.84)) +
		scale_colour_manual(name=" ", values=c("black")) +
		scale_fill_manual(name="HRF model             ", values=c("#F8766D", "#00BFC4"))
dev.off()

