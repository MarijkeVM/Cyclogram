
# Packages
#install.packages(IntClust)
#install.packages(circlize)
library(IntClust)
library(circlize)

# Data
load("Script/Data/Fingerprints.RData")
head(Fingerprints)[,1:5]
#                   FP 1  FP 2  FP 3  FP 4  FP 5
# metformin        FALSE FALSE FALSE FALSE FALSE
# phenformin       FALSE FALSE FALSE FALSE FALSE
# phenyl biguanide FALSE FALSE FALSE FALSE FALSE
# estradiol        FALSE FALSE FALSE FALSE FALSE
# dexamethasone    FALSE FALSE FALSE FALSE FALSE
# verapamil        FALSE FALSE FALSE FALSE FALSE

load("Script/Data/TargetPredictions.RData")
head(TargetPredictions)[,1:5]
#                  TP 1 TP 2 TP 3 TP 4 TP 5
# metformin           0    0    0    0    0
# phenformin          0    0    0    0    0
# phenyl biguanide    0    0    0    0    0
# estradiol           0    0    1    0    0
# dexamethasone       0    0    1    0    0
# verapamil           0    0    0    0    0


# Cyclogram
CircularComparePlot<-function (List, nrclusters = NULL, cols = NULL, fusionsLog = FALSE, 
		WeightClust = FALSE, names = NULL, canvaslims=c(-1.0,1.0,-1.0,1.0), margins = c(8.1, 3.1, 
				3.1, 4.1), Highlight=NULL, plottype = "new", location = NULL) 
{
	plottypein <- function(plottype, location) {
		if (plottype == "pdf" & !(is.null(location))) {
			grDevices::pdf(paste(location, ".pdf", sep = ""))
		}
		if (plottype == "new") {
			grDevices::dev.new()
		}
		if (plottype == "sweave") {
		}
	}
	plottypeout <- function(plottype) {
		if (plottype == "pdf") {
			grDevices::dev.off()
		}
	}
	for (i in 1:length(List)) {
		if (attributes(List[[i]])$method == "Weighted" & WeightClust == 
				TRUE) {
			T = List[[i]]$Clust
			attr(T, "method") = "Single Clustering"
			List[[i]] = T
		}
	}
	MatrixColors = ReorderToReference(List, nrclusters, fusionsLog, 
			WeightClust, names)
	Names = ColorsNames(MatrixColors, cols)
	colnames(Names)=colnames(MatrixColors)
	nobs = dim(MatrixColors)[2]
	nmethods = dim(MatrixColors)[1]
	if (is.null(names)) {
		for (j in 1:nmethods) {
			names[j] = paste("Method", j, sep = " ")
		}
	}
	plottypein(plottype, location)
	
	circlize::circos.initialize(factors =c(1:ncol(MatrixColors)) , xlim = c(0, ncol(MatrixColors)))
	for(i in 1:length(List)){
		circlize::circos.trackPlotRegion(factors = c(1:ncol(MatrixColors)), ylim = c(0,1),track.height = 0.05)
	}
	track=c(length(List):1)
	for(i in 1:nrow(MatrixColors)){
		for(j in 1:ncol(MatrixColors)){
			circlize::highlight.sector(sector.index=j, track.index = track[i],col=Names[i,j])
		}
	}
	hc=List[[1]]$Clust
	#labels=substr(hc$order.lab,1,5)
	print(hc$order.lab)
	labels=hc$order.lab
	if(!is.null(Highlight)){
		for(h in 1:length(Highlight)){
			Name=names(Highlight)[h]
			#HL=which(labels%in%substr(Highlight[[h]],1,5))
			HL=which(labels%in%Highlight[[h]])
			
			Sims=c()
			for(i in 1:length(List)){
				Values=List[[i]]$DistM[lower.tri(List[[i]]$DistM)]
				Sims=c(Sims,as.numeric(1-List[[i]]$DistM[Highlight[[h]],Highlight[[h]]][lower.tri(List[[i]]$DistM[Highlight[[h]],Highlight[[h]]])]))
			}
			MedSim=round(stats::median(Sims),2)
			
			circlize::draw.sector(circlize::get.cell.meta.data("cell.start.degree", sector.index = min(HL)),
					circlize::get.cell.meta.data("cell.end.degree", sector.index = max(HL)),
					rou1 = 1, col = "#00000020")
			
			circlize::highlight.sector(sector.index=c(min(HL):max(HL)), track.index = 1, text = paste(Name,": ",MedSim,sep=""),
					facing = "bending.inside", niceFacing = TRUE, text.vjust = -1.5)
			
		}
	}
	circlize::circos.clear()
	graphics::par(new = TRUE)
	hc=List[[1]]$Clust
	max_height=max(hc$height)
	dend=stats::as.dendrogram(hc)
	#labels=substr(hc$order.lab,1,5)
	labels=hc$order.lab
	#ct=cutree(dend,6)
	circlize::circos.par("canvas.xlim" = c(canvaslims[1], canvaslims[2]), "canvas.ylim" = c(canvaslims[3], canvaslims[4]))
	circlize::circos.initialize(factors =1 , xlim = c(0, ncol(MatrixColors)))
	circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.4,bg.border=NA,
			panel.fun = function(x, y) {
				for(i in seq_len(ncol(MatrixColors))) {
					circlize::circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5),
							facing = "clockwise", niceFacing = TRUE,
							col = Names[1,colnames(MatrixColors)[i]], 
							cex = 0.6,font=2)
				}
			})
	circlize::circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,track.height = 0.4, panel.fun = function(x, y) {
				circlize::circos.dendrogram(dend, max_height = max_height)})
	circlize::circos.clear()
	
	
	
	
	plottypeout(plottype)
}


# Example
FP_Clust=Cluster(Data=Fingerprints,type="data",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=20,StopRange=TRUE) 
TP_Clust=Cluster(Data=TargetPredictions,type="data",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=20,StopRange=TRUE) 
Weighted_Clust=


