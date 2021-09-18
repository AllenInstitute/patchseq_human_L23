#############################################################################
## Setup and library loading
setwd("Z:\\hct\\HCT_RNAseq\\Jeremy\\patchseq_analysis\\human_IVSCC_paper_July2020\\human_cell_sizes\\")
library(RImageJROI)
library(msir)
library(ggplot2)

# Find the cases
cases = dir()
cases = cases[grepl("case",cases)]

# Pixel to um conversion (100um ~ 95 pixels)
convert = 100/95

#############################################################################
## Read in the ROIs and layers per case and plot them as a sanity check

rois <- layer <- list()
cols = setNames(c("blue","green","purple","red","cyan"),cases)

pdf("image_rois.pdf",height=10,width=7)
for (case in cases){
  # Read in the layer boundaries
  layer[[case]] <- read.ijzip(paste0(case,"/lines.zip"))
  plot(layer[[case]],asp=TRUE,main=case,col=cols[case])
  
  # Read in the rois
  fn <- dir(case)[substr(dir(case),1,1)=="L"]
  rois[[case]] <- list()
  for (f in fn){
    rois[[case]][[f]] <- read.ijzip(paste0(case,"/",f), verbose = FALSE)
    plot(rois[[case]][[f]],TRUE,col=cols[case])
  }
}
dev.off()


#############################################################################
## Calculate matrix of L1-4 coordinates for each x position

layerMat <- list()
for (case in cases){
  # Determine the range and values to assess layer boundaries
  xVals <- yVals <- xRange <- NULL
  len <- length(layer[[case]])
  for (i in 1:len) xVals <- rbind(xVals,c(layer[[case]][[i]]$x1,layer[[case]][[i]]$x2))
  xRange <- range(xVals)
  xRange <- round(c(xRange[1]-0.05*xRange[2],xRange[2]+0.05*xRange[2]))
  for (i in 1:len) yVals <- rbind(yVals,c(layer[[case]][[i]]$y1,layer[[case]][[i]]$y2))
  
  # Assess layer boundaries
  layerTmp <- matrix(0,nrow=diff(xRange)+1,ncol=5)
  colnames(layerTmp) <- c("x","L1","L2","L3","L4")
  layerTmp[,"x"] <- xRange[1]:xRange[2]
  m  <- (yVals[,2]-yVals[,1])/(xVals[,2]-xVals[,1])
  b  <- xVals[,1]
  y1 <- yVals[,1]
  layerTmp[,2:5] <- t(apply(layerTmp,1,function(x,j) round((x[1]-b[j])*m[j]+y1[j]), 1:4))
  if(length(b)>4) for (k in 2:round(length(b)/4)){
    xMin <- mean(b[(1:4)+4*(k-1)])
	kp   <- layerTmp[,1]>xMin 
	layerTmp[kp,2:5] <- t(apply(layerTmp[kp,],1,function(x,j) round((x[1]-b[j])*m[j]+y1[j]), (1:4)+4*(k-1)))
  }
  layerTmp <- as.data.frame(layerTmp)
  
  # Calculate layer depths, etc. at each x position for scaling
  layerTmp$L23_depth <- layerTmp$L4-layerTmp$L2
  layerMat[[case]]  <- layerTmp
  
  # Print average depth for comparison to mouse
  print(paste("Average L23 depth for",case,"=",round(mean(layerMat[[case]]$L23_depth)*convert)))
}
#[1] "Average L23 depth for case_072 = 1058"
#[1] "Average L23 depth for case_246 = 1293"
#[1] "Average L23 depth for case_316 = 1387"
#[1] "Average L23 depth for case_360 = 1073"
#[1] "Average L23 depth for case_732 = 1344"


#############################################################################
## For each ROI/case, determine the cell area as well as the x and y medians 

area_from_roi <- function(coord){
  area = 0
  for (x in unique(coord[,1])){
    area = area+diff(range(coord[coord[,1]==x,2]))
  }
  area
}


roiInfo <- list()
for (case in cases){
  fn  <- dir(case)[substr(dir(case),1,1)=="L"]
  roiInfo[[case]] <- NULL
  for (f in fn){
    roi    <- rois[[case]][[f]]
    roiTmp <- matrix(0,nrow=length(roi),ncol=4)
    colnames(roiTmp) <- c("x","depth","area","scaled_depth")
	for (i in 1:length(roi)){
	  coord <- roi[[i]]$coords
	  roiTmp[i,1:2] <- colMeans(coord)
	  roiTmp[i,1]   <- round(roiTmp[i,1])
	  roiTmp[i,3]   <- area_from_roi(coord)*convert*convert
	  lay <- layerMat[[case]][layerMat[[case]][,1]==roiTmp[i,1],]
	  roiTmp[i,4]   <- (lay$L2-roiTmp[i,2])/lay$L23_depth
    }
	roiInfo[[case]] <- rbind(roiInfo[[case]], roiTmp)
  }
  roiInfo[[case]][,4] <- pmin(pmax(roiInfo[[case]][,4],-1),0)
}
# We now have everything we need to do the depth comparison analysis!


#############################################################################
## Calculate and plot the average and standard deviation of cell size at each scaled depth for each donor

fits <- list()
l2   <- NULL
for (case in cases){
  fits[[case]] <- loess.sd(as.data.frame(roiInfo[[case]][,4:3]))
  l2 <- c(l2,-mean((layerMat[[case]]$L3-layerMat[[case]]$L2))/mean(layerMat[[case]]$L23_depth))
}

pdf("cell_size_vs_depth_human.pdf",height=7,width=7)
par(mar=c(6,6,6,6))
plot(roiInfo[[1]][,4:3], main = "Cortical depth vs. cell size (human)", ylim=c(0,800), xlim=c(-1,0),
  xlab="Scaled cortical depth (-1=L3/L4, ~0.18=L2/L3, 0=L1/L2)",ylab="Area of cell (um^2)",col="white")
for(case in cases){ 
  lines(fits[[case]], col = cols[case], lwd=2)
  lines(fits[[case]]$x, fits[[case]]$upper, lty=2, col = cols[case])
  lines(fits[[case]]$x, fits[[case]]$lower, lty=2, col = cols[case])
}
abline(v=0)

for(case in cases){ 
  plot(roiInfo[[case]][,4:3], main = case, ylim=c(0,800), xlim=c(-1,00),
    xlab="Scaled cortical depth (>-1=L4, [-1,0]=L3, >0=L2)",ylab="Area of cell (# pixels in ROI)",col="white")
  points(roiInfo[[case]][,4:3],pch=19,cex=0.2, col = cols[case])
  lines(fits[[case]], col = cols[case], lwd=2)
  lines(fits[[case]]$x, fits[[case]]$upper, lty=2, col = cols[case])
  lines(fits[[case]]$x, fits[[case]]$lower, lty=2, col = cols[case])
  abline(v=0)
}
dev.off()

#############################################################################
## Calculate cell density and heterogeneity vs. depth

# NOTE: I CAN ADD THE MOUSE DATA TO THIS PLOT DIRECTLY
# I could probably do the same thing for the cell count above, maybe...

breaks <- -(20:0)*(1/20)  #-(10:0)*(1/10)
b2     <- (breaks[2:length(breaks)]+breaks[1:(length(breaks)-1)])/2
dens   <- NULL
cellsPerMMh <- NULL
for(case in cases){
  rng  <- roiInfo[[case]][,"scaled_depth"]
  a <- h <- NULL
  for (b in 1:(length(breaks)-1)){
    tmp <- roiInfo[[case]][(rng>=breaks[b])&(rng<breaks[b+1]),"area"]
	a   <- c(a,mean(tmp))
	h   <- c(h,sd(tmp))
	q9  <- c(q9,as.numeric(quantile(tmp,0.95)))
  }
  dens <- rbind(dens,cbind(hist(rng, breaks=breaks, plot=FALSE)$density,h,a,q9,b2,case))
  
  cpm  <- diff(range(roiInfo[[case]][,"x"]))*diff(range(roiInfo[[case]][roiInfo[[case]][,"x"]<=200,"depth"]))
  cpm  <- 1000000*dim(roiInfo[[case]])[1]/(cpm*convert*convert)
  cellsPerMMh <- c(cellsPerMMh,cpm)
}
colnames(dens) <- c("dens","het","area","q90","group","case")
dens <- as.data.frame(dens)
dens$dens <- as.numeric(as.character(dens$dens))
dens$area <- as.numeric(as.character(dens$area))
dens$q90  <- as.numeric(as.character(dens$q90))
dens$het  <- as.numeric(as.character(dens$het))
dens$group <- factor(dens$group,levels=as.character(b2))
dens$species <- "human"

## SAVE THE HUMAN DATA FOR COMPARISON WITH MOUSE
densH <- dens
l2H   <- l2
roiInfoH <- roiInfo
save(densH, l2H, cellsPerMMh, breaks, roiInfoH, file="../mouse_cell_sizes/human_cell_data.RData")
save(densH, l2H, cellsPerMMh, breaks, roiInfoH, file="../mouse_cell_sizes_TEa/human_cell_data.RData")