#############################################################################
## Setup and library loading
setwd("Z:\\hct\\HCT_RNAseq\\Jeremy\\patchseq_analysis\\human_IVSCC_paper_July2020\\mouse_cell_sizes_TEa\\")
library(RImageJROI)
library(msir)
library(ggplot2)
library(dplyr)

# Find the cases
cases = dir()
cases = cases[grepl("case",cases)]

# Pixel to um conversion (100um ~ 165 pixels) -- THIS NEEDS TO BE UPDATED!!!
convert = setNames(c(100/330,100/330,100/330),cases)  # I measured ~300, but likely 330 for consistency with mouse VISp images

#############################################################################
## Read in the ROIs and layers per case and plot them as a sanity check

rois <- layer <- list()
cols = setNames(c("brown","orange","darkgrey"),cases)  # "blue","green","purple","red","cyan"

pdf("image_rois_TEa.pdf",height=10,width=7)
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
  colnames(layerTmp) <- c("x","L1","L2","L4","L6b")
  # L1 and L2 are swapped for the first case, so swap them
  if (case==cases[1])
    colnames(layerTmp) <- c("x","L2","L1","L4","L6b")
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
  print(paste("Average L23 depth for",case,"=",round(mean(layerMat[[case]]$L23_depth)*convert[case])))
}
#VISp[1] "Average L23 depth for case_set_1 = 281"
#VISp[1] "Average L23 depth for case_set_2 = 223"
#VISp[1] "Average L23 depth for case_set_3 = 257"
#TEa [1] "Average L23 depth for case_526976 = 167"
#TEa [1] "Average L23 depth for case_527216 = 156"
#TEa [1] "Average L23 depth for case_527859 = 174"


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
	  roiTmp[i,3]   <- area_from_roi(coord)*convert[case]*convert[case]
	  lay <- layerMat[[case]][layerMat[[case]][,1]==roiTmp[i,1],]
	  roiTmp[i,4]   <- (lay$L2-roiTmp[i,2])/lay$L23_depth
    }
	roiInfo[[case]] <- rbind(roiInfo[[case]], roiTmp)
	roiInfo[[case]][,4] <- pmin(pmax(roiInfo[[case]][,4],-1),0)
  }
}
# We now have everything we need to do the depth comparison analysis!

## To ensure we only include parts of the tissue where the complete L2/3 is 
##  included, only keep the cells between 500 and 1500 in the x direction
for (case in cases){
 kp = (roiInfo[[case]][,"x"]<1500)&(roiInfo[[case]][,"x"]>=500)
 roiInfo[[case]] = roiInfo[[case]][kp,]
 kp = (layerMat[[case]][,"x"]<1500)&(layerMat[[case]][,"x"]>=500)
 layerMat[[case]] <- layerMat[[case]][kp,c("x","L1","L2","L4","L6b","L23_depth")]
}

#############################################################################
## Calculate and plot the average and standard deviation of cell size at each scaled depth for each donor

fits <- list()
pdf("cell_size_vs_depth_mouse_TEa.pdf",height=7,width=7)
par(mar=c(6,6,6,6))
plot(roiInfo[[1]][,4:3], main = "Cortical depth vs. cell size (mouse)", ylim=c(0,800), xlim=c(-1,0),
  xlab="Scaled cortical depth (-1=L3/L4, ~0.25=L2/L3, 0=L1/L2)",ylab="Area of cell (# pixels in ROI)",col="white")
for(case in cases){ 
  lines(fits[[case]], col = cols[case], lwd=2)
  lines(fits[[case]]$x, fits[[case]]$upper, lty=2, col = cols[case])
  lines(fits[[case]]$x, fits[[case]]$lower, lty=2, col = cols[case])
}
abline(v=0)

for(case in cases){ 
  plot(roiInfo[[case]][,4:3], main = case, ylim=c(0,800), xlim=c(-1,0),
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

## LOAD THE HUMAN MTG AND MOUSE VISp DATA FOR COMPARISON WITH MOUSE TEa
load(file="human_cell_data.RData")
load(file="mouse_cell_data.RData")

# breaks carried over from human analysis
b2     <- (breaks[2:length(breaks)]+breaks[1:(length(breaks)-1)])/2
dens   <- NULL
cellsPerMM <- NULL

for(case in cases){
  rng  <- roiInfo[[case]][,"scaled_depth"]
  a <- h <- q9 <- NULL
  for (b in 1:(length(breaks)-1)){
    tmp <- roiInfo[[case]][(rng>=breaks[b])&(rng<breaks[b+1]),"area"]
	a   <- c(a,mean(tmp))
	h   <- c(h,sd(tmp))
	q9  <- c(q9,as.numeric(quantile(tmp,0.95)))
  }
  dens <- rbind(dens,cbind(hist(rng, breaks=breaks, plot=FALSE)$density,h,a,q9,b2,case))
  
  cpm  <- diff(range(roiInfo[[case]][,"x"]))*diff(range(roiInfo[[case]][roiInfo[[case]][,"x"]<=600,"depth"]))  # Was 200 in other scripts
  cpm  <- 1000000*dim(roiInfo[[case]])[1]/(cpm*convert[case]*convert[case])
  cellsPerMM <- c(cellsPerMM,cpm)
}
colnames(dens) <- c("dens","het","area","q90","group","case")
dens <- as.data.frame(dens)
dens$dens <- as.numeric(as.character(dens$dens))
dens$area <- as.numeric(as.character(dens$area))
dens$q90  <- as.numeric(as.character(dens$q90))
dens$het  <- as.numeric(as.character(dens$het))
dens$group <- factor(dens$group,levels=as.character(b2))
dens$species <- "mouse_TEa"
densT <- dens

# Merge human and mouse
dens <- rbind(densT,densM,densH)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}


# OUTPUT SPECIES DIFFERENCES IN CROSS-SECTIONAL AREA
dpt   <- as.numeric(as.character(dens$group))
l2dex <- data_summary(dens[dpt >= -0.15,],"area",c("species"))
l3dex <- data_summary(dens[(dpt >= -0.95)&(dpt < -0.7),],"area",c("species"))
print(paste("Average human / mouse cross-sectional area (L2):",signif(l2dex[1,"area"]/l2dex[3,"area"],2)))
print(paste("Average human / mouse cross-sectional area (deep L3):",signif(l3dex[1,"area"]/l3dex[3,"area"],2)))
#VISp[1] "Average human / mouse cross-sectional area (L2): 2.1"
#VISp[1] "Average human / mouse cross-sectional area (deep L3): 3.0"
#TEa [1] "Average human / mouse cross-sectional area (L2): 2.9"
#TEa [1] "Average human / mouse cross-sectional area (deep L3): 5"

l2dex <- data_summary(dens[dpt >= -0.15,],"het",c("species"))
l3dex <- data_summary(dens[(dpt >= -0.95)&(dpt < -0.7),],"het",c("species"))
print(paste("SD human / mouse cross-sectional area (L2):",signif(l2dex[1,"het"]/l2dex[3,"het"],2)))
print(paste("SD human / mouse cross-sectional area (deep L3):",signif(l3dex[1,"het"]/l3dex[3,"het"],2)))
#VISp[1] "SD human / mouse cross-sectional area (L2): 2.4"
#VISp[1] "SD human / mouse cross-sectional area (deep L3): 4.0"
#TEa [1] "SD human / mouse cross-sectional area (L2): 3.9"
#TEa [1] "SD human / mouse cross-sectional area (deep L3): 6.2"

# PLOT CELL DENSITY
df2 <- data_summary(dens,"dens",c("group","species"))
df2 <- df2[!is.element(df2$group,as.character(c(b2[1],b2[length(b2)]))),]  # Remove bookends

p <- ggplot(df2, aes(x=group, y=dens, group=species, color=species)) + 
  geom_line(size=2.5) +
  geom_point()+
  scale_color_manual(values=c("#4D8000","#A84920","grey")) +  # ADD grey for mouse TEa
  geom_errorbar(aes(ymin=dens-sd, ymax=dens+sd), width=0.25,
                 position=position_dodge(0.05),size=2)

p + labs(title="Density", x="Scaled depth", y = "Density") + 
    ylim(0,3) + 
	theme_classic() + 
	geom_hline(yintercept=0)

ggsave("cell_density_vs_depth_TEa.pdf",height=7,width=7)


# PLOT CELL AREA
df3 <- data_summary(dens,"area",c("group","species"))
df3 <- df3[!is.element(df3$group,as.character(c(b2[1],b2[length(b2)]))),]  # Remove bookends

p <- ggplot(df3, aes(x=group, y=area, group=species, color=species)) + 
  geom_line(size=2.5) +
  geom_point()+
  scale_color_manual(values=c("#4D8000","#A84920","grey")) +
  geom_errorbar(aes(ymin=area-sd, ymax=area+sd), width=0.25,
                 position=position_dodge(0.05),size=2)

p + labs(title="Cell size", x="Scaled depth", y = "Cell area") + 
    ylim(0,500) + 
	theme_classic() + 
	geom_hline(yintercept=0)

ggsave("cell_area_vs_depth_TEa.pdf",height=7,width=7)



# PLOT CELL HETERGENIETY (e.g., sd of cell area)
df3 <- data_summary(dens,"het",c("group","species"))
df3 <- df3[!is.element(df3$group,as.character(c(b2[1],b2[length(b2)]))),]  # Remove bookends

p <- ggplot(df3, aes(x=group, y=het, group=species, color=species)) + 
  geom_line(size=2.5) +
  geom_point()+
  scale_color_manual(values=c("#4D8000","#A84920","grey")) +
  geom_errorbar(aes(ymin=het-sd, ymax=het+sd), width=0.25,
                 position=position_dodge(0.05),size=2)

p + labs(title="Cell heterogeneity", x="Scaled depth", y = "Cell heterogeneity") + 
    ylim(0,250) + 
	theme_classic() + 
	geom_hline(yintercept=0)
ggsave("cell_heterogeneity_vs_depth_TEa.pdf",height=7,width=7)


## Print the absolute density
dh = 40*cellsPerMMh/1000 
dm = 40*cellsPerMMm/1000
dt = 40*cellsPerMM/1000
print(paste("Mouse density (mm^3) =",signif(mean(dt),3),"+/-",signif(sd(dt),3)))
#VISp[1] "Mouse density (mm^3) = 165 +/- 24.9"
#TEa [1] "Mouse density (mm^3) = 261 +/- 42.2"
print(paste("Human density (mm^3) =",signif(mean(dh),3),"+/-",signif(sd(dh),3)))
#[1] "Human density (mm^3) = 27.7 +/- 4.46"

## Plot the absolute density
cpm <- data.frame(species=c(rep("Human",5),rep("Mouse",3),rep("Mouse_TEa",3)),density=c(dh,dm,dt))

cpm.summary <- cpm %>%
  group_by(species) %>%
  dplyr::summarise(
    sd = sd(density, na.rm = TRUE),
    density = mean(density)
	
  )

p <- ggplot(cpm, aes(y=density, x=species)) + 
  geom_bar(data=cpm.summary, aes(fill = species), color="black", stat="identity") +
  xlab("") + ylab("Cell density (thousands of cells mm^3)") +
  scale_fill_manual(values=c("#4D8000","#A84920","grey")) +
  geom_point(color = "black", shape="-", size=7) + 
  ylim(0,310) 
p

ggsave("absolute_cell_density_TEa.pdf",height=7,width=4)


## Plot a histogram of human and mouse cell areas.  

pdf("cell_area_histograms_TEa.pdf",height=10,width=5)
par(mfrow=c(8,1))
par(mar=c(2,3,0,0))
for (i in 1:3){
  hist(roiInfo[[i]][,"area"],breaks=(0:50)*20,xlab="",main="",ylab="Frequency",col="red")
  abline(v=mean(roiInfo[[i]][,"area"])) 
}
for (i in 1:5){
  hist(roiInfoH[[i]][,"area"],breaks=(0:50)*20,xlab="",main="",ylab="Frequency",col="green")
  abline(v=mean(roiInfoH[[i]][,"area"])) 
  print(quantile(roiInfoH[[i]][,"area"]))  
}
dev.off()

p=c(0,0.05,0.25,0.5,0.75,0.9,0.99,1)
a<-NULL
for (i in 1:3)  a <- c(a,roiInfo[[i]][,"area"])
print(quantile(a,p))
#        0%         5%        25%        50%        75%        90%        99%       100% 
#  4.499541  18.999082  44.719927  66.115702  96.464646 130.394858 178.681359 244.628099
head(-sort(-a))
#[1] 244.6281 216.7126 214.1414 210.4683 195.4086 194.3067

pdf("cell_diameter_histograms_TEa.pdf",height=10,width=5)
par(mfrow=c(8,1))
par(mar=c(2,3,0,0))
for (i in 1:3){
  d = 2*sqrt(roiInfo[[i]][,"area"]/pi)
  hist(d,breaks=0:40,xlab="",main="",ylab="Frequency",col="red")
  abline(v=mean(d)) 
}
for (i in 1:5){
  d = 2*sqrt(roiInfoH[[i]][,"area"]/pi)
  hist(d,breaks=0:40,xlab="",main="",ylab="Frequency",col="green")
  abline(v=mean(d)) 
}
dev.off()

a<-NULL
for (i in 1:5)  a <- c(a,roiInfoH[[i]][,"area"])
print(quantile(a,p))
#       0%        5%       25%       50%       75%       90%       99%      100% 
# 13.29640  62.04986 107.47922 168.42105 252.63158 343.49030 542.30471 886.42659
head(-sort(-a))
#[1] 886.4266 885.3186 824.3767 787.8116 745.7064 722.4377

