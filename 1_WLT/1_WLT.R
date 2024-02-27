### set general project parameters
project <- "Emmerich_TH"
###

# load required libraries
library(dplyr)
library(sp)
library(raster)
library(rgeos)
library(mapview)
library(raster)
library(forecast)
library(TTR)
library(zoo)
library(data.table)
library(rgdal)
library(stringr)
library(ggplot2)
library(geometry)
library(retistruct)
library(matrixStats)
library(matlabr)
library(purrr)
library(PEIP)
library(pracma)
library(stats)
library(warbleR)
library(dtt)
library(profvis)
library(gridExtra)
library(grid)
library(gtools)
library(Orcs)

# change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# set WLT parameters (optional)
mother_wavelet <- "MORLET"  # c('MORLET','MEXICANHAT')
wlt_parameter <- 6  # c(6,7,8,9,10) or c(2,6)
delta_freq <- 0.05  # c(0.01:0.5)

# list all files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                           file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)
measurements <- unique(BEP_list$folder)

# creating dataframes to store results
wlt_windows <- data.frame()
wlt_peaks <- data.frame()

# execute WLT over all files with different significance levels
for(signif_level in seq(0.3,0.9,0.3)){
  print(paste0("signif_level: ", signif_level))
  for(m in measurements){
    print(paste0("measurement: ",m))
    BEP_list_d <- BEP_list[BEP_list$folder==m, ]
    for(p in 1:nrow(BEP_list_d)){
      BEP <- read.table(BEP_list_d[p,2], header = F)
      ID <- unique(BEP[ ,4])
      source("wlt_functions.R")
      WLT_list <- WLT(BEP[ ,1:2], mother_wavelet, wlt_parameter, delta_freq, 
                          signif_level, ID)
      if(length(WLT_list)==1){
        peaks_all <- data.frame(WLT_list[1], BEP=ID, measurement=m, signif_level)
        wlt_peaks <- rbind(wlt_peaks, peaks_all)
      } else {
        peaks_sig <- data.frame(WLT_list[1], BEP=ID, measurement=m, signif_level)
        peaks_all <- data.frame(WLT_list[2], BEP=ID, measurement=m, signif_level)
        wlt_peaks <- rbind(wlt_peaks, peaks_all)
        wlt_windows <- rbind(wlt_windows, peaks_sig)
      }
    }
  }
}

# create result folders
dir.create(paste0("results"))
dir.create(paste0("results/wlt"))
dir.create(paste0("results/wlt/", project))
write.table(wlt_windows, paste0("results/wlt/", project,"/wlt_windows.txt"), row.names = F)

# Calculate ranges of respected wavelengths 
sd1=round(sd(wlt_windows$x[wlt_windows$BEP%in%c(8:14)&wlt_windows$layer==1]))
sd2=round(sd(wlt_windows$x[wlt_windows$BEP%in%c(8:14)&wlt_windows$layer==2]))
m1=round(mean(wlt_windows$x[wlt_windows$BEP%in%c(8:14)&wlt_windows$layer==1]))
m2=round(mean(wlt_windows$x[wlt_windows$BEP%in%c(8:14)&wlt_windows$layer==2]))

# Visualize results
plot(wlt_windows[ ,c(4,1)], xlim=c(1,16), ylim=c(0,30), pch=20, xlab="BEP", ylab="wavelength [m]")
grid()
polygon(x=c(8,14,14,8), y=c(m1+2*sd1, m1+2*sd1, m1-2*sd1, m1-2*sd1), col=rgb(1,0,0,0.5), border=NA)
polygon(x=c(8,14,14,8), y=c(m2+2*sd2, m2+2*sd2, m2-2*sd2, m2-2*sd2), col=rgb(0,0,1,0.5), border=NA)
points(wlt_windows[wlt_windows$layer==1,c(4,1)], col="red", pch=19)
points(wlt_windows[wlt_windows$layer==2,c(4,1)], col="blue", pch=19)

