### set parameters
project <- "Emmerich_TH"
s_parameter <- 24
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
    
# create result folders
dir.create("results")
dir.create(paste0("results/detrending/"))
dir.create(paste0("results/detrending/",project))

# list all files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                       file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)
measurements <- unique(BEP_list$folder)


for(measurement in measurements){
  print(measurement)
  BEP_list_m <- BEP_list[BEP_list$folder==measurement, ]
  dir.create(paste0("results/detrending/",project,"/",measurement))
  for(p in 1:nrow(BEP_list_m)){
    
    BEP <- read.table(BEP_list_m[p,2], header = F)
    ID <- unique(BEP[ ,4])
    
    source("detrending_functions.R")
    signal_data <- BEP[ ,2]  
    
    # defining the limits of the s-parameter
    s_param_power <- c(-12:-1, seq(from = 0, to = 60, by = 2))
    
    # creating df for smooth versions of the signal
    smooth_data <- zeros(length(signal_data), length(s_param_power))
    
    m <- mean(signal_data)
    sd <- std(signal_data)
    nsignal_data <- (signal_data - m)/sd
    
    # Actual value 12.5# of max amplitude and using normalized data
    param <- 0.125 * (max(nsignal_data) - min(nsignal_data))
    mintab <- peakdet(nsignal_data, param)
    mintab <- data.frame(mintab[2])
    troughSignal <- signal_data[mintab[, 1]]
    
    # calculating the smoothed signals for selected s-parameter
    Z <- SMOOTHN(signal_data, 2^(s_parameter))  # Base of the S Parameter
    delta <- 0
    troughSmooth <- data.frame()
    m <- mean(Z)
    sd <- std(Z)
    nsignal_data <- (Z - m)/sd
    param <- 0.125 * (max(nsignal_data) - min(nsignal_data))
    mintab <- peakdet(nsignal_data, param)
    mintab <- data.frame(mintab[2])
    if (nrow(mintab) > 0) {
      troughSmooth <- Z[mintab[, 1]]
    }
    if (length(troughSignal) == length(troughSmooth)) {
      delta <- mean(abs(troughSmooth - troughSignal))
    }
    Z <- Z - delta
    smooth_data <- Z
    
    if(all(is.na(smooth_data[1, ]))) next
    
    plot(BEP [ ,1:2], type="l", xlab="x [m]", ylab="z[m]", main=paste0("BEP ", ID))
    lines(data.frame(x=BEP[ ,1], z=smooth_data[1, ]), col="red", lwd=2)
    
    # detrending the BEP
    BEP[ ,2] <- BEP[,2] - smooth_data[1, ]
    #plot(BEP[ ,1:2], type="l", xlab="x [m]", ylab="z[m]")
    
    write.table(BEP,paste0("results/detrending/", project, "/", measurement,"/de_BEP_",ID,".txt"), row.names = F, col.names = F)
  }  
}



  
  
