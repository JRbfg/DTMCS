### set parameters
project <- "Emmerich_TH"
simulation_title <- "run_1"
n_iterations <- 100
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
library(moments)
library(sf)
library(terra)

# change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create result folders
dir.create(paste0("results"))
dir.create(paste0("results/",project))
dir.create(paste0("results/",project,"/",simulation_title))
  
# list all BEP-files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                       file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)

# read input parameter settings
inputs_df <- read.table("inputs.csv", header = T, sep=";")

# create random values based on input parameter settings
if(unique(inputs_df$W1_min!=inputs_df$W1_max)) random_vec_w1 <- sample(seq(unique(inputs_df$W1_min),unique(inputs_df$W1_max),(unique(inputs_df$W1_max)-unique(inputs_df$W1_min))/(n_iterations-1)), n_iterations)
if(unique(inputs_df$zc_min!=inputs_df$zc_max)) random_vec_zc <- sample(seq(unique(inputs_df$zc_min),unique(inputs_df$zc_max),(unique(inputs_df$zc_max)-unique(inputs_df$zc_min))/(n_iterations-1)), n_iterations)
if(unique(inputs_df$W1_min==inputs_df$W1_max)) random_vec_w1 <- rep(unique(inputs_df$W1_min), n_iterations)
if(unique(inputs_df$zc_min==inputs_df$zc_max)) random_vec_zc <- rep(unique(inputs_df$zc_min), n_iterations)
if(!is.na(inputs_df$W2_min)){
  if(unique(inputs_df$W2_min==inputs_df$W2_max)) random_vec_w2 <- rep(unique(inputs_df$W2_min), n_iterations)
  if(unique(inputs_df$W2_min!=inputs_df$W2_max)) random_vec_w2 <- sample(seq(unique(inputs_df$W2_min),unique(inputs_df$W2_max),(unique(inputs_df$W2_max)-unique(inputs_df$W2_min))/(n_iterations-1)), n_iterations)
}else{
  random_vec_w2 <- rep(NA, n_iterations)
}


inputs <- data.frame(W1=random_vec_w1, W2=random_vec_w2, zc=random_vec_zc)

params <- data.frame(W1_min=floor(min(inputs$W1)), W1_max=ceiling(max(inputs$W1)),
                     W2_min=floor(min(inputs$W2, na.rm=T)), W2_max=ceiling(max(inputs$W2, na.rm=T)),
                     zc_min=min(inputs$zc), zc_max=max(inputs$zc))
if(is.na(unique(inputs$W2))) params[1,3:4] <- NaN

# export input parameter settings  
write.table(params,paste0("results/",project,"/",simulation_title,"/params.txt"),row.names = F)
  
# execute zerocrossing over all availlable BEPs
for(measurement in unique(BEP_list$folder)){
  # define path to profiles 
  print(measurement)
  BEP_list_d <- BEP_list[BEP_list$folder==measurement, ]
  
  for(p in 1:nrow(BEP_list_d)){

    # read BEP
    BEP <- read.table(BEP_list_d[p,2], header = F)
    ID <- unique(BEP[ ,4])
    print(paste0("BEP: ", ID))
      
    # create result folders
    dir.create(paste0("results/", project,"/", simulation_title,"/", measurement))
    path <- paste0("results/", project,"/", simulation_title,"/", measurement,"/",ID)
    dir.create(path)
    dir.create(paste0(path, "/baselines"))
    dir.create(paste0(path, "/statistics"))
    
    # calculate additional parameters
    n_layers <- length(inputs[1,1:2][!is.na(inputs[1,1:2])])
    sample_frequency <- round(median(ave(BEP[ ,1], FUN=function(x) c(0,diff(x)))),2)

    inputs_p <- inputs
    inputs_p$W1 <- round(inputs$W1/sample_frequency)
    inputs_p$W2 <- round(inputs$W2/sample_frequency)

    # for-loop over all MCS iterations
    for (i in 1:nrow(inputs_p)) {

      # define break-off criteria
      if (n_layers > 1) {
        if (ceiling(inputs_p[i,1]) >= floor(inputs_p[i,2])) 
          next
        if (ceiling(inputs_p[i,2]) >= nrow(BEP)) 
          next
      }

      print(paste0("zerocrossing: iteration=", i))
      # run zerocrossing procedure
      source(paste0("zc_functions.R"))
      zerocrossing_scale_1(input=inputs_p[i, ], iteration=i, BEP, 
                                   n_layers, path, min_modus=F, sample_frequency)
    }
    
    # exoprt simulation settins
    settings <- data.frame(BEP=ID, iteration=rownames(inputs_p), inputs_p[ ,1:2]*sample_frequency, inputs_p[ ,3])
    write.table(settings, paste0(path, "/settings.txt"), row.names = F)
  }
}



