### set parameters
project <- "Emmerich_detail"
simulation_title <- "run_1"
n_layers <- 2
n_iterations <- 100
# define maximnal lag [m] 
lag_max <- 50
# define sample frequency of BEPs
sample_freq <- 0.1
# define resolution for cc 
resolution <- 0.1
# define porosity and grain density for bedload estimation
porosity <- 0.3
density <- 2603
###

# load required libraries
library(gtools)
library(geometry)

# change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read table with information about time offsets between measurements
measurement_table <- read.table(paste0("../0_DATA/",project,"/measurements.csv"), header = T, sep=";")

# list all BEP-files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                       file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)

# create empty dataframes to store the results
migration_BEP <- data.frame()
migration_L1 <- data.frame()
migration_L2 <- data.frame()

lag_max <- round(lag_max/resolution)

# for-loop over all constellations of measurements
for(m in 1:(nrow(measurement_table))){

  print(paste("measurement pair: ",m))
    
  # define t1, t2 and delta_t
  t1 <- measurement_table[m,1]
  t2 <- measurement_table[m,2]
  delta_t <- measurement_table[m,3]
    
  # read BEP list
  BEP_list_t1 <- BEP_list[BEP_list$folder==t1, ]
  BEP_list_t2 <- BEP_list[BEP_list$folder==t2, ]
    
  BEPs_1 <- list.dirs(paste0("../2_ZC/results/",project,"/",simulation_title,"/",t1), recursive = F, full.names = F)
  BEPs_2 <- list.dirs(paste0("../2_ZC/results/",project,"/",simulation_title,"/",t2), recursive = F, full.names = F)
  if(length(BEPs_1)==0) print(paste0("WARNING: no zc-results available for measurement ", t1))
  
  BEPs <- BEPs_1[BEPs_1%in%BEPs_2]
  n_BEPs <- length(BEPs)
    
  for(p in 1:n_BEPs){
    
    # read BEPs for t1
    BEP_t1 <- read.table(BEP_list_t1[p,2])
    ID <- BEP_t1[1,4] 
    print(paste("BEP:",ID))
    n_sections <- length(unique(BEP_t1[ ,3]))
    # interpolate BEP 1 corresponding to defined resolution
    int <- approx(BEP_t1[,1:2],xout = seq(BEP_t1[1,1], BEP_t1[nrow(BEP_t1),1],resolution))      
    BEP_t1 <- data.frame(int$x,int$y)
      
    # read BEP for t2
    BEP_t2 <- read.table(BEP_list_t2[p,2])
    # interpolate BEP 2 corresponding to defined resolution
    int <- approx(BEP_t2[,1:2],xout = seq(BEP_t2[1,1], BEP_t2[nrow(BEP_t2),1],resolution))
    BEP_t2 <- data.frame(int$x,int$y)
      
    L_t1 <- BEP_t1[nrow(BEP_t1),1]-BEP_t1[1,1]
    L_t2 <- BEP_t2[nrow(BEP_t2),1]-BEP_t2[1,1]
    L_BEP <- mean(c(L_t1, L_t2))
    
    # for-loop over all iterations (MCS)
    for(i in 1:n_iterations){
      #print(i)
        
      # path to baselines
      base_path_t1 <- paste0("../2_ZC/results/",project,"/",simulation_title,"/", t1,"/",ID,"/baselines/")
      base_path_t2 <- paste0("../2_ZC/results/",project,"/",simulation_title,"/", t2,"/",ID,"/baselines/")
      
      # define break-off criteria
      check_1 <- list.files(base_path_t1, full.names = F)
      check_1 <- length(check_1[grepl(paste0("iteration_",i,".txt"), check_1)])-1
      check_2 <- list.files(base_path_t2, full.names = F)
      check_2 <- length(check_2[grepl(paste0("iteration_",i,".txt"), check_2)])-1
        
      if(check_1!=n_layers|check_1!=n_layers) next
        
      # Read baselines
      bl1t1 <- read.table(paste0(base_path_t1, "baseline_scale_1_iteration_", i, ".txt"), header = T)
      bl1t2 <- read.table(paste0(base_path_t2, "baseline_scale_1_iteration_", i, ".txt"), header = T)
          
      # Interpolate baselines according to defined resolution
      int <- approx(bl1t1[,1:2],xout = seq(bl1t1[1,1], bl1t1[nrow(bl1t1),1],resolution))
      bl1t1 <- data.frame(int$x,int$y)
      int <- approx(bl1t2[,1:2],xout = seq(bl1t2[1,1], bl1t2[nrow(bl1t2),1],resolution))
      bl1t2 <- data.frame(int$x,int$y)
        
      BEP_t1[ ,1] <- round(BEP_t1[ ,1], 3)
      BEP_t2[ ,1] <- round(BEP_t2[ ,1], 3)
      bl1t1[ ,1] <- round(bl1t1[ ,1], 3)
      bl1t2[ ,1] <- round(bl1t2[ ,1], 3)
      
      # if necessary, set both profiles to the same X values
      minus_1 <- nrow(BEP_t1)-nrow(bl1t1)
      if(minus_1>0) BEP_t1 <- BEP_t1[BEP_t1[ ,1]%in%bl1t1[ ,1], ]
      if(minus_1<0) bl1t1 <- bl1t1[bl1t1[ ,1]%in%BEP_t1[ ,1], ]
      minus_2 <- nrow(BEP_t2)-nrow(bl1t2)
      if(minus_2>0) BEP_t2 <- BEP_t2[BEP_t2[ ,1]%in%bl1t2[ ,1], ]
      if(minus_2<0) bl1t2 <- bl1t2[bl1t2[ ,1]%in%BEP_t2[ ,1], ]
      
      poly_points_t1 <- rbind(BEP_t1[ ,1:2], bl1t1[order(bl1t1[ ,1], decreasing = T),1:2], BEP_t1[1,1:2])
      poly_points_t2 <- rbind(BEP_t2[ ,1:2], bl1t2[order(bl1t2[ ,1], decreasing = T),1:2], BEP_t2[1,1:2])
      
      area_L1 <- mean(c(abs(polyarea(poly_points_t1[ ,1], poly_points_t1[ ,2])), abs(polyarea(poly_points_t2[ ,1], poly_points_t2[ ,2]))))
      
      # calculate bedfroms of layer 1 by detrending
      bedforms_t1_L1 <- BEP_t1[ ,2] - bl1t1[ ,2]
      bedforms_t2_L1 <- BEP_t2[ ,2] - bl1t2[ ,2]
      bedforms_t1_L1 <- data.frame(X=BEP_t1[ ,1], Z=bedforms_t1_L1)
      bedforms_t2_L1 <- data.frame(X=BEP_t2[ ,1], Z=bedforms_t2_L1)
      bedforms_t1_L1 <- bedforms_t1_L1[!is.na(bedforms_t1_L1$Z),]
      bedforms_t2_L1 <- bedforms_t2_L1[!is.na(bedforms_t2_L1$Z),]

      
      if(n_layers>1){
        
        # Read baselines
        bl2t1 <- read.table(paste0(base_path_t1, "baseline_scale_2_iteration_", i, ".txt"), header = T)
        bl2t2 <- read.table(paste0(base_path_t2, "baseline_scale_2_iteration_", i, ".txt"), header = T)
        
        # Interpolate baselines according to defined resolution
        int <- approx(bl2t1[,1:2],xout = seq(bl2t1[1,1], bl2t1[nrow(bl2t1),1],resolution))
        bl2t1 <- data.frame(int$x,int$y)
        int <- approx(bl2t2[,1:2],xout = seq(bl2t2[1,1], bl2t2[nrow(bl2t2),1],resolution))
        bl2t2 <- data.frame(int$x,int$y)
        
        bl2t1[ ,1] <- round(bl2t1[ ,1], 3)
        bl2t2[ ,1] <- round(bl2t2[ ,1], 3)
        
        # if necessary, set both to same X values
        minus_1 <- nrow(bl1t1)-nrow(bl2t1)
        if(minus_1>0) bl1t1 <- bl1t1[bl1t1[ ,1]%in%bl2t1[ ,1], ]
        if(minus_1<0) bl2t1 <- bl2t1[bl2t1[ ,1]%in%bl1t1[ ,1], ]
        minus_2 <- nrow(bl1t2)-nrow(bl2t2)
        if(minus_2>0) bl1t2 <- bl1t2[bl1t2[ ,1]%in%bl2t2[ ,1], ]
        if(minus_2<0) bl2t2 <- bl2t2[bl2t2[ ,1]%in%bl1t2[ ,1], ]
        
        # if necessary, set both to same X values
        minus_1 <- nrow(BEP_t1)-nrow(bl2t1)
        if(minus_1>0) BEP_t1 <- BEP_t1[BEP_t1[ ,1]%in%bl2t1[ ,1], ]
        if(minus_1<0) bl2t1 <- bl2t1[bl2t1[ ,1]%in%BEP_t1[ ,1], ]
        minus_2 <- nrow(bl1t2)-nrow(bl2t2)
        if(minus_2>0) BEP_t2 <- BEP_t2[BEP_t2[ ,1]%in%bl2t2[ ,1], ]
        if(minus_2<0) bl2t2 <- bl2t2[bl2t2[ ,1]%in%BEP_t2[ ,1], ]
        
        poly_points_t1 <- rbind(bl1t1[ ,1:2], bl2t1[order(bl2t1[ ,1], decreasing = T),1:2], bl1t1[1,1:2])
        poly_points_t2 <- rbind(bl1t2[ ,1:2], bl2t2[order(bl2t2[ ,1], decreasing = T),1:2], bl1t2[1,1:2])
        
        area_L2 <- mean(c(abs(polyarea(poly_points_t1[ ,1], poly_points_t1[ ,2])), abs(polyarea(poly_points_t2[ ,1], poly_points_t2[ ,2]))))
      }
      
      # calculate cross correlation for BEPs t1 and t2      
      ccf_BEP <- ccf(BEP_t2[ ,2], BEP_t1[ ,2], lag.max=lag_max, plot=F)
      ccf_BEP <- data.frame(lag = ccf_BEP$lag, corr = ccf_BEP$acf, lag.max = lag_max)
      max_lag_BEP <- data.frame(ccf_BEP[which.max(ccf_BEP$corr), ])
      v_BEP <- round(max_lag_BEP$lag*resolution/delta_t,2)
      corr_BEP <- data.frame(m, BEP=ID, iteration=i, t1=t1, t2=t2, dt=delta_t, 
                             lag=max_lag_BEP$lag*resolution, v_BEP, corr=max_lag_BEP$corr)
      migration_BEP <- rbind(migration_BEP, corr_BEP)
        
      # calculate cross correlation for bedforms layer 1 t1 and t2    
      ccf_L_1 <- ccf(bedforms_t2_L1[ ,2], bedforms_t1_L1[ ,2], lag.max=lag_max, plot=F)
      ccf_L_1 <- data.frame(lag = ccf_L_1$lag, corr = ccf_L_1$acf, lag.max = lag_max)
      max_lag_L1 <- data.frame(ccf_L_1[which.max(ccf_L_1$corr), ])
      v_L1 <- round(max_lag_L1$lag*resolution/delta_t,2)
      bedload_L1 <- (v_L1/3600)*(area_L1)*(1/L_BEP)*(1-porosity)*density*1000
      corr_L1 <- data.frame(m, BEP=ID, iteration=i, t1=t1, t2=t2, dt=delta_t, 
                            lag=max_lag_L1$lag*resolution, v_L1, bedload_L1, corr=max_lag_L1$corr)
      migration_L1 <- rbind(migration_L1, corr_L1)
        
      if(n_layers>1){
        # calculate cross correlation for bedforms layer 2 t1 and t2
        ccf_L_2 <- ccf(bl1t2[,2], bl1t1[ ,2], lag.max=lag_max, plot=F)
        ccf_L_2 <- data.frame(lag = ccf_L_2$lag, corr = ccf_L_2$acf, lag.max = lag_max)
        max_lag_L2 <- data.frame(ccf_L_2[which.max(ccf_L_2$corr), ])
        v_L2 <- round(max_lag_L2$lag*resolution/delta_t,2)
        bedload_L2 <- v_L2/3600*area_L2*1/L_BEP*(1-porosity)*density*1000
        corr_L2 <- data.frame(m, BEP=ID, iteration=i, t1=t1, t2=t2, dt=delta_t, 
                              lag=max_lag_L2$lag*resolution, v_L2, bedload_L2, corr=max_lag_L2$corr)
        migration_L2 <- rbind(migration_L2, corr_L2)
      }
    }
  }
}
  
# export results
dir.create(paste0("results/"))
dir.create(paste0("results/",project))
dir.create(paste0("results/",project,"/",simulation_title))
cc_results <- list(migration_BEP,migration_L1,migration_L2)
saveRDS(cc_results, paste0("results/",project,"/",simulation_title,"/cc_results.rds"))

