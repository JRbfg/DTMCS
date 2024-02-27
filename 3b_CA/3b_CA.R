### set parameters
project <- "Emmerich_detail"
simulation_title <- "run_1"
n_layers <- 2
n_iterations <- 100
# define porosity and grain density for bedload estimation
porosity <- 0.3
density <- 2603
###

# load required libraries
library(pracma)
library(sf)
library(rgeos)
library(dplyr)

# change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# read table with information about time offsets between measurements
measurement_table <- read.table(paste0("../0_DATA/",project,"/measurements.csv"), header = T, sep=";")

# list all BEP-files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                       file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)

# create empty dataframes to store the results
ca_results_L1 <- data.frame()
ca_results_L2 <- data.frame()
ca_means_L1 <- data.frame()
ca_means_L2 <- data.frame()

# for-loop over all meauserements
for(m in 1:nrow(measurement_table)){
  print(paste0("measurement pair: ",m))
  
  BEPs_1 <- BEP_list[BEP_list$folder==measurement_table[m,1], ]
  BEPs_2 <- BEP_list[BEP_list$folder==measurement_table[m,2], ]
  delta_t <- measurement_table[m,3]
  
    for(p in 1:nrow(BEPs_1)){
      # for-loop over all layers
      for(layer in 1:n_layers){
        # for-loop over all iterations (MCS)
        for(it in 1:n_iterations){
          
          BEP_1 <- read.table(BEPs_1[p,2])
          BEP_2 <- read.table(BEPs_2[p,2])
          ID <- unique(BEP_1[ ,4])
          BEP_1 <- BEP_1[ ,c(1,2,3)]
          BEP_2 <- BEP_2[ ,c(1,2,3)]
          
          check_dirs <- list.dirs(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,1]), recursive = F, full.names = F)
          if(length(check_dirs)==0) print(paste0("WARNING: no zc-results available for measurement ", measurement_table[m,1]))
          
          if(layer==1){
          bl1 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,1],"/",ID,"/baselines/baseline_scale_1_iteration_",it,".txt"), header = T)
          bl2 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,2],"/",ID,"/baselines/baseline_scale_1_iteration_",it,".txt"), header = T)
          }
          
          # in case of layer 2 baseline of layer 1 acts as BEP 
          if(layer==2){
            BEP_1 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,1],"/",ID,"/baselines/baseline_scale_1_iteration_",it,".txt"), header = T)
            BEP_2 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,2],"/",ID,"/baselines/baseline_scale_1_iteration_",it,".txt"), header = T)
            bl1 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,1],"/",ID,"/baselines/baseline_scale_2_iteration_",it,".txt"), header = T)
            bl2 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",measurement_table[m,2],"/",ID,"/baselines/baseline_scale_2_iteration_",it,".txt"), header = T)
          }
          
          # calculate bedform attributes for t1
          bedforms_1 <- data.frame()
          bf_params <- data.frame()
          
          for(j in 2:nrow(bl1)){
            sel <- BEP_1[BEP_1[ ,1]>=bl1[j-1,1]&BEP_1[ ,1]<bl1[j, 1], ]
            if(nrow(sel)<3)next
            p_max <- sel[which.max(sel[ ,2]), ]
            poly_points <- (rbind(sel[,1:2],sel[1,1:2]))
            poly <- Polygon(poly_points)
            poly <- sp::Polygons(list(poly), ID = "A")
            poly <- sp::SpatialPolygons(list(poly))
            center <- gCentroid((poly))
            center_x <- center@coords[1]
            center_y <- center@coords[2]
            area <- abs(polyarea(poly_points[ ,1], poly_points[ ,2]))
            l <- max(sel[ ,1])-min(sel[ ,1])
            center_km <- mean(sel[ ,3])
            bedform <- data.frame(x=sel[,1], y=sel[ ,2], ID=j-1,  h=p_max[ ,2]-mean(c(bl1[j,2], bl1[j-1,2])))
            bedforms_1 <- rbind(bedforms_1, bedform)
            bf_params <- rbind(bf_params, data.frame(p, j=j-1, t=1, h=p_max[ ,2]-mean(c(bl1[j,2], bl1[j-1,2])), area=area, center_x, center_y, center_km, l))
          }
          
          #exclude bedforms smaller than 5 cm
          bf_params <- bf_params[bf_params$h>=0.05, ]
          
          # calculate bedform attributes for t1
          bedforms_2 <- data.frame()
          
          for(j in 2:nrow(bl2)){
            sel <- BEP_2[BEP_2[ ,1]>=bl2[j-1,1]&BEP_2[ ,1]<bl2[j, 1], ]
            if(nrow(sel)<3) next
            p_max <- sel[which.max(sel[ ,2]), ]
            poly_points <- (rbind(sel[,1:2],sel[1,1:2]))
            poly <- Polygon(poly_points)
            poly <- sp::Polygons(list(poly), ID = "A")
            poly <- sp::SpatialPolygons(list(poly))
            center <- gCentroid((poly))
            center_x <- center@coords[1]
            center_y <- center@coords[2]
            area <- abs(polyarea(poly_points[ ,1], poly_points[ ,2]))
            center_km <- mean(sel[ ,3])
            l <- max(sel[ ,1])-min(sel[ ,1])
            dune <- data.frame(x=sel[,1], y=sel[ ,2], ID=j-1,  h=p_max[ ,2]-mean(c(bl2[j,2], bl2[j-1,2])))
            bedforms_2 <- rbind(bedforms_2, dune)
            bf_params <- rbind(bf_params, data.frame(p, j=j-1, t=2, h=p_max[ ,2]-mean(c(bl2[j,2], bl2[j-1,2])), area=area, center_x, center_y, center_km, l))
          }
          
          #exclude bedforms smaller than 5 cm
          bf_params <- bf_params[bf_params$h>=0.05, ]
          bedforms_1 <- bedforms_1[bedforms_1$h>=0.05, ]
          bedforms_2 <- bedforms_2[bedforms_2$h>=0.05, ]
          
          # reset bedform IDs
          for(i in unique(bf_params$t)){
            bf_params$j[bf_params$t==i] <- c(1:nrow(bf_params[bf_params$t==i, ]))
          }
          
          # identify corresponding bedforms (minimum downstream distance)
          for(j in unique(bf_params$j)){
            check_1 <- bf_params$center_x[bf_params$j==j&bf_params$t==1]
            if(length(check_1)==0) next
            
            bf_params$help <- NA
            bf_params$help[bf_params$t==2] <- abs(bf_params$center_x[bf_params$t==2] - bf_params$center_x[bf_params$t==1&bf_params$j==j])
            bf_params$j[which.min(bf_params$help)] <- j
          }
          
          # calculate attributes of corresponding bedforms
          
          # create dataframe to store interim results
          ca_temp <- data.frame()
          for(j in unique(bf_params$j)){
            check_1 <- bf_params$center_x[bf_params$j==j&bf_params$t==1]
            check_2 <- bf_params$center_x[bf_params$j==j&bf_params$t==2]
            
            if(length(check_1)==0|length(check_2)==0) next
            
            lag <- bf_params$center_x[bf_params$j==j&bf_params$t==2] - bf_params$center_x[bf_params$j==j&bf_params$t==1]
            L_ratio <- bf_params$l[bf_params$j==j&bf_params$t==2]/bf_params$l[bf_params$j==j&bf_params$t==1]
            area_ratio <- bf_params$area[bf_params$j==j&bf_params$t==2]/bf_params$area[bf_params$j==j&bf_params$t==1]
            area_mean <- mean(bf_params$area[bf_params$j==j])
            L_mean <- mean(bf_params$l[bf_params$j==j])
            xpos <-  mean(bf_params$center_x[bf_params$j==j])
            rate <- lag/delta_t
            L_ges <- mean(c(bl1[nrow(bl1),1]-bl1[1,1], bl2[nrow(bl2),1]-bl2[1,1]))
            km <- mean(bf_params$center_km[bf_params$j==j])
            H_mean <- mean(bf_params$h[bf_params$j==j])
            ca_j <- data.frame(m, BEP=ID, iteration=it, dt=delta_t, layer=layer, j, lag, migration=rate, L_ratio, area_ratio, L_mean, H_mean, area_mean, xpos, L_ges, km)
            ca_temp <- rbind(ca_temp, ca_j)
          }
          
          # filter similarity of corresponding bedforms
          ca_temp <- ca_temp[ca_temp$lag>=0&abs(ca_temp$L_ratio-1)<=0.25&abs(ca_temp$area_ratio-1)<=0.25, ]
            
          # remove outliers
          sd <- sd(ca_temp$migration, na.rm=T)
          med <- median(ca_temp$migration, na.rm=T)
          check <- quantile(ca_temp$migration, 0.95, na.rm=T)
          if(nrow(ca_temp)>1) if(sd>3*med) check <- sd
          ca_temp <- ca_temp[ca_temp$migration<check, ]
   
          # store results
          if(layer==1) ca_results_L1 <- rbind(ca_results_L1, ca_temp)
          if(layer==2) ca_results_L2 <- rbind(ca_results_L2, ca_temp)
        
          # calculate auxiliary variables
          ca_temp$h1 <- ca_temp$migration*ca_temp$area_mean
          ca_temp$h2 <- ca_temp$migration*ca_temp$L_mean
          ca_temp$h3 <- ca_temp$lag*ca_temp$L_mean
            
          # calculate weighted average results
          bedload <- sum(ca_temp$h1, na.rm=T)*density*1000*(1-porosity)*(1/3600)*(1/sum(ca_temp$L_mean, na.rm=T))
          migration <- sum(ca_temp$h2, na.rm=T)*(1/sum(ca_temp$L_mean, na.rm=T))
          lag <- sum(ca_temp$h3, na.rm=T)*(1/sum(ca_temp$L_mean, na.rm = T))
          ca_mean_temp <- data.frame(m, BEP=ID, iteration=it, t1=measurement_table[m,1], t2=measurement_table[m,2], dt=delta_t, lag, migration, bedload, count=nrow(ca_temp))
            
          # store weighted average results
          if(layer==1) ca_means_L1 <- rbind(ca_means_L1, ca_mean_temp)
          if(layer==2) ca_means_L2 <- rbind(ca_means_L2, ca_mean_temp)
        }
      }
  }
}


# export results
dir.create("results")
dir.create(paste0("results/", project))
dir.create(paste0("results/", project, "/", simulation_title))
saveRDS(list(ca_means_L1, ca_means_L2, ca_results_L1, ca_results_L2), paste0("results/", project, "/", simulation_title, "/ca_results.rds"))

