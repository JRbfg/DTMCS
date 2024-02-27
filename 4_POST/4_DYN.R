### set parameters
project <- "Emmerich_detail"
simulation_title <- "run_1"
n_layers <- 2
n_iterations <- 100
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

# load results
cc_results <- readRDS(paste0("../3a_CC/results/", project, "/", simulation_title, "/cc_results.rds"))
ca_results <- readRDS(paste0("../3b_CA/results/", project, "/", simulation_title, "/ca_results.rds"))

cc_BEP <- cc_results[[1]]
cc_l1 <- cc_results[[2]]
cc_l2 <- cc_results[[3]]

# filter by correlation
cc_BEP <- cc_BEP[cc_BEP$corr>0.5&cc_BEP$lag>0, ]
cc_l1 <- cc_l1[cc_l1$corr>0.5&cc_l1$lag>0, ]
cc_l2 <- cc_l2[cc_l2$corr>0.5&cc_l2$lag>0, ]

ca_l1 <- ca_results[[1]]
ca_l2 <- ca_results[[2]]
ca_l1 <- ca_l1[ca_l1$lag>0, ]


IDs <- unique(cc_BEP$BEP)

# create empty dataframe to store results
bf_dyn <- data.frame()

# for-loop over all meauserements
for(m in 1:nrow(measurement_table)){
  print(paste0("measurement pair: ",m))
  
  BEPs_1 <- BEP_list[BEP_list$folder==measurement_table[m,1], ]
  BEPs_2 <- BEP_list[BEP_list$folder==measurement_table[m,2], ]
  delta_t <- measurement_table[m,3]
  
  for(ID in IDs){

    cc_lag_BEP <- median(cc_BEP[cc_BEP$m==m&cc_BEP$BEP==ID,7], na.rm=T)
    cc_mig_BEP <- median(cc_BEP[cc_BEP$m==m&cc_BEP$BEP==ID,8], na.rm=T)

    cc_lag_l1 <- median(cc_l1[cc_l1$m==m&cc_l1$BEP==ID,7], na.rm=T)
    cc_mig_l1 <- median(cc_l1[cc_l1$m==m&cc_l1$BEP==ID,8], na.rm=T)
    cc_bl_l1 <- median(cc_l1[cc_l1$m==m&cc_l1$BEP==ID,9], na.rm=T)
    
    cc_lag_l2 <- median(cc_l2[cc_l2$m==m&cc_l2$BEP==ID,7], na.rm=T)
    cc_mig_l2 <- median(cc_l2[cc_l2$m==m&cc_l2$BEP==ID,8], na.rm=T)
    cc_bl_l2 <- median(cc_l2[cc_l2$m==m&cc_l2$BEP==ID,9], na.rm=T)
    
    
    ca_lag_l1 <- median(ca_l1[ca_l1$m==m&ca_l1$BEP==ID,7], na.rm=T)
    ca_n_l1 <- median(ca_l1[ca_l1$m==m&ca_l1$BEP==ID,10], na.rm=T)
    ca_mig_l1 <- median(ca_l1[ca_l1$m==m&ca_l1$BEP==ID,8], na.rm=T)
    ca_bl_l1 <- median(ca_l1[ca_l1$m==m&ca_l1$BEP==ID,9], na.rm=T)

    ca_lag_l2 <- median(ca_l2[ca_l2$m==m&ca_l2$BEP==ID,7], na.rm=T)
    ca_n_l2 <- median(ca_l2[ca_l2$m==m&ca_l2$BEP==ID,10], na.rm=T)
    ca_mig_l2 <- median(ca_l2[ca_l2$m==m&ca_l2$BEP==ID,8], na.rm=T)
    ca_bl_l2 <- median(ca_l2[ca_l2$m==m&ca_l2$BEP==ID,9], na.rm=T)
    
    bf_dyn <- rbind(bf_dyn, data.frame(m, BEP=ID, t1=measurement_table[m,1], t2=measurement_table[m,2], dt=measurement_table[m,3],
                                       cc_lag_BEP, cc_mig_BEP, cc_lag_l1, cc_mig_l1, cc_bl_l1, cc_lag_l2, cc_mig_l2, cc_bl_l2,
                                       ca_lag_l1, ca_n_l1, ca_mig_l1, ca_bl_l1, ca_lag_l2, ca_n_l2, ca_mig_l2, ca_bl_l2))
  }
}


### create plots

par(mar=c(5,5,3,5))
plot(NULL, xlim=c(1,45), ylim=c(0,2.5), xlab="measurement pair", ylab="migraion rate [m/h]", main=" bedform migration cc")
abline(v=seq(0,50), h=seq(0,3,0.5), col="grey", lty=2)
points(cc_l1[cc_l1$corr>0.5,c(1,8)], pch=4, cex=1, col="red")
points(cc_l2[cc_l2$corr>0.5,c(1,8)], pch=4, cex=1, col="blue")
points(bf_dyn$m, bf_dyn$cc_mig_BEP, pch=19, cex=2)
points(bf_dyn$m, bf_dyn$cc_mig_l1, pch=19, cex=2, col="red")
points(bf_dyn$m, bf_dyn$cc_mig_l2, pch=19, cex=2, col="blue")
par(new=T)
plot(NULL, xlim=c(1,45), ylim=c(0,50), xlab="", ylab="", xaxt="none", yaxt="none")
axis(4)
lines(bf_dyn$m, bf_dyn$dt, col="darkgreen")
mtext("dt [h]", side=4, line=3)


par(mar=c(5,5,3,5))
plot(NULL, xlim=c(1,45), ylim=c(0,5), xlab="measurement pair", ylab="migraion rate [m/h]", main=" bedform migration ca")
abline(v=seq(0,50), h=seq(0,5,0.5), col="grey", lty=2)
points(ca_l1[ ,c(1,8)], pch=4, cex=1, col="red")
points(ca_l2[ ,c(1,8)], pch=4, cex=1, col="blue")
points(bf_dyn$m, bf_dyn$ca_mig_l1, pch=19, cex=2, col="red")
points(bf_dyn$m, bf_dyn$ca_mig_l2, pch=19, cex=2, col="blue")
par(new=T)
plot(NULL, xlim=c(1,45), ylim=c(0,50), xlab="", ylab="", xaxt="none", yaxt="none")
axis(4)
lines(bf_dyn$m, bf_dyn$dt, col="darkgreen")
mtext("dt [h]", side=4, line=3)


par(mar=c(5,5,3,5))
plot(NULL, xlim=c(1,50), ylim=c(0,1), xlab="dt [h]", ylab="correlation [-]", main="cross correlation")
abline(v=seq(0,50), h=seq(0,1,0.2), col="grey", lty=2)
points(cc_results[[2]][ ,c(6,10)], pch=4, cex=1, col="red")
points(cc_results[[3]][ ,c(6,10)], pch=4, cex=1, col="blue")
points(cc_results[[1]][ ,c(6,9)], pch=4, cex=1, col="black")


par(mar=c(5,5,3,5))
plot(NULL, xlim=c(1,45), ylim=c(0,70), xlab="measurement pair", ylab="number of bedforms [-]", main="")
abline(v=seq(0,50), h=seq(0,100,5), col="grey", lty=2)
points(ca_l1$m, ca_l1$count, col="red", pch=4, cex=0.5)
points(ca_l2$m, ca_l2$count, col="blue", pch=4, cex=0.5)
points(bf_dyn$m, bf_dyn$ca_n_l1, pch=19, cex=2, col="red")
points(bf_dyn$m, bf_dyn$ca_n_l2, pch=19, cex=2, col="blue")


max_cc_bl1 <- quantile(bf_dyn$cc_bl_l1[bf_dyn$dt<2], probs=0.95, na.rm=T)
max_ca_bl1 <- quantile(bf_dyn$ca_bl_l1[bf_dyn$dt<2], probs=0.95, na.rm=T)
max_cc_bl2 <- quantile(bf_dyn$cc_bl_l2[bf_dyn$dt>5], probs=0.95, na.rm=T)
max_ca_bl2 <- quantile(bf_dyn$ca_bl_l2[bf_dyn$dt>5], probs=0.95, na.rm=T)
min_cc_bl1 <- quantile(bf_dyn$cc_bl_l1[bf_dyn$dt<2], probs=0.05, na.rm=T)
min_ca_bl1 <- quantile(bf_dyn$ca_bl_l1[bf_dyn$dt<2], probs=0.05, na.rm=T)
min_cc_bl2 <- quantile(bf_dyn$cc_bl_l2[bf_dyn$dt>5], probs=0.05, na.rm=T)
min_ca_bl2 <- quantile(bf_dyn$ca_bl_l2[bf_dyn$dt>5], probs=0.05, na.rm=T)
med_cc_bl1 <- quantile(bf_dyn$cc_bl_l1[bf_dyn$dt<2], probs=0.5, na.rm=T)
med_ca_bl1 <- quantile(bf_dyn$ca_bl_l1[bf_dyn$dt<2], probs=0.5, na.rm=T)
med_cc_bl2 <- quantile(bf_dyn$cc_bl_l2[bf_dyn$dt>5], probs=0.5, na.rm=T)
med_ca_bl2 <- quantile(bf_dyn$ca_bl_l2[bf_dyn$dt>5], probs=0.5, na.rm=T)

max_cc_total <- max_cc_bl1+max_cc_bl2
max_ca_total <- max_ca_bl1+max_ca_bl2
min_cc_total <- min_cc_bl1+min_cc_bl2
min_ca_total <- min_ca_bl1+min_ca_bl2
med_cc_total <- med_cc_bl1+med_cc_bl2
med_ca_total <- med_ca_bl1+med_ca_bl2

cc_df <- data.frame(med_cc_bl1, min_cc_bl1, max_cc_bl1, med_cc_bl2, min_cc_bl2, max_cc_bl2, med_cc_total, min_cc_total, max_cc_total)
ca_df <- data.frame(med_ca_bl1, min_ca_bl1, max_ca_bl1, med_ca_bl2, min_ca_bl2, max_ca_bl2, med_ca_total, min_ca_total, max_ca_total)

par(oma = c(0,0,2,0))

par(mar=c(6,5,3,3))
plot(NULL, xlim=c(0.5,2.5), ylim=c(0,180), xaxt="none", xlab="", ylab="", cex.axis=1.3, cex.lab=1.5)
abline(h=seq(0,200,25), lty=2, col="lightgrey")
mtext(side=2, text="bedload transport rate [g/s*m]", cex=1.5, line=3.5)
axis(1, at=c(1:2), labels = c( "cross correlation analysis", "centroid analysis"), las=1, cex.axis=1.5, cex.lab=1.5, mgp=c(1,1,0))
polygon(x=c(0.7,0.9,0.9,0.7), y=c(cc_df[8],cc_df[8], cc_df[9], cc_df[9]), col="grey", lty=2)
lines(x=c(0.7,0.9), y=rep(cc_df[7],2), lwd=2)
polygon(x=c(1.7,1.9,1.9,1.7), y=c(ca_df[8],ca_df[8], ca_df[9], ca_df[9]), col="grey", lty=2)
lines(x=c(1.7,1.9), y=rep(ca_df[7],2), lwd=2)
polygon(x=c(0.9,1.1,1.1,0.9), y=c(cc_df[2],cc_df[2], cc_df[3], cc_df[3]), col="red", lty=2)
lines(x=c(0.9,1.1), y=rep(cc_df[1],2), lwd=2)
polygon(x=c(1.9,2.1,2.1,1.9), y=c(ca_df[2],ca_df[2], ca_df[3], ca_df[3]), col="red", lty=2)
lines(x=c(1.9,2.1), y=rep(ca_df[1],2), lwd=2)
polygon(x=c(1.1,1.3,1.3,1.1), y=c(cc_df[5],cc_df[5], cc_df[6], cc_df[6]), col="blue", lty=2)
lines(x=c(1.1,1.3), y=rep(cc_df[4],2), lwd=2)
polygon(x=c(2.1,2.3,2.3,2.1), y=c(ca_df[5],ca_df[5], ca_df[6], ca_df[6]), col="blue", lty=2)
lines(x=c(2.1,2.3), y=rep(ca_df[4],2), lwd=2)

