### set parameters
project <- "Emmerich_TH"
simulation_title <- "run_1"
n_layers <- 2
n_iterations <- 100
# define subsections [km]
start <- 860.0
end <- 860.5
l_section <- 0.1
###
sections <- seq(start, end, l_section)

# change working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# list all BEP-files in project folder
BEP_list <- data.frame(folder=sub("/.*", "", list.files(paste0("../0_DATA/", project), full.names = F, pattern=".txt", recursive = T)), 
                       file=list.files(paste0("../0_DATA/", project), full.names = T, pattern=".txt", recursive = T), stringsAsFactors = F)
measurements <- unique(BEP_list$folder)

# create empty dataframe to store results
bf_geom <- data.frame()

# for-loop over all meauserements
for(m in measurements[1]){
  print(m)
  BEP_list_d <- BEP_list[BEP_list$folder==m, ]
  BEP_list_d <- BEP_list_d[c(15,16,2:6), ]
  for(p in 1:nrow(BEP_list_d)){
    
    BEP <- read.table(BEP_list_d[p,2])
    ID <- unique(BEP[ ,4])
    print(ID)
    # for-loop over all iterations (MCS)
    for(it in 1:n_iterations){
      bl <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",m,"/",ID,"/baselines/baseline_scale_1_iteration_",it,".txt"), header = T)
      if(n_layers==2) bl <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",m,"/",ID,"/baselines/baseline_scale_2_iteration_",it,".txt"), header = T)
      bl_int <- data.frame(approx(bl[ ,1:2], xout = BEP[ ,1]))
      
      stats_l1 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",m,"/",ID,"/statistics/statistics_scale_1_iteration_",it,".txt"), header = T)
      if(n_layers==2) stats_l2 <- read.table(paste0("../2_ZC/results/",project,"/",simulation_title,"/",m,"/",ID,"/statistics/statistics_scale_2_iteration_",it,".txt"), header = T)
      
      bedforms <- data.frame()
      for(j in 2:nrow(bl)){
        bf <- BEP[BEP[ ,1]>=bl[j-1,1]&BEP[ ,1]<bl[j,1], ]
        if(nrow(bf)<3) next
        trend <- approx(bf[c(1,nrow(bf)),1:2], xout = bf[ ,1])
        bf[ ,2] <- bf[ ,2]-trend$y
        h_total <- max(bf[ ,2], na.rm=T)
        l_total <- bf[nrow(bf),1]-bf[1,1,]
        center <- mean(bf[ ,3], na.rm=T)
        bedforms <- rbind(bedforms, data.frame(h_total, l_total, center))
      }
      
      for(s in 2:length(sections)){
        BEP_s <- BEP[BEP[ ,3]>=sections[s-1]&BEP[ ,3]<sections[s], ]
        bl_s <- bl_int[bl_int$x%in%BEP_s[ ,1], ]
        T90 <- as.numeric(quantile((BEP_s[ ,2]-bl_s$y), probs=0.9, na.rm=T))
        
        stats_l1_s <- stats_l1[stats_l1$x_pos>=sections[s-1]&stats_l1$x_pos<sections[s], ]
        stats_l2_s <- stats_l2[stats_l2$x_pos>=sections[s-1]&stats_l2$x_pos<sections[s], ]
        bedforms_s <- bedforms[bedforms$center>=sections[s-1]&bedforms$center<sections[s], ]
        
        l_total <- mean(bedforms_s$l_total, na.rm=T)
        h_total <- mean(bedforms_s$h_total, na.rm=T)
        l_1 <- mean(stats_l1_s$l, na.rm=T)
        h_1 <- mean(stats_l1_s$h, na.rm=T)
        l_2 <- mean(stats_l2_s$l, na.rm=T)
        h_2 <- mean(stats_l2_s$h, na.rm=T)
        
                
        params_s <- data.frame(m, ID, it, section=sections[s-1], T90, l_total, h_total, l_1, h_1, l_2, h_2)
        bf_geom <- rbind(bf_geom, params_s)
      }
    }
  }
}
for(ID in unique(bf_geom$ID)){
  for(i in 1:n_iterations){
    bf_geom$h_total[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$h_total[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$l_total[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$l_total[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$T90[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$T90[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$h_1[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$h_1[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$h_2[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$h_2[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$l_1[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$l_1[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
    bf_geom$l_2[bf_geom$ID==ID&bf_geom$it==i] <- mean(bf_geom$l_2[bf_geom$ID==ID&bf_geom$it==i], na.rm=T)
  }
}


### create plots

# H_total
par(mar=c(5,5,3,3))
plot(NULL, xlim=c(8,14), ylim=c(0,0.5), xlab="BEP", ylab="H_total [m]")
abline(v=seq(1,20,1), h=seq(0,1,0.1), lty=2, col="grey")
df <- data.frame()
for(ID in unique(bf_geom$ID)){
  max <- max(bf_geom$h_total[bf_geom$ID==ID], na.rm=T)
  min <- min(bf_geom$h_total[bf_geom$ID==ID], na.rm=T)
  med <- median(bf_geom$h_total[bf_geom$ID==ID], na.rm=T)
  df <- rbind(df, data.frame(ID, min, max, med))
}
df <- df[order(df$ID), ]
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min, rev(df$max)), col="grey", border=NA)
lines(df$ID, df$med)


# T90
par(mar=c(5,5,3,3))
plot(NULL, xlim=c(8,14), ylim=c(0,0.5), xlab="BEP", ylab="T90 [m]")
abline(v=seq(1,20,1), h=seq(0,1,0.1), lty=2, col="grey")
df <- data.frame()
for(ID in unique(bf_geom$ID)){
  max <- max(bf_geom$T90[bf_geom$ID==ID], na.rm=T)
  min <- min(bf_geom$T90[bf_geom$ID==ID], na.rm=T)
  med <- median(bf_geom$T90[bf_geom$ID==ID], na.rm=T)
  df <- rbind(df, data.frame(ID, min, max, med))
}
df <- df[order(df$ID), ]
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min, rev(df$max)), col="grey", border=NA)
lines(df$ID, df$med)


# H_1,2
par(mar=c(5,5,3,3))
plot(NULL, xlim=c(8,14), ylim=c(0,0.5), xlab="BEP", ylab="H_1,2 [m]")
abline(v=seq(1,20,1), h=seq(0,1,0.1), lty=2, col="grey")
df <- data.frame()
for(ID in unique(bf_geom$ID)){
  max1 <- max(bf_geom$h_1[bf_geom$ID==ID], na.rm=T)
  min1 <- min(bf_geom$h_1[bf_geom$ID==ID], na.rm=T)
  med1 <- median(bf_geom$h_1[bf_geom$ID==ID], na.rm=T)
  max2 <- max(bf_geom$h_2[bf_geom$ID==ID], na.rm=T)
  min2 <- min(bf_geom$h_2[bf_geom$ID==ID], na.rm=T)
  med2 <- median(bf_geom$h_2[bf_geom$ID==ID], na.rm=T)
  df <- rbind(df, data.frame(ID, min1, max1, med1, min2, max2, med2))
}
df <- df[order(df$ID), ]
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min1, rev(df$max1)), col="red", border=NA)
lines(df$ID, df$med1)
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min2, rev(df$max2)), col="blue", border=NA)
lines(df$ID, df$med2)


# L_1,2
par(mar=c(5,5,3,3))
plot(NULL, xlim=c(8,14), ylim=c(0,30), xlab="BEP", ylab="L_1,2 [m]")
abline(v=seq(1,20,1), h=seq(0,100,5), lty=2, col="grey")
df <- data.frame()
for(ID in unique(bf_geom$ID)){
  max1 <- max(bf_geom$l_1[bf_geom$ID==ID], na.rm=T)
  min1 <- min(bf_geom$l_1[bf_geom$ID==ID], na.rm=T)
  med1 <- median(bf_geom$l_1[bf_geom$ID==ID], na.rm=T)
  max2 <- max(bf_geom$l_2[bf_geom$ID==ID], na.rm=T)
  min2 <- min(bf_geom$l_2[bf_geom$ID==ID], na.rm=T)
  med2 <- median(bf_geom$l_2[bf_geom$ID==ID], na.rm=T)
  df <- rbind(df, data.frame(ID, min1, max1, med1, min2, max2, med2))
}
df <- df[order(df$ID), ]
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min1, rev(df$max1)), col="red", border=NA)
lines(df$ID, df$med1)
polygon(x=c(df$ID, rev(df$ID)), y=c(df$min2, rev(df$max2)), col="blue", border=NA)
lines(df$ID, df$med2)

