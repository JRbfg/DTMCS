# add basepoints layer 1

add_basepoints_l1 <- function(intersections_2, profile_plus_baseline, baseline_extended){
  
  bp_add <- data.frame()
  colnames(intersections_2) <- c("X","Z")
  baseline_extended <- baseline_extended[!duplicated(baseline_extended$X), ]
  baseline_extended$X <- round(baseline_extended$X,3)
  baseline_extended <- baseline_extended[order(baseline_extended$X),]
  intersections_2[,1] <- round(intersections_2[,1],3)  
  
  for(i in 1:(nrow(baseline_extended)-1)){
    
    intersections_i <- intersections_2[intersections_2$X>baseline_extended[i,1]&intersections_2$X<baseline_extended[i+1,1],]
    if(nrow(intersections_i)==0&&(i!=nrow(baseline_extended)-1)) next
    
    # x-values
    x1 <- baseline_extended[i,1]
    # y-values
    y1 <- baseline_extended[i,2]
    # slope
    m1 <- (baseline_extended[i+1,2] - baseline_extended[i,2])/(baseline_extended[i+1,1] - baseline_extended[i,1])
    # intersection with y-axis
    b1 <- y1 - m1 * x1
    
    filter2 <- profile_plus_baseline[which(profile_plus_baseline[,1]>baseline_extended[i,1]&profile_plus_baseline[,1]<baseline_extended[i+1,1]),]
    
    # calculate the distance of each selected profile point to straight
    # profile point with largest distance equals the local minimum
    distances <- data.frame()
    for (j in 1:nrow(filter2)) {
      # x-values
      x2 <- filter2[j, 1]
      # y-values of baseline1/2
      y2 <- filter2[j, 2]
      # y-value of profile/baseline1
      y3 <- filter2[j, 3]
      # slope
      m2 <- (-1)/m1
      # intersection with y-axis
      b2 <- y2 - m2 * x2
      
      # intersection with connecting straight
      x_intersect <- (b2 - b1)/(m1 - m2)
      y_intersect <- m1 * x_intersect + b1
      
      dx <- x_intersect - x2
      dy <- y_intersect - y2
      
      dist <- abs(sqrt(dx^2 + dy^2))
      temp <- data.frame(x2, y2, y3, dist)
      distances <- rbind(distances, temp)
    }
    distances <- distances[which((distances$y3-distances$y2)>0),]
    if(nrow(distances)==0) next
    select_order <- distances[order(distances$dist), ]
    min_point <- select_order[nrow(select_order), 1:3]
    minimum <- min_point[, 1:2]
    colnames(minimum) <- colnames(intersections_2)
    bp_add <- rbind(bp_add, minimum)
    
    if(nrow(minimum)==0)next
    
    colnames(bp_add) <- colnames(baseline_extended)  
    #new_baseline <- baseline_extended[!duplicated(baseline_extended$X), ]  
    new_baseline <- rbind(baseline_extended[i,],baseline_extended[i+1,], bp_add)
    new_baseline <- new_baseline[order(new_baseline$X), ]
    new_baseline_int <- approx(new_baseline[1:nrow(new_baseline), ], xout = filter2[,1])
    new_baseline_int <- data.frame(X=new_baseline_int$x, Z=new_baseline_int$y)
    new_baseline_int <- new_baseline_int[new_baseline_int[,1]>filter2[1,1],]
    base_point <- new_baseline[new_baseline$X>minimum$x, ][1,]
    filter2b <- filter2[filter2[,1]>=new_baseline_int[1,1], ]
    # recalculate intersections for newly added basepoint with next intersection point
    # if there are still intersections another basepoint has to be added
    colnames(minimum) <- colnames(intersections_2)
    
    while(any(new_baseline_int[,2]>filter2b[,2])){
      x1 <- filter2b[1,1]
      # y-values
      y1 <- filter2b[1,3]
      # slope
      m1 <- (filter2b[nrow(filter2b),2] - filter2b[1,2])/(filter2b[nrow(filter2b),1]-  filter2b[1,1])
      # intersection with y-axis
      b1 <- y1 - m1 * x1
      
      # calculate the distance of each selected profile point to straight
      # profile point with largest distance equals the local minimum
      distances <- data.frame()
      for (j in 1:nrow(filter2b)) {
        # x-values
        x2 <- filter2b[j, 1]
        # y-values
        y2 <- filter2b[j, 2]
        # y-value of new baseline
        y3 <- new_baseline_int[j,2]
        # slope
        m2 <- (-1)/m1
        # intersection with y-axis
        b2 <- y2 - m2 * x2
        
        # intersection with connecting straight
        x_intersect <- (b2 - b1)/(m1 - m2)
        y_intersect <- m1 * x_intersect + b1
        
        dx <- x_intersect - x2
        dy <- y_intersect - y2
        
        dist <- abs(sqrt(dx^2 + dy^2))
        temp <- data.frame(x2, y2, y3, dist)
        distances <- rbind(distances, temp)
      }
      distances <- distances[which((distances$y3-distances$y2)>0),]
      if(nrow(distances)==0) break
      select_order <- distances[order(distances$dist), ]
      min_point <- select_order[nrow(select_order), 1:3]
      minimum <- min_point[, 1:2]
      colnames(minimum) <- colnames(intersections_2)
      bp_add <- rbind(bp_add, minimum)
      
      colnames(minimum) <- colnames(intersections_2)
      colnames(bp_add) <- colnames(baseline_extended)
      
      
      new_baseline <- rbind(baseline_extended[i,],baseline_extended[i+1,], bp_add)
      new_baseline <- new_baseline[order(new_baseline$X), ]
      new_baseline_int <- approx(new_baseline[1:nrow(new_baseline), ], xout = filter2[,1])
      new_baseline_int <- data.frame(X=new_baseline_int$x, Z=new_baseline_int$y)
      new_baseline_int <- new_baseline_int[new_baseline_int[,1]>filter2[1,1],]
    }
  }
  return(bp_add)
}



# add basepoints layer 2

add_basepoints_l2 <- function(intersections_2, profile_plus_baseline, baseline_extended, min_modus){
  
  bp_add <- data.frame()
  colnames(intersections_2) <- c("X","Z")
  baseline_extended <- baseline_extended[!duplicated(baseline_extended$X), ]
  baseline_extended$X <- round(baseline_extended$X,3)
  baseline_extended <- baseline_extended[order(baseline_extended$X),]
  intersections_2[,1] <- round(intersections_2[,1],3)  
  
  for(i in 1:(nrow(baseline_extended)-1)){
    
    intersections_i <- intersections_2[intersections_2$X>baseline_extended[i,1]&intersections_2$X<baseline_extended[i+1,1],]
    if(nrow(intersections_i)==0&&(i!=nrow(baseline_extended)-1)) next
    
    # x-values
    x1 <- baseline_extended[i,1]
    # y-values
    y1 <- baseline_extended[i,2]
    # slope
    m1 <- (baseline_extended[i+1,2] - baseline_extended[i,2])/(baseline_extended[i+1,1] - baseline_extended[i,1])
    # intersection with y-axis
    b1 <- y1 - m1 * x1
    
    filter2 <- profile_plus_baseline[which(profile_plus_baseline[,1]>baseline_extended[i,1]&profile_plus_baseline[,1]<baseline_extended[i+1,1]),]
    
    # calculate the distance of each selected profile point to straight
    # profile point with largest distance equals the local minimum
    distances <- data.frame()
    for (j in 1:nrow(filter2)) {
      # x-values
      x2 <- filter2[j, 1]
      # y-values of baseline1/2
      y2 <- filter2[j, 2]
      # y-value of profile/baseline1
      y3 <- filter2[j, 3]
      # slope
      m2 <- (-1)/m1
      # intersection with y-axis
      b2 <- y2 - m2 * x2
      
      # intersection with connecting straight
      x_intersect <- (b2 - b1)/(m1 - m2)
      y_intersect <- m1 * x_intersect + b1
      
      dx <- x_intersect - x2
      dy <- y_intersect - y2
      
      dist <- abs(sqrt(dx^2 + dy^2))
      temp <- data.frame(x2, y2, y3, dist)
      distances <- rbind(distances, temp)
    }
    distances <- distances[which((distances$y3-distances$y2)>0),]
    if(nrow(distances)==0) next
    select_order <- distances[order(distances$dist), ]
    min_point <- select_order[nrow(select_order), 1:3]
    minimum <- min_point[, 1:2]
    if(i==1&min_modus==T){
      minimum <- filter2[which.min(filter2[,2]),1:2]
      colnames(minimum) <- colnames(intersections_2)
      bp_add <- rbind(bp_add, minimum)
      next
    }
    colnames(minimum) <- colnames(intersections_2)
    bp_add <- rbind(bp_add, minimum)
    
    if(nrow(minimum)==0)next
    
    colnames(bp_add) <- colnames(baseline_extended)  
    #new_baseline <- baseline_extended[!duplicated(baseline_extended$X), ]  
    new_baseline <- rbind(baseline_extended[i,],baseline_extended[i+1,], bp_add)
    new_baseline <- new_baseline[order(new_baseline$X), ]
    new_baseline_int <- approx(new_baseline[1:nrow(new_baseline), ], xout = filter2[,1])
    new_baseline_int <- data.frame(X=new_baseline_int$x, Z=new_baseline_int$y)
    new_baseline_int[new_baseline_int[,1]>filter2[1,1],]
    base_point <- new_baseline[new_baseline$X>minimum$x, ][1,]
    filter2b <- filter2[filter2[,1]>=new_baseline_int[1,1], ]
    # recalculate intersections for newly added basepoint with next intersection point
    # if there are still intersections another basepoint has to be added
    colnames(minimum) <- colnames(intersections_2)
    
    while(any(round(new_baseline_int[,2],3)>filter2b[,2])){
      x1 <- new_baseline[1,1]
      # y-values
      y1 <- new_baseline[1,2]
      # slope
      m1 <- (new_baseline[1,2] - new_baseline[(nrow(new_baseline)-1),2])/(new_baseline[1,1]-  new_baseline[(nrow(new_baseline)-1),1])
      # intersection with y-axis
      b1 <- y1 - m1 * x1
      
      # calculate the distance of each selected profile point to straight
      # profile point with largest distance equals the local minimum
      distances <- data.frame()
      filter2b_sel <- filter2b[which((filter2b[,2]-new_baseline_int[,2])<0),]
      for (jj in 1:nrow(filter2b_sel)) {
        # x-values
        j <- as.numeric(rownames(filter2b_sel[jj,]))
        x2 <- filter2b[j, 1]
        # y-values
        y2 <- filter2b[j, 2]
        # y-value of new baseline
        y3 <- new_baseline_int[j,2]
        # slope
        m2 <- (-1)/m1
        # intersection with y-axis
        b2 <- y2 - m2 * x2
        
        # intersection with connecting straight
        x_intersect <- (b2 - b1)/(m1 - m2)
        y_intersect <- m1 * x_intersect + b1
        
        dx <- x_intersect - x2
        dy <- y_intersect - y2
        
        dist <- abs(sqrt(dx^2 + dy^2))
        temp <- data.frame(x2, y2, y3, dist)
        distances <- rbind(distances, temp)
      }
      distances <- distances[which((distances$y3-distances$y2)>0),]
      if(nrow(distances)==0) break
      select_order <- distances[order(distances$dist), ]
      min_point <- select_order[nrow(select_order), 1:3]
      minimum <- min_point[, 1:2]
      #minimum <- filter2b[which.min(filter2b[,2]),1:2]
      colnames(minimum) <- colnames(intersections_2)
      bp_add <- rbind(bp_add, minimum)
      
      colnames(minimum) <- colnames(intersections_2)
      colnames(bp_add) <- colnames(baseline_extended)
      
      
      new_baseline <- rbind(baseline_extended[i,],baseline_extended[i+1,], bp_add)
      new_baseline <- new_baseline[order(new_baseline$X), ]
      new_baseline_int <- approx(new_baseline, xout = filter2[,1])
      new_baseline_int <- data.frame(X=new_baseline_int$x, Z=new_baseline_int$y)
      #filter2b <- filter2[filter2[,1]>new_baseline[(nrow(new_baseline)-1),1], ]
    }
  }
  return(bp_add)
}






### zerocrossing for scale 1

zerocrossing_scale_1 <- function(input, iteration, BEP, count_scales, 
    result_path, min_modus, sample_frequency) {
    window_scale_1 <- input[1]
    window_scale_2 <- input[2]
    zc_threshold <- input[3]
    profile_original <- BEP
    profile <- BEP[ ,1:2]
    if (nrow(profile) == 0) {
        baseline_scale_1 <- NULL
    }
    if (is.na(window_scale_1)) {
        baseline_scale_1 <- NULL
    } else {
        colnames(profile) <- c("X", "Z")

        # define increment for moving average (equals windowlength of scale 1)
        increment <- as.integer(window_scale_1)

        # rename rows
        rownames(profile) <- seq(length = nrow(profile))
        
        # calculate moving average over the profile dividing profile into an
        # inner area as well as a left and right margin
        profile_inner <- profile [which(as.numeric(rownames(profile)) > 
            ceiling(0.5 * increment) & as.numeric(rownames(profile)) <= 
            nrow(profile) - 0.5 * increment), ]
        vector_inner <- rep(1/(increment + 1), increment + 1)
        filter_inner <- stats::filter(profile$Z, vector_inner, sides = 2)
        mov_avg_inner <- data.frame(profile_inner$X, filter_inner[(ceiling(0.5 * 
            increment + 1)):(nrow(profile) - 0.5 * increment)])
        # calculate moving average for left margin with constantly decreasing
        # increment
        mov_avg_left <- data.frame()
        if (increment <= 2) {
            mov_avg_left <- NULL
        } else {
            for (i in 2:(ceiling(0.5 * increment))) {
                window_left <- i - 1
                profile_left <- profile[(i - window_left):(i + window_left), 
                  ]
                mean <- mean.default(profile_left$Z)
                x <- profile[i, 1]
                df <- data.frame(x, mean)
                mov_avg_left <- rbind(mov_avg_left, df)
            }
            colnames(mov_avg_left) <- c("X", "Z_avg")
        }
        
        # calculate moving average for right margin with constantly decreasing
        # increment
        mov_avg_right <- data.frame()
        for (i in ((nrow(profile) - ceiling(0.5 * increment)) + 1):(nrow(profile))) {
            window_right <- nrow(profile) - i
            profile_right <- profile[(i - window_right):(i + window_right), 
                ]
            mean <- mean.default(profile_right$Z)
            x <- profile[i, 1]
            df <- data.frame(x, mean)
            mov_avg_right <- rbind(mov_avg_right, df)
        }
        
        # moving average for very first point equals profile
        mov_avg_1 <- profile[1, ]
        colnames(mov_avg_1) <- colnames(mov_avg_right) <- colnames(mov_avg_inner) <- c("X", 
            "Z_avg")
        
        # writing all moving averages in one vector
        mov_avg_all <- rbind(mov_avg_1, mov_avg_left, mov_avg_inner, mov_avg_right)
        
        # writing moving average in a dataframe with detrend profile
        profile_plus_avg <- data.frame(profile$X, profile$Z, mov_avg_all$Z_avg)
        colnames(profile_plus_avg) <- c("X", "Z", "Z_avg")
        profile_plus_avg_x <- profile_plus_avg[, 1]
        profile_plus_avg_z <- profile_plus_avg[, 2]
        profile_plus_avg_z_avg <- profile_plus_avg[, 3]
        profile_plus_avg <- profile_plus_avg[order(profile_plus_avg[, 1]), 
            ]

        # calculating straights bewtween profile points
        profile_straights_fun <- function(i) {
            # x-values
            x <- profile_plus_avg_x[i]
            # z-values
            y <- profile_plus_avg_z[i]
            # slope
            m <- (profile_plus_avg_z[i + 1] - profile_plus_avg_z[i])/(profile_plus_avg_x[i + 
                1] - profile_plus_avg_x[i])
            # intercection with z-axis
            b <- y - m * x
            
            profile_straight <- matrix(c(x, m, b), 1, 3)
            return(profile_straight)
        }
        
        vec <- c(1:nrow(profile_plus_avg))
        profile_straights <- sapply(vec, profile_straights_fun)
        profile_straights <- data.frame(t(profile_straights))
        colnames(profile_straights) <- c("X", "m", "b")
        
        # calculating straights bewtween moving average points
        movavg_straights_fun <- function(i) {
            # x-values
            x <- profile_plus_avg_x[i]
            # z-values
            y <- profile_plus_avg_z_avg[i]
            # slope
            m <- (profile_plus_avg_z_avg[i + 1] - profile_plus_avg_z_avg[i])/(profile_plus_avg_x[i + 
                1] - profile_plus_avg_x[i])
            # intercection with z-axis
            b <- y - m * x
            
            movavg_straight <- matrix(c(x, m, b), 1, 3)
            return(movavg_straight)
        }
        
        vec <- c(1:nrow(profile_plus_avg))
        movavg_straights <- sapply(vec, movavg_straights_fun)
        movavg_straights <- data.frame(t(movavg_straights))
        colnames(movavg_straights) <- c("x", "m", "b")
        
        # calculate intersections bewtween both straights dataframes (for
        # determining the zero crossings)
        intersections <- data.frame(x = rep(NA, nrow(movavg_straights)), 
                                           y = rep(NA, nrow(movavg_straights)))
        condition_2 <- (profile_straights[, 3] - movavg_straights[, 3])/(movavg_straights[, 
                                                                                          2] - profile_straights[, 2])
        condition_1a <- movavg_straights[, 3] == 0
        condition_1b <- is.na(profile_straights[, 2])
        condition_1c <- is.na(movavg_straights[, 2])
        condition_1d <- movavg_straights[, 2] - profile_straights[, 2] == 
          0
        
        # for-loop over all straights
        for (i in 1:nrow(movavg_straights)) {
          if (condition_1b[i] || condition_1c[i]) {
            next
          } else {
            if(condition_1d[i]){
              check <- profile_plus_avg[profile_plus_avg$X==movavg_straights[i,1],]

              if(round(check[,2],3)==round(check[,3],3)){
                if(i==1){
                  intersections[i, 1] <- check[,1]
                  intersections[i, 2] <- check[,2]
                }else{
                  if(is.na(intersections[i-1,1])){
                    intersections[i, 1] <- check[,1]
                    intersections[i, 2] <- check[,2]
                  }
                }
                next
              }else{
                next
              }
              
            }
            
            
            # At the point of intersection y1=y2 and x1=x2 -> m1x+b1 = m2x+b2 ->
            # converting to x check if there is an intersection in the current
            # interval
            if (round(condition_2[i],3) >= round(profile_plus_avg[i, 
                                                                         1], digits = 3) && round(condition_2[i],3) <= round(profile_plus_avg[(i + 
                                                                                                                                                      1), 1], digits = 3)) {
              x <- (profile_straights[i, 3] - movavg_straights[i, 3])/(movavg_straights[i, 
                                                                                        2] - profile_straights[i, 2])
              y <- profile_straights[i, 2] * x + profile_straights[i, 
                                                                   3]
              intersections[i, 1] <- x
              intersections[i, 2] <- y
            }
          }
        }
        
        # filter out NA-values
        intersections$x <- round(intersections$x,3)
        intersections <- filter(intersections, !is.na(intersections$x))
        intersections <- intersections[!duplicated(intersections$x),]

        # create dataframe for minima/maxima
        minmax <- data.frame(X = rep(NA, (nrow(intersections) - 1)), Z = rep(NA, 
            (nrow(intersections) - 1)), label = rep(NA, (nrow(intersections) - 
            1)))
        
        # skip if less than 3 intersections 
        if(nrow(intersections)<3){
          print("error - window size 1 or zc_threshold too large")
          return("skip")
        }
        
        # for-loop over all intersections selecting all profile points between
        # two intersections and order by z-value first point corresponds to a
        # potential minimum, last point corresponds to a potential maximum if
        # point is above the moving average it is a maximum otherwise it is a
        # minimum adding label (1=minimum, 2=maximum)

        for (i in 1:(nrow(intersections) - 1)) {
          select <- profile_plus_avg[which(profile_plus_avg$X >= 
                                                    intersections[i, 1] & profile_plus_avg$X <= 
                                                    intersections[(i + 1), 1]), ]
          if(nrow(select)==0) next
          
      
          if (!is.na(select[1, 1])&&(nrow(select)>2)) {
            
            if(min_modus==T){
              if(all(select[2:(nrow(select)-1),3]>=select[2:(nrow(select)-1),2])){
                extremum <- select[which.min(select[,2]), ] 
              } else {
                extremum <- select[which.max(select[,2]), ] 
              }
            } else {
              # alternative way to detrmine the minimum by detrending the currently selected section
              # leads to same results as procedure in zerocrossing.V2
              # detrend_select <- data.frame(z=select$Z-select$Z_avg)
              # detrend_select$id <- as.numeric(rownames(detrend_select))
              # detrend_select <- detrend_select[order(detrend_select$z),]
              # minimum <- detrend_select[1,2]
              
              # alternative way to detrmine the minimum by creating a connection line between first and last point of moving average curve
              # of currently selected section (select[ ,3])
              # Calculating distance between each currently selected profile point and connection line and order by distance
              # This procedure is more suitable than V2 because local fluctuations of the moving average curve are eliminated by
              # creating the connection line. This way the actual dune troughs are represented
              #straight_connection <- coords2Lines(cbind(c(select[1,1],select[nrow(select),1]),c(select[1,3],select[nrow(select),3])),ID=1)
              m <- (select[nrow(select),3]-select[1,3])/(select[nrow(select),1]-select[1,1])
              b <- select[1,3]-m*select[1,1]
              m_inverse <- -1/m
              #print(straight_connection)
              #print(select)
              
              for(j in 1:nrow(select)){
                
                if(m!=0){
                  b_j <- select[j,2]-m_inverse*select[j,1]
                  y1 <- -1000
                  y2 <- 1000
                  x1 <- (y1-b_j)/m_inverse
                  x2 <- (y2-b_j)/m_inverse
                  x_cross <- (b_j-b)/(m-m_inverse)
                  y_cross <- m*x_cross+b
                  #int_2 <- data.frame(x_cross,y_cross)
                  #cross_line <- coords2Lines(cbind(c(x1,x2),c(y1,y2)),ID=2)
                  #intersection <- as.data.frame(coordinates(raster::intersect(straight_connection,cross_line)))
                  
                  #if(nrow(intersection)>0){
                  distance <- sqrt((select[j,2]-y_cross)^2 + (select[j,1]-x_cross)^2)
                  #} else {
                  #  distance <- 0
                  # }
                } else {
                  distance <- abs(select[j,2]-select[j,3])
                }
                select$distance[j] <- distance
              }
              
              order <- select[order(select$distance, decreasing=T),]
              extremum <- order[1, ]
            }
            
          } else if (nrow(select)<=2){
            if(all(select[,3]>=select[,2])){
              extremum <- select[which.min(select[,2]), ] 
            } else {
              extremum <- select[which.max(select[,2]), ] 
            }
          }
            if (extremum[1, 2] > extremum[1, 3]) {
              {
                point <- extremum[1, 1:2]
                label <- 2  #'max_pot'
                minmax[i, 1] <- point$X
                minmax[i, 2] <- point$Z
                minmax[i, 3] <- label
              }
            } else {
              point <- extremum[1, 1:2]
              label <- 1  #'min_pot'
              minmax[i, 1] <- point$X
              minmax[i, 2] <- point$Z
              minmax[i, 3] <- label
            }
        } 
        
        minmax <- filter(minmax, !is.na(minmax$X))
        
        # converting minmax to spatial points
        sp <- minmax[, 1:3]
        sp_points <- SpatialPointsDataFrame(sp, data = data.frame(label = sp$label))
        
        # converting moving average to spatial line
        MEAN_Lines <- data.frame(profile_plus_avg$X, profile_plus_avg$Z_avg)
        colnames(MEAN_Lines) <- c("X", "Z")
        MEAN_Lines <- filter(MEAN_Lines, abs(profile_plus_avg$Z_avg) > 
            0)
        
        L1 <- Line(MEAN_Lines)
        S1 <- Lines(list(L1), ID = "a")
        Sl <- SpatialLines(list(S1))
        
        # define Buffer: extension corresponds to zc_threshold
        Buffer <- gBuffer(Sl, width = zc_threshold, capStyle = "FLAT")
        
        # select points inside buffer in order to delete them
        sp_points_sel <- sp_points[Buffer, ]
        points_sel <- data.frame(coordinates(sp_points_sel))
        
        # create new dataframe that contains only min/max values outside buffer
        # area
        minmax_filtered <- anti_join(sp, points_sel, by = c("X", "Z", "label"))
        
        # create baseline: contains only minima if the maximum between two
        # minima was deleted in the previous step minimum is deleted, too.
        
        if (is.null(minmax_filtered[1, 1]) | is.na(minmax_filtered[1, 1])) {
            baseline_scale_1 <- NULL
            print("error - window size 1 or zc_threshold too small")
            return("skip")
        } else {
            baseline <- minmax_filtered[which(minmax_filtered$label == 
                1), 1:2]
            if (is.null(baseline[1, 1])) {
                baseline_scale_1 <- NULL
                print("error - window size 1 or zc_threshold too small")
                return("skip")
            } else {
                colnames(baseline) <- c("X", "Z")
                
                # extend baseline by adding the boundary points of the profile
                baseline_extended <- rbind(profile[1, 1:2], baseline, profile[nrow(profile), 
                  1:2])
                
                # interpolate
                interpolated <- approx(baseline_extended$X, baseline_extended$Z, 
                  xout = profile_plus_avg$X)
                interpolated_df <- data.frame(interpolated$x, interpolated$y)
                
                # create dataframe containing x-values, profile points, and base points
                profile_plus_baseline <- data.frame(profile_plus_avg$X, 
                  profile_plus_avg$Z, interpolated_df$interpolated.y)
                
                # the base line must not cross the profile in any point therefore all
                # intersections are calculated then, all profile points bewteen two
                # intersection are selected and the minimum is calaculated the minimum
                # is added to the base line
                
                # calcuate straights between all base points
                bp_straights <- data.frame(x = rep(NA, nrow(baseline_extended)), 
                  m = rep(NA, nrow(baseline_extended)), b = rep(NA, nrow(baseline_extended)))
                BPE_x <- baseline_extended[, 1]
                BPE_y <- baseline_extended[, 2]
                
                # for-loop over all profile points and calculate straight between each
                # two points
                for (i in 1:nrow(baseline_extended)) {
                  
                  # x-values
                  x <- BPE_x[i]
                  # y-values
                  y <- BPE_y[i]
                  # slope
                  m <- (BPE_y[i + 1] - BPE_y[i])/(BPE_x[i + 1] - BPE_x[i])
                  # intersection with y-axis
                  b <- y - m * x
                  
                  bp_straights[i, 1] <- x
                  bp_straights[i, 2] <- m
                  bp_straights[i, 3] <- b
                }
                
                # calculate intersections between profile_straights and bp_straights
                profile_straights <- profile_straights[which(!is.na(profile_straights$m)), 
                  ]
                profile_straights_x <- profile_straights[, 1]
                bp_straights <- bp_straights[which(!is.na(bp_straights$m)), 
                  ]
                bp_straights_x <- bp_straights[, 1]
                bp_straights_m <- bp_straights[, 2]
                bp_straights_b <- bp_straights[, 3]
                intersections_2 <- data.frame()

                # for-loop over all straights
                for (i in 1:(nrow(bp_straights) - 1)) {
                  # filter all profile_straights between two bp_straights
                  filter1 <- profile_straights[which(profile_straights_x >= 
                    bp_straights_x[i] & profile_straights_x <= bp_straights_x[i + 
                    1]), ]
                  if(nrow(filter1)==0) next
                  for (j in 1:(nrow(filter1) - 1)) {
                    # check for NA-values
                    if (bp_straights_m[i] - filter1[j, 2] != 0) {
                      
                      # At the point of intersection y1=y2 and x1=x2 -> m1x+b1 = m2x+b2 ->
                      # converting to x
                      condition_1 <- (filter1[j, 3] - bp_straights_b[i])/(bp_straights_m[i] - 
                        filter1[j, 2])
                      # check if there is an intersection in the current interval if TRUE:
                      # condition_1==x_value and mx+b==y-value
                      if (condition_1 >= round(filter1[j, 1], digits = 3) && 
                        condition_1 <= round(filter1[(j + 1), 1], digits = 3)) {
                        x <- (filter1[j, 3] - bp_straights_b[i])/(bp_straights_m[i] - 
                          filter1[j, 2])
                        y <- filter1[j, 2] * x + filter1[j, 3]
                        dtemp3 <- data.frame(x, y)
                        intersections_2 <- rbind(intersections_2, dtemp3)
                      }
                    }
                  }
                }
                
                # check if there are any intersections at all
                if (nrow(intersections_2) == 0) {
                    bp_add <- NULL
                    baseline_new <- baseline_extended
                    baseline_new_order <- baseline_new[order(baseline_new$X),]
                } else {
                
                # removing double points (distance to base points must be larger than >
                # 0.001m)
                intersections_2_new <- data.frame(x = rep(NA, nrow(intersections_2) - 
                  1), y = rep(NA, nrow(intersections_2) - 1), diff_x = rep(NA, 
                  nrow(intersections_2) - 1), diff_y = rep(NA, nrow(intersections_2) - 
                  1))
                for (i in 1:(nrow(intersections_2) - 1)) {
                  diff_x <- abs(intersections_2[i + 1, 1] - intersections_2[i, 
                    1])
                  diff_y <- abs(intersections_2[i + 1, 2] - intersections_2[i, 
                    2])
                  intersections_2_new[i, 1] <- intersections_2[i, 1]
                  intersections_2_new[i, 2] <- intersections_2[i, 2]
                  intersections_2_new[i, 3] <- diff_x
                  intersections_2_new[i, 4] <- diff_y
                }
                
                intersections_2_new_filtered <- filter(intersections_2_new, 
                  diff_x > 0.001 | diff_y > 0.001)
                intersections_2_new_filtered <- intersections_2_new_filtered[, 
                  1:2]
                colnames(intersections_2_new_filtered) <- colnames(intersections_2)
                intersections_2 <- rbind(intersections_2_new_filtered, 
                  intersections_2[nrow(intersections_2), ])
                
                bp_add <- add_basepoints_l1(intersections_2, profile_plus_baseline, baseline_extended)
            
                if (nrow(bp_add) > 0) {
                  colnames(bp_add) <- c("X", "Z")
                
                # filter out double points
                bp_try_select <- data.frame()
                for (i in 1:(nrow(baseline_extended) - 1)) {
                  bp_add_try <- bp_add[which(bp_add[, 1] > baseline_extended[i, 
                    1] & bp_add[, 1] < baseline_extended[i + 1, 1]), ]
                  bp_try_order <- bp_add_try[order(bp_add_try[, 2]), ]
                  bp_try_select <- rbind(bp_try_select, bp_try_order)
                }
                
                # filter out NA-values
                bp_try_select <- bp_try_select[which(!is.na(bp_try_select[, 
                  1])), ]
                bp_add <- bp_try_select
                
                # adding new base points to base line
                baseline_new <- rbind(baseline_extended, bp_add)
                baseline_new_order <- baseline_new[order(baseline_new$X), 
                  ]
                } else {
                    baseline_new <- baseline_extended
                    baseline_new_order <- baseline_new[order(baseline_new$X),]
                }
                }
                # interpolate new base points over profile
                baseline_new_int <- approx(baseline_new_order$X, baseline_new_order$Z, 
                  xout = profile_plus_baseline$profile_plus_avg.X)
                Basislinie_new_int_df <- data.frame(baseline_new_int$x, 
                  baseline_new_int$y)
                
                # in some points the base line might still cross the profile in these
                # points the baseline is projected on the profile
                difference_vec <- profile_plus_baseline[, 2] - Basislinie_new_int_df[, 
                  2]
                base_x_vec <- profile_plus_baseline[, 1]
                base_y_vec <- profile_plus_baseline[, 2]
                baseline_new_2 <- data.frame(X = base_x_vec, Z = base_y_vec, 
                  Diff = difference_vec)
                baseline_new_2 <- baseline_new_2[which(baseline_new_2$Diff < 
                  0), 1:2]
                
                # the projected points are added to the baseline
                baseline_scale_1 <- rbind(baseline_new_2, baseline_extended, 
                  bp_add)
                baseline_scale_1 <- baseline_scale_1[order(baseline_scale_1$X), 
                  ]

      plot(BEP[ ,1:2], type="l", xlab="x [m]", ylab="z [m]")#, ylim=c(0.8,1.4))
      lines(baseline_scale_1, col="red")
      legend("topright", legend=c("baseline layer 1"), col="red", lty=1, cex=0.7)
             
      #call zerocrossing function for scale 2
        zerocrossing_scale_2(profile, baseline_scale_1, window_scale_1, 
            window_scale_2, zc_threshold, iteration, 
            count_scales, result_path, sample_frequency, profile_original, min_modus)
            }
        }
    }
}



### zerocrossing for scale 2

zerocrossing_scale_2 <- function(BEP, baseline_scale_1, window_scale_1, 
    window_scale_2, zc_threshold, iteration, 
     count_scales, result_path, sample_frequency, profile_original, min_modus) {

    if (is.null(window_scale_2) || is.na(window_scale_2)) {
        baseline_scale_2 <- baseline_scale_1
    } else {
        if (is.null(baseline_scale_1) || all(is.na(baseline_scale_1[ ,2])) || nrow(baseline_scale_1)==0) {
            baseline_scale_1 <- BEP
        }
        colnames(BEP) <- c("X", "Z")
        
        # define increment for moving average (equals windowlength of scale 2)
        increment <- as.integer(window_scale_2)
        # interpolate base line scale 1 over profile points
        bl1_interpolated <- approx(baseline_scale_1[, 1], baseline_scale_1[, 
            2], xout = BEP[, 1])
        profile <- data.frame(bl1_interpolated$x, bl1_interpolated$y)
        #profile <- BEP
        colnames(profile) <- c("X", "Z")
        
        # calculate moving average over base line 1
        profile_inner <- filter(profile, as.numeric(rownames(profile)) > 
            ceiling(0.5 * increment) & as.numeric(rownames(profile)) <= 
            nrow(profile) - 0.5 * increment)
        profile_left <- filter(profile[1:increment, 1:2])
        
        # moving average for inner area
        vector_inner <- rep(1/(increment + 1), increment + 1)
        filter_inner <- stats::filter(profile$Z, vector_inner, sides = 2)
        mov_avg_inner <- data.frame(profile_inner$X, filter_inner[ceiling((0.5 * 
            increment + 1)):(nrow(profile) - 0.5 * increment)])
        
        # calculate moving average for left margin with constantly decreasing
        # increment
        mov_avg_left <- data.frame()
        
        for (i in 2:(ceiling(0.5 * increment))) {
            window_left <- i - 1
            profile_l <- profile[(i - window_left):(i + window_left), ]
            mean <- mean.default(profile_l$Z)
            x <- profile[i, 1]
            df <- data.frame(x, mean)
            mov_avg_left <- rbind(mov_avg_left, df)
        }
        
        # moving average for very first point equals profile
        mov_avg_1 <- mov_avg_left[1, ]
        
        # calculate moving average for right margin with constantly decreasing
        # increment
        mov_avg_right <- data.frame()
        
        for (i in ((nrow(profile) - ceiling(0.5 * increment) + 1)):(nrow(profile))) {
            # i=nrow(profile)
            window_right <- nrow(profile) - i
            profile_r <- profile[(i - window_right):(i + window_right), 
                ]
            mean <- mean.default(profile_r$Z)
            x <- profile[i, 1]
            df <- data.frame(x, mean)
            mov_avg_right <- rbind(mov_avg_right, df)
        }
        
        # writing all moving averages in one vector
        colnames(mov_avg_1) <- colnames(mov_avg_left) <- colnames(mov_avg_right) <- colnames(mov_avg_inner) <- c("X", 
            "Z_avg")
        mov_avg_all <- rbind(mov_avg_1, mov_avg_left, mov_avg_inner, mov_avg_right)

        # writing moving average in a dataframe with detrend profile
        profile_plus_avg_scale2 <- data.frame(profile$X, profile$Z, mov_avg_all$Z_avg)
        colnames(profile_plus_avg_scale2) <- c("X", "Z", "Z_avg")
        profile_plus_avg_scale2_x <- profile_plus_avg_scale2[, 1]
        profile_plus_avg_scale2_z <- profile_plus_avg_scale2[, 2]
        profile_plus_avg_scale2_z_avg <- profile_plus_avg_scale2[, 3]
        profile_plus_avg_scale2 <- profile_plus_avg_scale2[order(profile_plus_avg_scale2[, 
            1]), ]
        
        # calculating straights bewtween profile points
        profile_straights_fun <- function(i) {
            # x-values
            x <- profile_plus_avg_scale2_x[i]
            # z-values
            y <- profile_plus_avg_scale2_z[i]
            # slope
            m <- (profile_plus_avg_scale2_z[i + 1] - profile_plus_avg_scale2_z[i])/(profile_plus_avg_scale2_x[i + 
                1] - profile_plus_avg_scale2_x[i])
            # intercection with z-axis
            b <- y - m * x
            
            profile_straight <- matrix(c(x, m, b), 1, 3)
            return(profile_straight)
        }
        
        vec <- c(1:nrow(profile_plus_avg_scale2))
        profile_straights <- sapply(vec, profile_straights_fun)
        profile_straights <- data.frame(t(profile_straights))
        colnames(profile_straights) <- c("X", "m", "b")
        
        # calculating straights bewtween moving average points
        movavg_straights_fun <- function(i) {
            # x-values
            x <- profile_plus_avg_scale2_x[i]
            # z-values
            y <- profile_plus_avg_scale2_z_avg[i]
            # slope
            m <- (profile_plus_avg_scale2_z_avg[i + 1] - profile_plus_avg_scale2_z_avg[i])/(profile_plus_avg_scale2_x[i + 
                1] - profile_plus_avg_scale2_x[i])
            # intercection with z-axis
            b <- y - m * x
            
            movavg_straight <- matrix(c(x, m, b), 1, 3)
            return(movavg_straight)
        }
        
        vec <- c(1:nrow(profile_plus_avg_scale2))
        movavg_straights <- sapply(vec, movavg_straights_fun)
        movavg_straights <- data.frame(t(movavg_straights))
        colnames(movavg_straights) <- c("x", "m", "b")
        
        # calculate intersections bewtween both straights dataframes (for
        # determining the zero crossings)
        intersections_scale2 <- data.frame(x = rep(NA, nrow(movavg_straights)), 
            y = rep(NA, nrow(movavg_straights)))
        condition_2 <- (profile_straights[, 3] - movavg_straights[, 3])/(movavg_straights[, 
            2] - profile_straights[, 2])
        condition_1a <- movavg_straights[, 3] == 0
        condition_1b <- is.na(profile_straights[, 2])
        condition_1c <- is.na(movavg_straights[, 2])
        condition_1d <- movavg_straights[, 2] - profile_straights[, 2] == 
            0
        
        # for-loop over all straights
        for (i in 1:nrow(movavg_straights)) {
            if (condition_1b[i] || condition_1c[i]) {
                next
            } else {
              if(condition_1d[i]){
                check <- profile_plus_avg_scale2[profile_plus_avg_scale2$X==movavg_straights[i,1],]
                if(round(check[,2],3)==round(check[,3],3)){
                  if(i==1){
                    intersections_scale2[i, 1] <- check[,1]
                    intersections_scale2[i, 2] <- check[,2]
                  }else{
                    if(is.na(intersections_scale2[i-1,1])){
                      intersections_scale2[i, 1] <- check[,1]
                      intersections_scale2[i, 2] <- check[,2]
                    }
                  }
                  next
                }else{
                  next
                }
                
              }
              
              
                # At the point of intersection y1=y2 and x1=x2 -> m1x+b1 = m2x+b2 ->
                # converting to x check if there is an intersection in the current
                # interval
                if (round(condition_2[i],3) >= round(profile_plus_avg_scale2[i, 
                  1], digits = 3) && round(condition_2[i],3) <= round(profile_plus_avg_scale2[(i + 
                  1), 1], digits = 3)) {
                  x <- (profile_straights[i, 3] - movavg_straights[i, 3])/(movavg_straights[i, 
                    2] - profile_straights[i, 2])
                  y <- profile_straights[i, 2] * x + profile_straights[i, 
                    3]
                  intersections_scale2[i, 1] <- x
                  intersections_scale2[i, 2] <- y
                }
            }
        }
        
        # filter out NA-values
        intersections_scale2$x <- round(intersections_scale2$x,3)
        intersections_scale2 <- filter(intersections_scale2, !is.na(intersections_scale2$x))
        intersections_scale2 <- intersections_scale2[!duplicated(intersections_scale2$x),]
        
        if(nrow(intersections_scale2) < 2){
          baseline_scale_2 <- baseline_scale_1
        } else {
        
        # create dataframe for minima/maxima
        minmax <- data.frame(X = rep(NA, 2*((nrow(intersections_scale2) - 
            1))), Z = rep(NA, 2*((nrow(intersections_scale2) - 1))), label = rep(NA, 
            2*((nrow(intersections_scale2) - 1))))
        
        # for-loop over all intersections selecting all profile points between
        # two intersections and order by z-value first point corresponds to a
        # potential minimum, last point corresponds to a potential maximum if
        # point is above the moving average it is a maximum otherwise it is a
        # minimum adding label (1=minimum, 2=maximum)

        for (i in 1:(nrow(intersections_scale2) - 1)) {
            select <- profile_plus_avg_scale2[which(profile_plus_avg_scale2$X > 
                intersections_scale2[i, 1] & profile_plus_avg_scale2$X < 
                intersections_scale2[(i + 1), 1]), ]

            if(nrow(select)==0) next
            
            if (!is.na(select[1, 1])&&(nrow(select)>2)) {
              if(min_modus==T){
                if(all(select[,3]>=select[,2])){
                  extremum <- select[which.min(select[,2]), ] 
                } else {
                  extremum <- select[which.max(select[,2]), ] 
                }
              }else{
                # alternative way to detrmine the minimum by creating a connection line between first and last point of moving average curve
                # of currently selected section (select[ ,3])
                # Calculating distance between each currently selected profile point and connection line and order by distance
                # This procedure is more suitable than V2 because local fluctuations of the moving average curve are eliminated by
                # creating the connection line. This way the actual dune troughs are represented
                #straight_connection <- coords2Lines(cbind(c(select[1,1],select[nrow(select),1]),c(select[1,3],select[nrow(select),3])),ID=1)
                
                m <- (select[nrow(select),3]-select[1,3])/(select[nrow(select),1]-select[1,1])
                b <- select[1,3]-m*select[1,1]
                m_inverse <- -1/m
                
                for(j in 1:nrow(select)){
                  if(m!=0){
                    b_j <- select[j,2]-m_inverse*select[j,1]
                    y1 <- -1000
                    y2 <- 1000
                    x1 <- (y1-b_j)/m_inverse
                    x2 <- (y2-b_j)/m_inverse
                    x_cross <- (b_j-b)/(m-m_inverse)
                    y_cross <- m*x_cross+b
                    
                    distance <- sqrt((select[j,2]-y_cross)^2 + (select[j,1]-x_cross)^2)
                    #} else {
                    #  distance <- 0
                    # }
                  } else {
                    distance <- abs(select[j,2]-select[j,3])
                  }
                  select$distance[j] <- distance
                } 
                order <- select[order(select$distance, decreasing=T),]
                extremum <- order[1, ]
              }
              
                #plot(profile_plus_avg[,c(1,2)],type="l")#, xlim=c(80,90),ylim=c(76.30,76.35))
                #lines(profile_plus_avg_scale2[,c(1,2)],col="red", xlim=c(730,750))
                #lines(profile_plus_avg_scale2[,c(1,3)],col="green", xlim=c(730,750))
                #lines(select, col="blue", type="b")
                #lines(avg_selected, col="green",type="b")
                #points(select[minimum, ],pch=2)
                
                # alternative way to detrmine the minimum by detrending the currently selected section
                # leads to same results as procedure in zerocrossing.V2
                # detrend_select <- data.frame(z=select$Z-select$Z_avg)
                # detrend_select$id <- as.numeric(rownames(detrend_select))
                # detrend_select <- detrend_select[order(detrend_select$z),]
                # minimum <- detrend_select[1,2]
    
            } else if (nrow(select)<=2){
              if(all(select[,3]>=select[,2])){
                extremum <- select[which.min(select[,2]), ] 
              } else {
                extremum <- select[which.max(select[,2]), ] 
              }
            }
                if (extremum[1, 2] > extremum[1, 3]) {
                  {
                    point <- extremum[1, 1:2]
                    label <- 2  #'max_pot'
                    minmax[(i-1)*2+2, 1] <- point$X
                    minmax[(i-1)*2+2, 2] <- point$Z
                    minmax[(i-1)*2+2, 3] <- label
                  }
                } else {
                  point <- extremum[1, 1:2]
                  label <- 1  #'min_pot'
                  minmax[(i-1)*2+2, 1] <- point$X
                  minmax[(i-1)*2+2, 2] <- point$Z
                  minmax[(i-1)*2+2, 3] <- label
                }
            }
        minmax <- filter(minmax, !is.na(minmax$X))

        # converting minmax to spatial points
        sp <- minmax[, 1:3]
        sp_points <- SpatialPointsDataFrame(sp, data = data.frame(label = sp$label))

        # converting moving average to spatial line
        MEAN_Lines <- data.frame(profile_plus_avg_scale2$X, profile_plus_avg_scale2$Z_avg)
        colnames(MEAN_Lines) <- c("X", "Z")
        MEAN_Lines <- filter(MEAN_Lines, abs(profile_plus_avg_scale2$Z_avg) > 
            0)
        
        L1 <- Line(MEAN_Lines)
        S1 <- Lines(list(L1), ID = "a")
        Sl <- SpatialLines(list(S1))
        
        # define Buffer: extension corresponds to zc_threshold
        Buffer <- gBuffer(Sl, width = 0.0001, capStyle = "FLAT")
        
        # select points inside buffer in order to delete them
        sp_points_sel <- sp_points[Buffer, ]
        points_sel <- data.frame(coordinates(sp_points_sel))
        
        # create new dataframe that contains only min/max values outside buffer
        # area
        minmax_filtered_scale2 <- anti_join(sp, points_sel, by = c("X", 
            "Z", "label"))
        # create baseline: contains only minima if the maximum between two
        # minima was deleted in the previous step minimum is deleted, too.
        
        if (window_scale_1 == window_scale_2 | is.null(minmax_filtered_scale2[1, 
            1]) | is.na(minmax_filtered_scale2[1, 1])| nrow(minmax_filtered_scale2[minmax_filtered_scale2$label==1,])==0) {
            baseline_scale_2 <- baseline_scale_1
            print("no minima detected for layer 2")
        } else {
            baseline <- minmax_filtered_scale2[which(minmax_filtered_scale2$label == 
                1), 1:2]
            if (is.null(baseline[1, 1])) {
                baseline_scale_2 <- baseline_scale_1
            } else {
                colnames(baseline) <- c("X", "Z")
                
                # extend baseline by adding the boundary points of the profile
                baseline_extended <- rbind(profile[1, 1:2], baseline, profile[nrow(profile), 
                  1:2])
                
                # interpolate
                interpolated <- approx(baseline_extended$X, baseline_extended$Z, 
                  xout = profile_plus_avg_scale2$X)
                interpolated_df <- data.frame(interpolated$x, interpolated$y)
                
                # create dataframe containing x-values, profile points, and base points
                profile_plus_baseline <- data.frame(profile_plus_avg_scale2$X, 
                  profile_plus_avg_scale2$Z, interpolated_df$interpolated.y)
                
                # the base line must not cross the profile in any point therefore all
                # intersections are calculated then, all profile points bewteen two
                # intersection are selected and the minimum is calaculated the minimum
                # is added to the base line
                
                # calcuate straights between all base points
                bp_straights <- data.frame(x = rep(NA, nrow(baseline_extended)), 
                  m = rep(NA, nrow(baseline_extended)), b = rep(NA, nrow(baseline_extended)))
                BPE_x <- baseline_extended[, 1]
                BPE_y <- baseline_extended[, 2]
                
                # for-loop over all profile points and calculate straight between each
                # two points
                for (i in 1:nrow(baseline_extended)) {
                  
                  # x-values
                  x <- BPE_x[i]
                  # y-values
                  y <- BPE_y[i]
                  # slope
                  m <- (BPE_y[i + 1] - BPE_y[i])/(BPE_x[i + 1] - BPE_x[i])
                  # intersection with y-axis
                  b <- y - m * x
                  
                  bp_straights[i, 1] <- x
                  bp_straights[i, 2] <- m
                  bp_straights[i, 3] <- b
                }
                
                # calculate intersections between profile_straights and bp_straights
                profile_straights <- profile_straights[which(!is.na(profile_straights$m)), 
                  ]
                profile_straights_x <- profile_straights[, 1]
                bp_straights <- bp_straights[which(!is.na(bp_straights$m)), 
                  ]
                bp_straights_x <- bp_straights[, 1]
                bp_straights_m <- bp_straights[, 2]
                bp_straights_b <- bp_straights[, 3]
                intersections_2 <- data.frame()
                
                # for-loop over all straights
                for (i in 1:(nrow(bp_straights))) {
                   # filter all profile_straights between two bp_straights
                  if(nrow(bp_straights)==1){
                    filter1 <- profile_straights[which(profile_straights_x >= 
                                                         bp_straights_x[i] & profile_straights_x <= baseline_extended[nrow(baseline_extended),1]), ]
                  }else{
                    if(i==nrow(bp_straights)){
                      filter1 <- profile_straights[which(profile_straights_x >= 
                                                           bp_straights_x[i]), ]
                    }else{
                      filter1 <- profile_straights[which(profile_straights_x >= 
                                                           bp_straights_x[i] & profile_straights_x <= bp_straights_x[i + 
                                                                                                                       1]), ]  
                    }
                  }
                  if(nrow(filter1)<=1){
                      next
                  }
                  for (j in 1:(nrow(filter1) - 1)) {
                    # check for NA-values
                    if (bp_straights_m[i] - filter1[j, 2] != 0) {
                      
                      # At the point of intersection y1=y2 and x1=x2 -> m1x+b1 = m2x+b2 ->
                      # converting to x
                      condition_1 <- (filter1[j, 3] - bp_straights_b[i])/(bp_straights_m[i] - 
                        filter1[j, 2])
                      # check if there is an intersection in the current interval if TRUE:
                      # condition_1==x_value and mx+b==y-value
                      if (condition_1 >= round(filter1[j, 1], digits = 3) && 
                        condition_1 <= round(filter1[(j + 1), 1], digits = 3)) {
                        x <- (filter1[j, 3] - bp_straights_b[i])/(bp_straights_m[i] - 
                          filter1[j, 2])
                        y <- filter1[j, 2] * x + filter1[j, 3]
                        dtemp3 <- data.frame(x, y)
                        intersections_2 <- rbind(intersections_2, dtemp3)
                      }
                    }
                  }
                }

                # check if there are any intersections at all
                if (nrow(intersections_2) == 0) {
                    bp_add <- data.frame()
                } else {
                
                # removing double points (distance to base points must be larger than >
                # 0.001m)
                intersections_2_new <- data.frame(x = rep(NA, nrow(intersections_2) - 
                  1), y = rep(NA, nrow(intersections_2) - 1), diff_x = rep(NA, 
                  nrow(intersections_2) - 1), diff_y = rep(NA, nrow(intersections_2) - 
                  1))
                for (i in 1:(nrow(intersections_2) - 1)) {
                  diff_x <- abs(intersections_2[i + 1, 1] - intersections_2[i, 
                    1])
                  diff_y <- abs(intersections_2[i + 1, 2] - intersections_2[i, 
                    2])
                  intersections_2_new[i, 1] <- intersections_2[i, 1]
                  intersections_2_new[i, 2] <- intersections_2[i, 2]
                  intersections_2_new[i, 3] <- diff_x
                  intersections_2_new[i, 4] <- diff_y
                }
                intersections_2_new_filtered <- filter(intersections_2_new, 
                  diff_x > 0.001 | diff_y > 0.001)
                intersections_2_new_filtered <- intersections_2_new_filtered[, 
                  1:2]
                colnames(intersections_2_new_filtered) <- colnames(intersections_2)
                intersections_2 <- rbind(intersections_2_new_filtered, 
                  intersections_2[nrow(intersections_2), ])
    
                bp_add <- add_basepoints_l2(intersections_2, profile_plus_baseline, baseline_extended, min_modus=F)

                if (nrow(bp_add) > 0) {
                  colnames(bp_add) <- c("X", "Z")
                  
                  # filter out double points
                  bp_try_select <- data.frame()
                  for (i in 1:(nrow(baseline_extended) - 1)) {
                    bp_add_try <- bp_add[which(bp_add[, 1] > baseline_extended[i, 
                                                                               1] & bp_add[, 1] < baseline_extended[i + 1, 1]), ]
                    bp_try_order <- bp_add_try[order(bp_add_try[, 2]), ]
                    bp_try_select <- rbind(bp_try_select, bp_try_order)
                  }
                  
                  # filter out NA-values
                  bp_try_select <- bp_try_select[which(!is.na(bp_try_select[, 
                                                                            1])), ]
                  bp_add <- bp_try_select
                  
                  # adding new base points to base line
                  baseline_new <- rbind(baseline_extended, bp_add)
                  baseline_new_order <- baseline_new[order(baseline_new$X), 
                                                     ]
                }
                }
                if (nrow(bp_add) != 0) {
                  baseline_new <- rbind(baseline_extended, bp_add)
                } else {
                  baseline_new <- baseline_extended
                }
                baseline_new_order <- baseline_new[order(baseline_new$X), ]
                
      
                # interpolate new base points over profile
                baseline_new_int <- approx(baseline_new_order$X, baseline_new_order$Z, 
                  xout = BEP$X)
                baseline_new_int_df <- data.frame(baseline_new_int$x, baseline_new_int$y)
                
                baseline_new2 <- data.frame()
                
                # in some points the base line might still cross the profile in these
                # points the baseline is projected on the profile
                for (i in 1:nrow(profile)) {
                  Differenz <- profile[i, 2] - baseline_new_int_df[i, 2]
                  # if(is.na(Differenz[i])){next}
                  if (Differenz < 0) {
                    base_x <- profile[i, 1]
                    base_y <- profile[i, 2]
                    temp <- data.frame(base_x, base_y)
                    colnames(temp) <- c("X", "Z")
                    baseline_new2 <- rbind(baseline_new2, temp)
                  }
                }

                baseline_scale_2 <- rbind(baseline_new2, baseline_extended, 
                  bp_add)
                baseline_scale_2 <- baseline_scale_2[order(baseline_scale_2$X),]

                plot(BEP[,1:2], type="l", xlab="x [m]", ylab="z [m]")
                lines(baseline_scale_2,col="blue")
                legend("topright", legend=c("baseline layer 1", "baseline layer 2"), col=c("red", "blue"), lty=1, cex=0.7)
            }
        }
        }
    }
      # write profile and baselines to ASCII-file
      
      baseline_scale_1$km[baseline_scale_1$X%in%profile_original[ ,1]] <- profile_original[profile_original[ ,1]%in%baseline_scale_1$X,3]
      baseline_scale_2$km[baseline_scale_2$X%in%profile_original[ ,1]] <- profile_original[profile_original[ ,1]%in%baseline_scale_2$X,3]
      
      write.table(profile_original, paste0(result_path, "/baselines/BEP_iteration_", iteration, ".txt"), row.names = F)
      write.table(baseline_scale_1, paste0(result_path, "/baselines/baseline_scale_1_iteration_", iteration, ".txt"), row.names = F)
      write.table(baseline_scale_2, paste0(result_path, "/baselines/baseline_scale_2_iteration_", iteration, ".txt"), row.names = F)
    
    # call zerocrossing function for scale 3 or function for calculating
    # bedform parameters

    calculate_bedform_parameters(profile_original, baseline_scale_1, 
    baseline_scale_2, window_scale_1, window_scale_2, 
    iteration, zc_threshold, count_scales, result_path, sample_frequency)
}





### calculating bedform parameters for each scale

calculate_bedform_parameters <- function(BEP, baseline_scale_1, 
    baseline_scale_2, window_scale_1, window_scale_2, 
    iteration, zc_threshold, count_scales, result_path, sample_frequency) {
    # if no bedfroms have been detected all bf_parameters are set to NULL
    if (identical(baseline_scale_1, BEP)) {
        baseline_scale_1 <- baseline_scale_2
    }
    if ((identical(baseline_scale_1, BEP) && identical(baseline_scale_2, 
        BEP)) || (is.null(baseline_scale_1))) {
        bf_parameters_scale_1 <- NULL
        bf_parameters_scale_2 <- NULL
    } else {
      
      bl1 <- baseline_scale_1
      bl2 <- baseline_scale_2

      bedforms_l1 <- data.frame()
      for(j in 2:nrow(bl1)){
        bf <- BEP[BEP[ ,1]>=bl1[j-1,1]&BEP[ ,1]<=bl1[j,1], ]
        if(nrow(bf)<3) next
        trend <- approx(bf[c(1,nrow(bf)),1:2], xout = bf[ ,1])
        bf[ ,2] <- bf[ ,2]-trend$y
        h1 <- max(bf[ ,2], na.rm=T)
        l1 <- bf[nrow(bf),1]-bf[1,1,]
        center <- mean(bf[ ,3], na.rm=T)
        npts <- nrow(bf)
        poly_points <- (rbind(bf[,1:2],bf[1,1:2]))
        area <- abs(polyarea(poly_points[ ,1], poly_points[ ,2]))
        shpf <- area/(0.5*h1*l1)
        bedforms_l1 <- rbind(bedforms_l1, data.frame(x_pos=center, l=l1, h=h1, points=npts, area=area, shape=shpf, layer=1, iteration, windowsize=window_scale_1*sample_frequency, zc_threshold))
      }
      bedforms_l1 <- bedforms_l1[bedforms_l1$h>0, ]
          
      if (!is.na(window_scale_2)) {
        bl1_int <- data.frame(approx(bl1[ ,1:2], xout=BEP[ ,1]))
        bl1_int_km <- approx(bl1[ ,c(1,3)], xout=BEP[ ,1])
        bl1_int$km <- bl1_int_km$y
        
        bedforms_l2 <- data.frame()
        for(j in 2:nrow(bl2)){
          bf <- bl1_int[bl1_int[ ,1]>=bl2[j-1,1]&bl1_int[ ,1]<=bl2[j,1], ]
          if(nrow(bf)<3) next
          trend <- approx(bf[c(1,nrow(bf)),1:2], xout = bf[ ,1])
          bf[ ,2] <- bf[ ,2]-trend$y
          h2 <- max(bf[ ,2], na.rm=T)
          l2 <- bf[nrow(bf),1]-bf[1,1,]
          center <- mean(bf[ ,3], na.rm=T)
          npts <- nrow(bf)
          poly_points <- (rbind(bf[,1:2],bf[1,1:2]))
          area <- abs(polyarea(poly_points[ ,1], poly_points[ ,2]))
          shpf <- area/(0.5*h2*l2)
          bedforms_l2 <- rbind(bedforms_l2, data.frame(x_pos=center, l=l2, h=h2, points=npts, area=area, shape=shpf, layer=2, iteration, windowsize=window_scale_2*sample_frequency, zc_threshold))
        }
        bedforms_l2 <- bedforms_l2[bedforms_l2$h>0, ]
      }
          write.table(as.data.frame(bedforms_l1), paste0(result_path, "/statistics/statistics_scale_1_iteration_", iteration, ".txt"), row.names = F)
          if (!is.na(window_scale_2)) {
              write.table(bedforms_l2, paste0(result_path, "/statistics/statistics_scale_2_iteration_", iteration, ".txt"), row.names = F)
          }
    }
}

