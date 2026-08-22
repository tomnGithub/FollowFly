my_kalman_filter<-function(data_acc){
  {
    dt <- 1/30
    u_x<-1
    u_y<- 1
    std_acc<- 10
    x_std_meas<- 0.1
    y_std_meas<- 0.1

    u <- matrix(c(u_x, u_y),nrow=2)
    x <- matrix(c(0, 0, 0, 0), ncol = 1)
    A <- matrix(c(1, 0, dt, 0, 0, 1, 0, dt,  0, 0, 1, 0, 0, 0, 0, 1), ncol=4, nrow = 4, byrow = TRUE)
    B <- matrix(c((dt^2)/2, 0, 0,(dt^2)/2,dt, 0,0, dt),ncol=2, nrow = 4, byrow = TRUE)
    H <- matrix(c(1, 0, 0, 0,
                  0, 1, 0, 0), ncol=4,nrow = 2, byrow = TRUE)
    Q <- matrix(c((dt^4)/4, 0, (dt^3)/2, 0,
                  0, (dt^4)/4, 0, (dt^3)/2,
                  (dt^3)/2, 0, dt^2, 0,
                  0, (dt^3)/2, 0, dt^2) * std_acc^2, ncol=4, nrow = 4, byrow = TRUE)
    R <- matrix(c(x_std_meas^2, 0,
                  0, y_std_meas^2), ncol = 2 ,nrow = 2, byrow = TRUE)
    P <- diag(nrow(A))
  }

  empty_df <- data.frame(matrix(nrow = nrow(data_acc), ncol = 6))

  for(kal_i in 1:nrow(data_acc)){

    if(!is.na(data_acc$x[kal_i])){

      z<-matrix(c(data_acc$x[kal_i], data_acc$y[kal_i]),nrow=2)
      S <- H %*% P %*% t(H) + R

      K <- P %*% t(H) %*% solve(S)

      x <- (x + K %*% (z - H %*% x))

      I <- diag(ncol(H))

      P <- (I - K %*% H) %*% P

      return_stuff2<-x[1:2, , drop = FALSE]
      estimated_x<-return_stuff2[1]
      estimated_y<-return_stuff2[2]
    }

    x <- A %*% x + B %*% u

    P <- A %*% P %*% t(A) + Q

    return_stuff<-x[0:2, , drop = FALSE]
    predicted_x<-return_stuff[1]
    predicted_y<-return_stuff[2]

    empty_df[kal_i,1]<-data_acc$x[kal_i]
    empty_df[kal_i,2]<-data_acc$y[kal_i]
    empty_df[kal_i,3]<-predicted_x
    empty_df[kal_i,4]<-predicted_y
    empty_df[kal_i,5]<-estimated_x
    empty_df[kal_i,6]<-estimated_y
  }
  {
    empty_df$frame<-data_acc$frame
    colnames(empty_df) <- c("x", "y", "predicted_x", "predicted_y", "estimated_x", "estimated_y", "time")
    return(empty_df)
  }
}

applyKalman_flytrack <- function(tracking_split) {

  if (nrow(tracking_split) > 1) {
    kal_col_names <- colnames(tracking_split)
    kal_num_col <- ncol(tracking_split)
    kal_exclude_insert <-
      c(
        "y",
        "x",
        "frame",
        "t",
        "nearest_neighbor",
        "isDuplicate",
        "distance",
        "xDiff",
        "yDiff",
        "predX",
        "predY",
        "isPred",
        "distance_ac",
        "SRC_ac",
        "named_this_frame",
        "reorder",
        "concatenated_names",
        "distanceNN",
        "inserted_prediction",
        "kalman_x",
        "kalman_y",
        "first_point"
      )

    kal_include_insert <-
      kal_col_names[!(kal_col_names %in% kal_exclude_insert)]
    kal_col_info <- tracking_split[1, ][, kal_include_insert]
    tracked_frame <- 1
    split_data_status <- 0

    while (split_data_status != 1) {
      split_occured <- 0
      split_frame_diff <-
        tracking_split$frame[tracked_frame + 1] - tracking_split$frame[tracked_frame]
      if (split_frame_diff > 1) {

        df_top <- tracking_split[1:(tracked_frame),]
        df_bottom <-
          tracking_split[((tracked_frame + 1):(nrow(tracking_split))),]
        new_rows <-
          data.frame(matrix(ncol = kal_num_col, nrow = (split_frame_diff - 1)))
        colnames(new_rows) <- kal_col_names
        new_rows$frame <-
          ((max(df_top$frame) + 1):(min(df_bottom$frame) - 1))

        new_rows[, kal_include_insert] <- kal_col_info
        new_rows$inserted_prediction <- 1
        new_rows$first_point <- 0
        x_split_diff <- df_bottom$x[1] - df_top$x[nrow(df_top)]
        y_split_diff <- df_bottom$y[1] - df_top$y[nrow(df_top)]
        x_split_dxdt <- x_split_diff / (split_frame_diff)
        y_split_dxdt <- y_split_diff / (split_frame_diff)

        new_rows$x <- as.double(df_top$x[nrow(df_top)])
        new_rows$y <- as.double(df_top$y[nrow(df_top)])
        for (estimate_i in 1:nrow(new_rows)) {
          new_rows$x[estimate_i] <-
            new_rows$x[estimate_i] + (x_split_dxdt * estimate_i)
          new_rows$y[estimate_i] <-
            new_rows$y[estimate_i] + (y_split_dxdt * estimate_i)
        }

        new_df_in <- rbind(df_top, new_rows, df_bottom)
        rownames(new_df_in) <- NULL
        split_occured <- 1
        tracking_split <- new_df_in

      }
      if (split_occured > 0) {
        tracked_frame <- tracked_frame + split_frame_diff
      } else{
        tracked_frame <- tracked_frame + 1
      }

      if (tracked_frame + 1 > nrow(tracking_split)) {
        kalman_results <- my_kalman_filter(tracking_split)
        tracking_split$kalman_x <- kalman_results$estimated_x
        tracking_split$kalman_y <- kalman_results$estimated_y

        tracking_split <- tracking_split %>%
          arrange(name, frame) %>%
          group_by(name) %>%
          mutate(
            x_distance = (x - lag(x, default = first(x))),
            y_distance = (y - lag(y, default = first(y))),
            euclidean_distance = sqrt((x - lag(
              x, default = first(x)
            )) ^ 2 + (y - lag(
              y, default = first(y)
            )) ^ 2),
            cumulative_distance = cumsum(euclidean_distance),
            kalman_x_distance = (kalman_x - lag(kalman_x, default = first(kalman_x))),
            kalman_y_distance = (kalman_y - lag(kalman_y, default = first(kalman_y))),
            kalman_euclidean_distance = sqrt((
              kalman_x - lag(kalman_x, default = first(kalman_x))
            ) ^ 2 + (
              kalman_y - lag(kalman_y, default = first(kalman_y))
            ) ^ 2),
            kalman_cumulative_distance = cumsum(kalman_euclidean_distance),
            kalman_angle = atan2(kalman_y_distance, kalman_x_distance) * (180 / pi)
          )

        tracking_split$first_point[1] <- 1
        split_data_status <- 1

      }

    }

  }

  return(tracking_split)
}
