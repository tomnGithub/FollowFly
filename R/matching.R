.point_distance <- function(p1, p2) {
  dx <- outer(p1$x, p2$x, "-")
  dy <- outer(p1$y, p2$y, "-")
  sqrt(dx^2 + dy^2)
}

calculate_3d_distance <- function(p1, p2) {
  sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2 + (p1$frame - p2$frame)^2)
}


matchPairs <-function(subsequent_frame_mp,
                      ghost_frame_mp,
                      ghost,
                      vial_max_y) {
  if (ghost == 1) {
    prediction_data <-ghost_frame_mp[, c("predX", "predY", "frame")]
    subsequent_frame_3d <- subsequent_frame_mp[, c("x", "y", "frame")]
    ghost_frame_3d <- ghost_frame_mp[, c("x", "y", "frame")]

    colnames(prediction_data) <- c("x", "y", "frame")

    prediction_data$frame <- max(prediction_data$frame)

    fly_dist_mp <-matrix(NA,
                         nrow = nrow(subsequent_frame_3d),
                         ncol = nrow(prediction_data))
    for (initial_frame_index in seq_len(nrow(subsequent_frame_3d))) {
      for (prediction_index in seq_len(nrow(prediction_data))) {
        fly_dist_mp[initial_frame_index, prediction_index] <-calculate_3d_distance(subsequent_frame_3d[initial_frame_index, ], prediction_data[prediction_index, ])
      }
    }

    fly_dist_mp_actual <-
      matrix(NA,
             nrow = nrow(subsequent_frame_3d),
             ncol = nrow(ghost_frame_3d))
    for (initial_frame_index_s in seq_len(nrow(subsequent_frame_3d))) {
      for (subsequent_frame_index in seq_len(nrow(ghost_frame_3d))) {
        fly_dist_mp_actual[initial_frame_index_s, subsequent_frame_index] <-calculate_3d_distance(subsequent_frame_3d[initial_frame_index_s, ], ghost_frame_3d[subsequent_frame_index, ])
      }
    }

  } else{
    fly_dist_mp <- .point_distance(
      subsequent_frame_mp[, c("x", "y")],
      ghost_frame_mp[, c("x", "y")]
    )
  }

  if (nrow(ghost_frame_mp) < 2) {
    subsequent_frame_mp$distance = fly_dist_mp
    subsequent_frame_mp$SRC_ID = ghost_frame_mp$name

  } else{
    if (nrow(subsequent_frame_mp) < 2) {
      subsequent_frame_mp$distance = min(fly_dist_mp)
      subsequent_frame_mp$SRC_ID = ghost_frame_mp$name[which.min(fly_dist_mp)]
    } else{
      min_dist_index_mp <- apply(fly_dist_mp, 1, which.min)
      subsequent_frame_mp$distance = fly_dist_mp[cbind(1:nrow(fly_dist_mp), min_dist_index_mp)]
      subsequent_frame_mp$SRC_ID = ghost_frame_mp$name[min_dist_index_mp]

      if (ghost == 1) {
        min_dist_index_mp <- apply(fly_dist_mp_actual, 1, which.min)
        subsequent_frame_mp$distance_ac <-fly_dist_mp_actual[cbind(1:nrow(fly_dist_mp_actual), min_dist_index_mp)]
        subsequent_frame_mp$SRC_ac <- ghost_frame_mp$name[min_dist_index_mp]

        for (match_iteration in 1:nrow(subsequent_frame_mp)) {
          do_not <- 0

          if (subsequent_frame_mp$SRC_ac[match_iteration] != subsequent_frame_mp$SRC_ID[match_iteration]) {
            if (subsequent_frame_mp$distance_ac[match_iteration] < subsequent_frame_mp$distance[match_iteration]) {
              subsequent_frame_mp$distance[match_iteration] <- subsequent_frame_mp$distance_ac[match_iteration]
              subsequent_frame_mp$SRC_ID[match_iteration] <-subsequent_frame_mp$SRC_ac[match_iteration]
              do_not <- 1
            }

          }
          if(sum(subsequent_frame_mp$SRC_ID == subsequent_frame_mp$SRC_ID[match_iteration])>1){

            dist_ghost_frame <- cbind(ghost_frame_mp, dist_to_point = fly_dist_mp[match_iteration, ])
            sub_identical_frames_mp <- subsequent_frame_mp[subsequent_frame_mp$SRC_ID == subsequent_frame_mp$SRC_ID[match_iteration],]

            original_X<-  ghost_frame_mp$x[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]]
            original_y<-  ghost_frame_mp$y[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]]

            sub_identical_frames_mp$original_X_diff<-abs(sub_identical_frames_mp$x- original_X)
            sub_identical_frames_mp$original_y_diff<-abs(sub_identical_frames_mp$y- original_y)

            sub_identical_frames_mp$weighted_diff<- (sub_identical_frames_mp$original_X_diff + 1 )^1.1  + (sub_identical_frames_mp$original_y_diff)
            sub_identical_frames_mp <- sub_identical_frames_mp[order(sub_identical_frames_mp$weighted_diff), ]

            smallest_row_mp<-row.names(sub_identical_frames_mp[which.min(sub_identical_frames_mp$weighted_diff),])
            all_row_names_mp <- row.names(sub_identical_frames_mp)
            largest_rows_mp <- all_row_names_mp[!(all_row_names_mp %in% smallest_row_mp)]
            for(rename_identical_it_mp in 1:length(largest_rows_mp)){

              dist_ghost_frame <- cbind(ghost_frame_mp, dist_to_point = fly_dist_mp[match_iteration, ])
              subsequent_frame_mp[largest_rows_mp[rename_identical_it_mp],]$reorder <- paste0("0_", sub_identical_frames_mp[smallest_row_mp,]$SRC_ID)
              ident_shortest_distance <- order(fly_dist_mp[match_iteration,])[1+rename_identical_it_mp]
              subsequent_frame_mp[largest_rows_mp[rename_identical_it_mp],]$SRC_ID <-ghost_frame_mp$name[ident_shortest_distance]
              subsequent_frame_mp[largest_rows_mp[rename_identical_it_mp],]$distance <- dist_ghost_frame$dist_to_point[dist_ghost_frame$name ==ghost_frame_mp$name[ident_shortest_distance]]

              subsequent_frame_mp[largest_rows_mp[rename_identical_it_mp],]$isDuplicate<-1

            }

          }

          if (subsequent_frame_mp$y[match_iteration] < (vial_max_y) &&
              (subsequent_frame_mp$y[match_iteration] - ghost_frame_mp$y[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]]) < 0 &&
              (subsequent_frame_mp$frame[match_iteration] - ghost_frame_mp$frame[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]]) > 5 &&
              do_not != 1 &&
              subsequent_frame_mp$isDuplicate[match_iteration]<1 ) {
            subsequent_frame_mp$reorder[match_iteration] <- paste0("1_", ghost_frame_mp$name[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]],"_",subsequent_frame_mp$reorder[match_iteration])
            shortest_distance <- order(fly_dist_mp[match_iteration,])[2]
            subsequent_frame_mp$SRC_ID[match_iteration] <- ghost_frame_mp$name[shortest_distance]
            subsequent_frame_mp$distance[match_iteration] <- fly_dist_mp[match_iteration, shortest_distance]

          }

          frames_are_diff <-subsequent_frame_mp$frame[match_iteration] - ghost_frame_mp$frame[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]]

          if (ghost_frame_mp$named_this_frame[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]] == 1 &&
              frames_are_diff > 3  &&
              subsequent_frame_mp$isDuplicate[match_iteration]<1) {
            dist_ghost_frame <- cbind(ghost_frame_mp, dist_to_point = fly_dist_mp[match_iteration, ])
            dist_ghost_frame_named<-subset(dist_ghost_frame, named_this_frame==0)
            if(nrow(dist_ghost_frame_named) > 0){
              shortest_distance <- order(dist_ghost_frame_named$dist_to_point)[1]
              subsequent_frame_mp$reorder[match_iteration] <- paste0("3_", ghost_frame_mp$name[ghost_frame_mp$name == subsequent_frame_mp$SRC_ID[match_iteration]], "_",subsequent_frame_mp$reorder[match_iteration])
              subsequent_frame_mp$SRC_ID[match_iteration]<-dist_ghost_frame_named$name[shortest_distance]
              subsequent_frame_mp$distance[match_iteration]<-dist_ghost_frame_named$dist_to_point[shortest_distance]
            }else{

            }

          }

        }
      }

    }

  }

  result <-list(
    subsequent_frame_mp$distance,
    subsequent_frame_mp$SRC_ID,
    subsequent_frame_mp$reorder
  )
  return(result)
}
