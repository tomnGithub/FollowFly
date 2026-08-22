flytrack_frame <-
  function(sub_vial_day_data = data.frame(),
           utilize_ghost = 1,
           ghost_list = list(),
           flynum = 0) {
    nam_day <- unique(sub_vial_day_data$day)[1]
    nam_vial <- unique(sub_vial_day_data$vial)[1]
    nam_pos <- unique(sub_vial_day_data$pos)[1]
    nam_batch <- unique(sub_vial_day_data$batch)[1]

    frame_split <-
      split(sub_vial_day_data, f = sub_vial_day_data$frame)
    vial_max_y <- max(sub_vial_day_data$y) * 0.95
    for (i_frame in (1:(length(frame_split) - 1))) {
      initial_frame <- frame_split[[i_frame]]

      if (utilize_ghost == 1 && i_frame > 5) {
        ghost_list <- c(ghost_list, list(initial_frame))

      }
      if (utilize_ghost == 1 && length(ghost_list) > 10) {
        ghost_list <- ghost_list[-1]
      }
      subsequent_frame <- frame_split[[i_frame + 1]]

      if (i_frame == 1) {
        for (ind_first_frame in (1:nrow(initial_frame))) {
          flynum <- flynum + 1
          if (flynum < 10) {
            initial_frame$name[ind_first_frame] <-
              paste0(
                "fly_0",
                as.character(flynum),
                "_",
                nam_day,
                nam_vial,
                nam_pos,
                nam_batch,
                paste0(sample(
                  c(0:9, letters, LETTERS), 6, replace = TRUE
                ), collapse = "")
              )
          } else{
            initial_frame$name[ind_first_frame] <-
              paste0(
                "fly_",
                as.character(flynum),
                "_",
                nam_day,
                nam_vial,
                nam_pos,
                nam_batch,
                paste0(sample(
                  c(0:9, letters, LETTERS), 6, replace = TRUE
                ), collapse = "")
              )
          }

          frame_split[[i_frame]]$name[ind_first_frame] <-
            initial_frame$name[ind_first_frame]
        }
      }

      if (i_frame < 10 || utilize_ghost == 0) {
        fly_dist <- .point_distance(
          subsequent_frame[, c("x", "y")],
          initial_frame[, c("x", "y")]
        )

        if (nrow(initial_frame) < 2) {

          subsequent_frame$distance = fly_dist
          subsequent_frame$SRC_ID = initial_frame$name

        } else{
          if (nrow(subsequent_frame) < 2) {
            subsequent_frame$distance = min(fly_dist)
            subsequent_frame$SRC_ID = initial_frame$name[which.min(fly_dist)]
          } else{
            min_dist_index <-
              apply(fly_dist, 1, which.min)
            subsequent_frame$distance = fly_dist[cbind(1:nrow(fly_dist), min_dist_index)]
            subsequent_frame$SRC_ID = initial_frame$name[min_dist_index]
          }

        }
      }
      if (utilize_ghost == 1 && i_frame >= 10) {
        ghost_df_all <-
          bind_rows(ghost_list)

        ghost_df <- ghost_df_all %>%
          group_by(name) %>%
          filter(frame == max(frame)) %>%
          ungroup()
        ghost_df <- ghost_df %>%
          group_by(name) %>%
          filter(distance == min(distance)) %>%
          ungroup()

        ghost_names <- unique(ghost_df$name)

        if (i_frame >= 10) {
          pred_ghost_df <- ghost_df
          for (pred_pos in 1:length(ghost_names)) {
            ghost_df_name_sub <-
              subset(ghost_df, name == ghost_names[pred_pos])

            recent_frame <- unique(initial_frame$frame)

            re_frame_diff <- recent_frame - ghost_df_name_sub$frame

            if (re_frame_diff > 0) {
              predict_x <-
                ghost_df_name_sub$x + (mean(ghost_df_all$xDiff[ghost_df_all$name == ghost_names[pred_pos]]) * re_frame_diff)
              predict_y <-
                ghost_df_name_sub$y + (mean(ghost_df_all$yDiff[ghost_df_all$name == ghost_names[pred_pos]]) * re_frame_diff)

              ghost_df$predX[ghost_df$name == ghost_names[pred_pos]] <-
                predict_x
              ghost_df$predY[ghost_df$name == ghost_names[pred_pos]] <-
                predict_y
              ghost_df$isPred[ghost_df$name == ghost_names[pred_pos]] <-
                1

            } else{
              ghost_df$predX[ghost_df$name == ghost_names[pred_pos]] <-
                ghost_df_name_sub$x
              ghost_df$predY[ghost_df$name == ghost_names[pred_pos]] <-
                ghost_df_name_sub$y
            }
            pred_ghost_df$xDiff[pred_ghost_df$name == ghost_names[pred_pos]] <-
              mean(ghost_df_all$xDiff[ghost_df_all$name == ghost_names[pred_pos]])
            pred_ghost_df$yDiff[pred_ghost_df$name == ghost_names[pred_pos]] <-
              mean(ghost_df_all$yDiff[ghost_df_all$name == ghost_names[pred_pos]])

          }
          pred_ghost_df$predX <- ghost_df$predX
          pred_ghost_df$predY <- ghost_df$predY
          pred_ghost_df$isPred <- ghost_df$isPred

        }
        if (i_frame < 11) {
          ghost_df <- as.data.frame(ghost_df)
          ghost_result <-
            matchPairs(subsequent_frame, ghost_df, 0, vial_max_y)

        } else{
          pred_ghost_df <- as.data.frame(pred_ghost_df)
          ghost_result <-
            matchPairs(subsequent_frame, pred_ghost_df, 1, vial_max_y)

        }

        subsequent_frame$distance <-
          ghost_result[[1]]
        subsequent_frame$SRC_ID <- ghost_result[[2]]

        subsequent_frame$reorder <-
          ghost_result[[3]]
        subsequent_frame$name <- subsequent_frame$SRC_ID

      }

      exceeded_distance_rows <-
        row.names(subsequent_frame[subsequent_frame$distance > 20, ])
      if (length(exceeded_distance_rows) > 0) {
        for (rename_exceeded_it in 1:length(exceeded_distance_rows)) {
          flynum <- flynum + 1
          if (flynum < 10) {
            subsequent_frame[exceeded_distance_rows[rename_exceeded_it], ]$SRC_ID <-
              paste0(
                "fly_0",
                as.character(flynum),
                "_",
                nam_day,
                nam_vial,
                nam_pos,
                nam_batch,
                paste0(sample(
                  c(0:9, letters, LETTERS), 6, replace = TRUE
                ), collapse = "")
              )
          } else{
            subsequent_frame[exceeded_distance_rows[rename_exceeded_it], ]$SRC_ID <-
              paste0(
                "fly_",
                as.character(flynum),
                "_",
                nam_day,
                nam_vial,
                nam_pos,
                nam_batch,
                paste0(sample(
                  c(0:9, letters, LETTERS), 6, replace = TRUE
                ), collapse = "")
              )
          }
          subsequent_frame[exceeded_distance_rows[rename_exceeded_it], ]$named_this_frame <-
            1
        }
      }

      for (rename_iteration in (1:nrow(subsequent_frame))) {
        if (sum(subsequent_frame$SRC_ID == subsequent_frame$SRC_ID[rename_iteration]) >
            1) {
          sub_identical_frames <-
            subsequent_frame[subsequent_frame$SRC_ID == subsequent_frame$SRC_ID[rename_iteration], ]

          if (i_frame < 10) {
            original_y <-
              initial_frame$y[initial_frame$name == subsequent_frame$SRC_ID[rename_iteration]]
            original_X <-
              initial_frame$x[initial_frame$name == subsequent_frame$SRC_ID[rename_iteration]]
          } else{
            original_X <-
              ghost_df$x[ghost_df$name == subsequent_frame$SRC_ID[rename_iteration]]
            original_y <-
              ghost_df$y[ghost_df$name == subsequent_frame$SRC_ID[rename_iteration]]
          }

          sub_identical_frames$original_X_diff <-
            abs(sub_identical_frames$x - original_X)
          sub_identical_frames$original_y_diff <-
            abs(sub_identical_frames$y - original_y)

          sub_identical_frames$weighted_diff <-
            (sub_identical_frames$original_X_diff + 1) ^ 1.1 + (sub_identical_frames$original_y_diff)
          smallest_row <-
            row.names(sub_identical_frames[which.min(sub_identical_frames$weighted_diff), ])
          all_row_names <- row.names(sub_identical_frames)
          largest_rows <-
            all_row_names[!(all_row_names %in% smallest_row)]
          for (rename_identical_it in 1:length(largest_rows)) {
            flynum <- flynum + 1
            if (flynum < 10) {
              subsequent_frame[largest_rows[rename_identical_it], ]$SRC_ID <-
                paste0(
                  "fly_0",
                  as.character(flynum),
                  "_",
                  nam_day,
                  nam_vial,
                  nam_pos,
                  nam_batch,
                  paste0(sample(
                    c(0:9, letters, LETTERS), 6, replace = TRUE
                  ), collapse = "")
                )
            } else{
              subsequent_frame[largest_rows[rename_identical_it], ]$SRC_ID <-
                paste0(
                  "fly_",
                  as.character(flynum),
                  "_",
                  nam_day,
                  nam_vial,
                  nam_pos,
                  nam_batch,
                  paste0(sample(
                    c(0:9, letters, LETTERS), 6, replace = TRUE
                  ), collapse = "")
                )
            }
            subsequent_frame[largest_rows[rename_identical_it], ]$named_this_frame <-
              1
          }

        }

        if (utilize_ghost == 1 &&
            i_frame >= 10 &&
            subsequent_frame$named_this_frame[rename_iteration] < 1) {
          frame_diff <-
            subsequent_frame$frame[rename_iteration] - ghost_df$frame[ghost_df$name == subsequent_frame$SRC_ID[rename_iteration]]
          subsequent_frame$xDiff[rename_iteration] = (subsequent_frame$x[rename_iteration] -
                                                        ghost_df$x[ghost_df$name == subsequent_frame$SRC_ID[rename_iteration]]) / frame_diff
          subsequent_frame$yDiff[rename_iteration] = (subsequent_frame$y[rename_iteration] -
                                                        ghost_df$y[ghost_df$name == subsequent_frame$SRC_ID[rename_iteration]]) / frame_diff
        }
      }

      frame_split[[i_frame + 1]]$name = subsequent_frame$SRC_ID
      frame_split[[i_frame + 1]]$distance = subsequent_frame$distance
      frame_split[[i_frame + 1]]$named_this_frame = subsequent_frame$named_this_frame
      frame_split[[i_frame + 1]]$reorder = subsequent_frame$reorder
      if (utilize_ghost == 1 && i_frame >= 10) {
        frame_split[[i_frame + 1]]$xDiff = subsequent_frame$xDiff
        frame_split[[i_frame + 1]]$yDiff = subsequent_frame$yDiff
      }

    }
    frames_tracked_df = bind_rows(frame_split)
    return(frames_tracked_df)
  }
