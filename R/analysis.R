lm_func <- function(data, c_of_i) {
  n <- nrow(data)
  window_size <- 30
  if (n < window_size) {
    window_size <- n
  }

  num_windows <- n - window_size + 1
  total_num_coi <- length(unique(c_of_i))
  slope_df <-
    data.frame(
      day = factor(),
      name = character(),
      vial = factor(),
      pos = factor(),
      batch = factor(),
      stringsAsFactors = FALSE
    )
  for (new_col_i in 1:total_num_coi) {
    slope_df[[paste0(c_of_i[new_col_i], "_startFrame")]] <- numeric()
    slope_df[[paste0(c_of_i[new_col_i], "_endFrame")]] <- numeric()
    slope_df[[paste0(c_of_i[new_col_i], "_slope")]] <- numeric()
    slope_df[[paste0(c_of_i[new_col_i], "_rsq")]] <- numeric()
  }

  for (num_coi in 1:(total_num_coi)) {

    for (i in 1:num_windows) {
      window_data <-
        data[i:(i + window_size - 1),]
      model <-
        stats::lm(paste0(c_of_i[num_coi], "~", "frame"), data = window_data)
      slope <- stats::coef(model)["frame"]
      rsquare <- summary(model)$r.squared
      startFrame <- data$frame[i]
      endFrame <- data$frame[(i + window_size - 1)]

      if (num_coi > 1) {
        slope_df[[paste0(c_of_i[num_coi], "_startFrame")]][i] <- startFrame
        slope_df[[paste0(c_of_i[num_coi], "_endFrame")]][i] <-
          endFrame
        slope_df[[paste0(c_of_i[num_coi], "_slope")]][i] <- slope
        slope_df[[paste0(c_of_i[num_coi], "_rsq")]][i] <- rsquare
      } else{
        row <-
          data.frame(
            day = window_data$day[1],
            name = unique(window_data$name),
            vial = window_data$vial[1],
            pos = unique(window_data$pos),
            batch = unique(window_data$batch),
            stringsAsFactors = FALSE
          )
        for (new_col_j in 1:total_num_coi) {
          row[[paste0(c_of_i[new_col_j], "_startFrame")]] <-
            as.numeric(NA)
          row[[paste0(c_of_i[new_col_j], "_endFrame")]] <-
            as.numeric(NA)
          row[[paste0(c_of_i[new_col_j], "_slope")]] <- as.numeric(NA)
          row[[paste0(c_of_i[new_col_j], "_rsq")]] <- as.numeric(NA)
        }
        row[[paste0(c_of_i[num_coi], "_startFrame")]] <- startFrame
        row[[paste0(c_of_i[num_coi], "_endFrame")]] <- endFrame
        row[[paste0(c_of_i[num_coi], "_slope")]] <- slope
        row[[paste0(c_of_i[num_coi], "_rsq")]] <- rsquare

        slope_df <- rbind(slope_df, row)
      }
    }

  }
  return(slope_df)
}

Flytrack_high_y <- function(framelimit_dataset_insert,
                            current_name) {
  y_limit <- framelimit_dataset_insert$y_limit[1]
  framelimit_dataset_insert$High_Y[framelimit_dataset_insert$name == current_name &
                                     framelimit_dataset_insert$y > y_limit] <- 1
  frame_limit <-
    framelimit_dataset_insert$frame[framelimit_dataset_insert$name == current_name &
                                      framelimit_dataset_insert$High_Y == 1][2]
  framelimit_dataset_insert$High_Y[framelimit_dataset_insert$name == current_name &
                                     framelimit_dataset_insert$frame > frame_limit] <- 1
  return(framelimit_dataset_insert)
}

FlyTrack_Analysis <-
  function(smoothed_tracks = data.frame(),
           usr_framelimit = 30,
           var_for_slope = c("kalman_y"),
           limit_late_flies = 0,
           frame_include_limit = 90,
           filter_negative_slope = 1,
           filter_rsq_slope = 1,
           rsq_limit = 0.8,
           quiet = FALSE) {
    framelimit_dataset_insert <- smoothed_tracks %>%
      group_by(name) %>%
      filter(max(frame) - min(frame) >= usr_framelimit)

    framelimit_dataset_insert <- framelimit_dataset_insert %>%
      group_by(day, vial) %>%
      mutate(y_limit = max(y) * 0.95)
    framelimit_dataset_insert <-
      as.data.frame(framelimit_dataset_insert)

    if (!quiet) message("Identifying where flies stopped climbing.")
    framelimit_dataset_insert$High_Y <- 0
    high_y_fun <- Flytrack_high_y
    framelimit_dataset_insert_highY <-
      foreach(
        name_i = 1:length(unique(framelimit_dataset_insert$name)),
        .combine = 'rbind',
        .packages = c("FlyTrackPACK"),
        .export = "high_y_fun"
      ) %dopar% {
        result <-
          high_y_fun(
            subset(
              framelimit_dataset_insert,
              name == unique(framelimit_dataset_insert$name)[name_i]
            ),
            unique(framelimit_dataset_insert$name)[name_i]
          )

      }
    if (!quiet) message("All stops identified.")

    data_y_and_frame_limited <-
      subset(framelimit_dataset_insert_highY, High_Y  < 1)

    if (limit_late_flies == 1) {
      data_y_and_frame_limited$y_late_limit <- 200

      data_y_and_frame_limited <- data_y_and_frame_limited %>%
        group_by(name) %>%
        filter(!(
          any(
            first_point == 1 & frame > frame_include_limit & y > y_late_limit
          )
        ))

      data_y_and_frame_limited <-
        as.data.frame(data_y_and_frame_limited)
      if (!quiet) message("Late flies filtered out.")
    }

    if (!quiet) message("Starting moving slope window.")
    slope_fun <- lm_func
    individual_slopes <-
      foreach(
        slope_i = 1:length(unique(data_y_and_frame_limited$name)),
        .combine = 'rbind',
        .packages = c("FlyTrackPACK"),
        .export = "slope_fun"
      ) %dopar% {
        result <-
          slope_fun(subset(
            data_y_and_frame_limited,
            name == unique(data_y_and_frame_limited$name)[slope_i]
          ),
          var_for_slope)

      }
    if (!quiet) message("Finished moving slope window.")

    max_slope_list <- list()

    for (slope_vars in 1:length(unique(var_for_slope))) {
      combined_subset <- cbind(individual_slopes[, 1:5],
                               individual_slopes[, (2 + (slope_vars * 4)):(5 + (slope_vars * 4))])
      max_slope_df <- combined_subset %>%
        group_by(day, name) %>%
        filter(eval(parse(text = paste0(
          var_for_slope[slope_vars], "_rsq"
        ))) == max((eval(
          parse(text = paste0(var_for_slope[slope_vars], "_rsq"))
        )))) %>%
        ungroup()

      if (filter_negative_slope == 1 && filter_rsq_slope == 1) {
        max_slope_df <-
          subset(max_slope_df, max_slope_df[[paste0(var_for_slope[slope_vars], "_slope")]] >= 0)
        max_slope_df <-
          subset(max_slope_df, max_slope_df[[paste0(var_for_slope[slope_vars], "_rsq")]] >= rsq_limit)
      } else if (filter_negative_slope == 1 && filter_rsq_slope == 0) {
        max_slope_df <-
          subset(max_slope_df, max_slope_df[[paste0(var_for_slope[slope_vars], "_slope")]] >= 0)
      } else if (filter_negative_slope == 0 && filter_rsq_slope == 1) {
        max_slope_df <-
          subset(max_slope_df, max_slope_df[[paste0(var_for_slope[slope_vars], "_rsq")]] >= rsq_limit)
      }
      max_slope_list[[paste0(var_for_slope[slope_vars], "_max_slope_df")]] <-
        max_slope_df

    }

    for (msl_entry_full in 1:length(max_slope_list)) {
      if (msl_entry_full == 1) {
        final_slope_dataset <- max_slope_list[[1]]
      } else{
        suppressWarnings(
          final_slope_dataset <-
            merge(
              final_slope_dataset,
              max_slope_list[[msl_entry_full]],
              by = "name",
              all.x = TRUE,
              all.y = TRUE,
              suffixes = c("", "")
            )
        )
        final_slope_dataset <-
          final_slope_dataset[,-c((2 + ((
            msl_entry_full * 4
          ))):(5 + ((
            msl_entry_full * 4
          ))))]

      }
    }
    return(final_slope_dataset)

  }
