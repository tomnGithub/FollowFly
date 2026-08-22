run_flytrack <- function(trackedDF,
                         workers = NULL,
                         seed = NULL,
                         var_for_slope = c("kalman_y"),
                         rsq_limit = 0.8,
                         filter_negative_slope = 1,
                         filter_rsq_slope = 1,
                         limit_late_flies = 0,
                         frame_include_limit = 90,
                         quiet = FALSE) {
  required_columns <- c("x", "y", "frame", "vialID", "day", "pos", "batch")
  missing_columns <- setdiff(required_columns, names(trackedDF))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  trackedDF$day <- as.numeric(gsub("\\D", "", trackedDF$day))
  trackedDF$nearest_neighbor <- NA
  trackedDF$isDuplicate <- 0
  trackedDF <- prepData_flytrack(trackedDF)

  if (is.null(workers)) {
    cores <- parallel::detectCores()
    if (is.na(cores)) cores <- 2L
    workers <- max(1L, cores - 1L)
  }
  workers <- max(1L, as.integer(workers))

  cl <- parallel::makeCluster(workers, outfile = "")
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, iseed = seed)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  doParallel::registerDoParallel(cl)

  start_time <- Sys.time()
  duplicate_fun <- duplicateFilter_Frame
  duplicate_groups <- split(trackedDF, interaction(trackedDF$day, trackedDF$vial, drop = TRUE))
  duplicate_filtered_trackedDF <- foreach::foreach(
    group = duplicate_groups,
    .combine = "rbind",
    .packages = c("dplyr", "FlyTrackPACK"),
    .export = "duplicate_fun"
  ) %dopar% {
    duplicate_fun(group)
  }
  rownames(duplicate_filtered_trackedDF) <- NULL

  duplicate_time <- Sys.time()
  if (!quiet) {
    message(sprintf(
      "Duplicate filtering complete (%.2f min).",
      as.numeric(difftime(duplicate_time, start_time, units = "mins"))
    ))
  }

  data <- prepData_flytrack(duplicate_filtered_trackedDF)
  track_fun <- flytrack_frame
  tracking_groups <- split(data, interaction(data$day, data$vial, drop = TRUE))
  final_base_trackedDF <- foreach::foreach(
    group = tracking_groups,
    .combine = "rbind",
    .packages = c("dplyr", "FlyTrackPACK"),
    .export = "track_fun"
  ) %dopar% {
    track_fun(group, 1)
  }
  rownames(final_base_trackedDF) <- NULL

  tracking_time <- Sys.time()
  if (!quiet) {
    message(sprintf(
      "Fly identification complete (%.2f min).",
      as.numeric(difftime(tracking_time, duplicate_time, units = "mins"))
    ))
  }

  final_base_trackedDF$inserted_prediction <- 0
  final_base_trackedDF$kalman_x <- final_base_trackedDF$x
  final_base_trackedDF$kalman_y <- final_base_trackedDF$y
  final_base_trackedDF$first_point <- 0

  kalman_fun <- applyKalman_flytrack
  kalman_groups <- split(final_base_trackedDF, final_base_trackedDF$name)
  final_Kalman_trackedDF <- foreach::foreach(
    group = kalman_groups,
    .combine = "rbind",
    .packages = c("dplyr", "FlyTrackPACK"),
    .export = "kalman_fun"
  ) %dopar% {
    kalman_fun(group)
  }
  rownames(final_Kalman_trackedDF) <- NULL

  kalman_time <- Sys.time()
  if (!quiet) {
    message(sprintf(
      "Kalman smoothing complete (%.2f min).",
      as.numeric(difftime(kalman_time, tracking_time, units = "mins"))
    ))
  }

  final_slope_dataset <- FlyTrack_Analysis(
    final_Kalman_trackedDF,
    var_for_slope = var_for_slope,
    limit_late_flies = limit_late_flies,
    frame_include_limit = frame_include_limit,
    filter_negative_slope = filter_negative_slope,
    filter_rsq_slope = filter_rsq_slope,
    rsq_limit = rsq_limit,
    quiet = quiet
  )

  if (!quiet) {
    end_time <- Sys.time()
    message(sprintf(
      "Slope analysis complete. Total runtime: %.2f min.",
      as.numeric(difftime(end_time, start_time, units = "mins"))
    ))
  }

  list(
    duplicate_filtered = duplicate_filtered_trackedDF,
    base_tracks = final_base_trackedDF,
    kalman_tracks = final_Kalman_trackedDF,
    slopes = final_slope_dataset
  )
}
