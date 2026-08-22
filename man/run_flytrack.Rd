\name{run_flytrack}
\alias{run_flytrack}
\title{Run the FlyTrackPACK workflow}
\usage{
run_flytrack(
  trackedDF,
  workers = NULL,
  seed = NULL,
  var_for_slope = c("kalman_y"),
  rsq_limit = 0.8,
  filter_negative_slope = 1,
  filter_rsq_slope = 1,
  limit_late_flies = 0,
  frame_include_limit = 90,
  quiet = FALSE
)
}
\arguments{
  \item{trackedDF}{A FreeClimber raw or concatenated tracking data frame.}
  \item{workers}{Number of parallel workers. Defaults to detected cores minus one.}
  \item{seed}{Optional random seed for reproducible fly ID suffixes.}
  \item{var_for_slope}{Columns used for moving-window slope analysis.}
  \item{rsq_limit}{Minimum R-squared value retained in slope results.}
  \item{filter_negative_slope}{Set to 1 to remove negative slopes.}
  \item{filter_rsq_slope}{Set to 1 to apply the R-squared threshold.}
  \item{limit_late_flies}{Set to 1 to enable late-fly filtering.}
  \item{frame_include_limit}{Frame threshold used by late-fly filtering.}
  \item{quiet}{Suppress progress messages when TRUE.}
}
\value{
A named list containing duplicate-filtered detections, base tracks, Kalman-smoothed tracks, and slope results.
}
\description{
Runs duplicate filtering, frame-to-frame fly identification, trajectory interpolation and Kalman smoothing, and moving-window climbing velocity analysis.
}
