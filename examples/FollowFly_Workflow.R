library(FlyTrackPACK)

trackedDF <- read.csv("path/to/your_data.raw.csv")
results <- run_flytrack(trackedDF, seed = 1)

final_base_trackedDF <- results$base_tracks
final_Kalman_trackedDF <- results$kalman_tracks
final_slope_dataset <- results$slopes
