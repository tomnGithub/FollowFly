duplicateFilter_Frame <- function(data = data.frame()) {
  for (frame_level in unique(data$frame)) {
    frame_data <-
      data[data$frame == frame_level,]

    for (i in 1:nrow(frame_data)) {
      current_row <- frame_data[i,]

      distances <-
        sqrt((current_row$x - frame_data$x) ^ 2 + (current_row$y - frame_data$y) ^
               2)

      distances[i] <- Inf

      nearest_index <- which.min(distances)

      nearest_neighbor <- row.names(frame_data[nearest_index, ])
      distance_to_NN <- min(distances)
      frame_data[i, "nearest_neighbor"] <- nearest_neighbor
      frame_data[i, "distanceNN"] <- distance_to_NN
      data_NN <- frame_data[nearest_neighbor,]
      current_id <- row.names(current_row)

      if (data_NN$isDuplicate < 1 &&
          !is.na(data_NN$nearest_neighbor)) {
        if (data_NN$nearest_neighbor == current_id &&
            distance_to_NN < 15) {
          frame_data[current_id, "isDuplicate"] <- 1
          frame_data[nearest_neighbor, "isDuplicate"] <- 2

          if (which.max(c(current_row$y, data_NN$y)) == 2) {
            data[current_id, "y"] <- data_NN$y
            data[current_id, "x"] <- data_NN$x
          }
          data[nearest_neighbor, "isDuplicate"] <- 2

        }
      }

    }

  }

  filtered_data <- data %>% filter(isDuplicate == 0)

  return(filtered_data)
}

prepData_flytrack <- function(data = data.frame()) {
  data$name <- NA
  data$vial <- data$vialID
  data$distance <- 0

  data$xDiff <- 0
  data$yDiff <- 0
  data$predX <- NA
  data$predY <- NA
  data$isPred <- 0
  data$distance_ac <- NA
  data$SRC_ac <- NA
  data$named_this_frame <- 0
  data$reorder <- 0
  data$reorder <- as.character(data$reorder)

  data$concatenated_names <- NA
  data$nearest_neighbor <- NA
  data$distanceNN <- NA
  return(data)
}
