imputeCol <- function(mat_tmp, coord_tmp, col_to_fix, direction = "both") {
  
  all_y_list <- sort(unique(coord_tmp[coord_tmp$x == col_to_fix, ]$y))
  for (i in all_y_list) {
    print(i)
    if (direction == "both") {
      idx_to_average <- c(which(coord_tmp$x == col_to_fix + 1 & coord_tmp$y == i),
                          which(coord_tmp$x == col_to_fix - 1 & coord_tmp$y == i))
    }
    
    if (direction == "left") {
      idx_to_average <- c(
        which(coord_tmp$x == col_to_fix - 1& coord_tmp$y == i),
        which(coord_tmp$x == col_to_fix - 2 & coord_tmp$y == i))
    }
    
    if (direction == "right") {
      idx_to_average <- c(which(coord_tmp$x == col_to_fix + 1 & coord_tmp$y == i),
                          which(coord_tmp$x == col_to_fix + 2 & coord_tmp$y == i))
    }
    
    idx_to_fix <- which(coord_tmp$x == col_to_fix & coord_tmp$y == i)
    if (length(idx_to_average) > 0) {
      mat_tmp[, idx_to_fix] <- round(rowMeans(mat_tmp[, idx_to_average, drop = FALSE]))
    }
    
  }
  return(mat_tmp)
  
}

imputeRow <- function(mat_tmp, coord_tmp, row_to_fix, direction = "both") {
  
  all_x_list <- sort(unique(coord_tmp[coord_tmp$y == row_to_fix, ]$x))
  for (i in all_x_list) {
    print(i)
    
    if (direction == "both") {
      idx_to_average <- c(which(coord_tmp$y == row_to_fix + 1 & coord_tmp$x == i),
                          which(coord_tmp$y == row_to_fix - 1 & coord_tmp$x == i))
    }
    
    if (direction == "down") {
      idx_to_average <- c(which(coord_tmp$y == row_to_fix - 1 & coord_tmp$x == i),
                          which(coord_tmp$y == row_to_fix - 2 & coord_tmp$x == i))
    }
    
    if (direction == "up") {
      idx_to_average <- c(which(coord_tmp$y == row_to_fix + 1 & coord_tmp$x == i),
                          which(coord_tmp$y == row_to_fix + 2& coord_tmp$x == i))
    }
    
    idx_to_fix <- which(coord_tmp$y == row_to_fix & coord_tmp$x == i)
    if (length(idx_to_average) > 0) {
      mat_tmp[, idx_to_fix] <- round(rowMeans(mat_tmp[, idx_to_average, drop = FALSE]))
    }
    
  }
  return(mat_tmp)
  
}
