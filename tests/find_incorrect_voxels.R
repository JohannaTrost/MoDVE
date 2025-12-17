find_incorrect_voxels <- function(seg_start, seg_end) {

  intersctd_voxels <- list()

  UniqueX <- c(ceiling(seg_start[1]), ceiling(seg_end[1]))
  UniqueY <- c(ceiling(seg_start[2]), ceiling(seg_end[2]))
  UniqueZ <- c(ceiling(seg_start[3]), ceiling(seg_end[3]))

  numX <- sum((UniqueX[2] - UniqueX[1]) > 0, na.rm=TRUE) + 1
  numY <- sum((UniqueY[2] - UniqueY[1]) > 0, na.rm=TRUE) + 1
  numZ <- sum((UniqueZ[2] - UniqueZ[1]) > 0, na.rm=TRUE) + 1

  for (x in seq_len(numX)) {
    xid <- UniqueX[x]

    for (y in seq_len(numY)) {
      yid <- UniqueY[y]

      for (z in seq_len(numZ)) {
        zid <- UniqueZ[z]
        intersctd_voxels <- c(intersctd_voxels, list(c(xid, yid, zid)))
      }
    }
  }
  return(intersctd_voxels)
}
