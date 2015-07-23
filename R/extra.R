#replicates matrix rows
replicate_rows <- function(mat, times)
{
  mat[rep(1:nrow(mat), times=times), ]
}
