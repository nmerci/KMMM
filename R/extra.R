#replicates matrix rows
replicate_rows <- function(mat, times)
{
  mat[rep(1:nrow(mat), times=times), ]
}

#calculates frequencies of each factor in vector
get_frequencies <- function(vec)
{
  as.data.frame(table(vec))$Freq
}

#forms an array (1,1,..,1, 2,2,..2, .., k,k,..,k)
replicate_numbers <- function(k, n)
{
  as.vector(sapply(1:k, rep, n))
}

#swaps choosen elements in two columns
swap_columns <- function(data, col1, col2, rows)
{
  data[rows, c(col1, col2)] <- data[rows, c(col2, col1)]
  data
}