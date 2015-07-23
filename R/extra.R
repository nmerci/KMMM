#replicates matrix rows
replicate_rows <- function(mat, times)
{
  mat[rep(1:nrow(mat), times=times), ]
}

#calculate frequencies of each factor in vector
get_frequencies <- function(vec)
{
  as.data.frame(table(vec))$Freq
}

#forms an array (1,1,..,1, 2,2,..2, .., n,n,..n)
replicate_numbers <- function(k, n)
{
  as.vector(sapply(1:k, rep, n))
}