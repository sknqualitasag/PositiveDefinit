# #' Bending of Matrix A
# #'
# #' @description
# #' The input matrix A is decomposed into its eigen-values and eigen-vectors,
# #' The negative eigen-values are projected into the range between zero and
# #' the smallest positive eigen-value.
# #'
# #' @param A input matrix
# #'
# #' @return Bended positive-definite matrix A
# #' @export makePD2
# #'
# #' @examples
# #' G<- matrix(c(100,80,20,6,80,50,10,2,20,10,6,1,6,2,1,1), ncol = 4, byrow=TRUE)
# #' makePD2(G)
makePD2 <- function(A){
  # compute eigenvalue-eigenvector decomposition
  D  <-  eigen(A)
  # assign separate values
  V <- D$values
  U <- D$vectors
  N <- nrow(A)
  # determine number of negative eigenvalues and sum twice the negative ev
  nneg <- sum(V < 0)
  vec_neg_ev <- V[which(V < 0)]
  sr <- 2 * sum(vec_neg_ev)
  # compute the weight factor wr and determine the smallest positive ev
  wr = (sr*sr*100+1)
  p = V[N - nneg]
  # correct the negative eigenvalues
  V[which(V < 0)] <- p * (sr - vec_neg_ev) * (sr - vec_neg_ev) / wr
  # reconstruct A from eval-evec-decomposition and return
  return(U %*% diag(V) %*% t(U))
}
