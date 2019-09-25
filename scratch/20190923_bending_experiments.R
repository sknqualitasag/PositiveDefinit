#' ---
#' title: Bending For Non-Positive Symmetric Matrices
#' date:  "`r Sys.Date()`"
#' ---
#'
#'
#' ## Disclaimer
#' This script attempts to understand how different approaches for bending non-positive definite matrices work.
#' We start with an approach that comes from Schaeffer that is implemented in function `makPD()` in this github
#' repository.
#'
#'
#' ## A First Approach using `makPD()`
#' The function `makPD()` is unrolled and we try to explain what happens. The only parameter that is passed to `makPD()` is
#' the matrix `A`. Therefore, we directly assign `A` to a test matrix that will be used.
#+ assign-test-matrix
(A <- matrix(c(100,80,20,6,
               80,50,10,2,
               20,10,6,1,
                6,2,1,1), ncol = 4, byrow=TRUE))

#' The first step in function `makPD()` is to call the function `eigen()` on the input matrix `A`. The result consisting of
#' the eigenvalues and the eigenvectors is assigned to `D`.
#+ eigen-val-vec
(D <- eigen(A))

#' From this result, we can check whether any of the eigenvalues is negative.
#+ neg-eigen-val
any(D$values < 0)

#' After we found any negative eigenvalues, the question is how many of them are negative.
#+ nr-neg-ev
sum(D$values < 0)

#' Because it is only one negative eigenvalue, it must be the last one, because they are ordered. The index vector of the
#' negative eigenvalues is found by
#+ idx-vec-neg-ev
which(D$values < 0)

#'
#' ### Variable Init
#' The following variables are initialized
#+ var-init
sr = 0
nneg = 0
V = D$values
U = D$vectors
N = nrow(A)

#' ### Summing over negative eigenvalues
#' In the following loop the number of negative eigenvalues is counted and they are summed up
#+ loop-sum-ev
for(k in 1:N){
  if(V[k] < 0){
    nneg = nneg + 1
    sr = sr + V[k] + V[k]
  }
}
nneg
sr

#' This loop can completely be avoided by the following to statements. This should be verified with an example where more
#' than one eigenvalue is negative.
#+ rep-loop1
nneg <- sum(V < 0)
sr <- 2 * sum(V[which(V < 0)])
nneg
sr

#' The sum of the negative eigenvalues is squred, multiplied by 100 and added to one leading to
#+ neg-ev-wr
(wr <- sr * sr * 100 + 1)

#' At this point it is not clear what this means.
#'
#' ### Correction Loop
#' In a second loop, for each negative eigenvalue, the difference of `sr` minus the respective negative eigenvalue is  squared
#' and multiplied with the smallest positive eigenvalue. That product is divided by `wr` and the result is taken as the
#' corrected eigenvalue.
#+ correction-loop
p = V[N - nneg]
V_old <- V;V_old
for(m in 1:N){
  if(V[m] < 0){
    c = V[m]
    V[m] = p*(sr-c)*(sr-c)/wr
  }
}
V

#' The above loop can be vectorized without special functions, we just have to assign the vector with all negative eigenvalues
#' to a variable name. Because `c` is a reserved word in R, it is not a good choice for a variable name. Hence, we replace it
#' with `vec_neg_ev`.
#+ vectorzize-corr-loop
V <- V_old
vec_neg_ev <- V[which(V < 0)]
V[which(V < 0)] <- p * (sr - vec_neg_ev) * (sr - vec_neg_ev) / wr
V

#' ### Reconstruction of corrected matrix
#' The matrix is reconstructed using the corrected eigenvalues with the eigenvectors-matrix from the original input matrix
#+ corr-mat-recon
U %*% diag(V) %*% t(U)


#' ## Vectorized Version of `makPD()`
#' We take the same approach as implemented in `makPD()`, but we try to come up with a vectorized version
#+ cleanup-before-vect-makepd2
rm(list = ls())
#+ vect-makepd2
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


#' ## Testing vectorized Version `makePD2`
#' We start with the initial matrix `G`
#+ test-vect-makepd2
G<- matrix(c(100,80,20,6,
             80,50,10,2,
             20,10,6,1,
             6,2,1,1), ncol = 4, byrow=TRUE)
makePD2(G)

#' ## Comparison with `makPD()`
#' As a reference, we include `makPD()` in here and the two results are compared.
#+ compare-makpd
makPD = function(A){
  D = eigen(A)
  sr = 0
  nneg = 0
  V = D$values
  U = D$vectors
  N = nrow(A)
  for(k in 1:N){
    if(V[k] < 0){
      nneg = nneg + 1
      sr = sr + V[k] + V[k]
    }
  }
  wr = (sr*sr*100+1)
  p = V[N - nneg]
  for(m in 1:N){
    if(V[m] < 0){
      c = V[m]
      V[m] = p*(sr-c)*(sr-c)/wr
    }
  }
  A = U %*% diag(V) %*% t(U)
  return(A)
}

### # Apply to G
makPD(G)

#' Difference between two results
#+ diff-results
makePD2(G) - makPD(G)

#' ## More Than One Negative Eigenvalue
#' Espicially the vectorized function `makePD2()` must be tested for a matrix that has more than one negative eigenvalue.
#' We start by constructing such a matrix using
#+ construct-test-mat-g2
eG <- eigen(G)
vec_eG_values <- eG$values
mat_eG_vectors <- eG$vectors
vec_eG_values;mat_eG_vectors
# determine number of negative eigenvalues
n_nr_neg <- sum(vec_eG_values < 0)
n_idx_last_post_ev <- length(vec_eG_values) - n_nr_neg
if (n_idx_last_post_ev < 1) stop(" * ERROR: Matrix does not have any positive eigenvalues")
# make that ev also negative
vec_eG_values[n_idx_last_post_ev] <- -vec_eG_values[n_idx_last_post_ev]
vec_eG_values
# construct the test matrix
G2 <- mat_eG_vectors %*% diag(vec_eG_values) %*% t(mat_eG_vectors)
G2

#' Testing the two functions leads to
#+ test-mat-g2
makPD(G2)
makePD2(G2)
makPD(G2) - makePD2(G2)

#' ## Parametrisation
#' The current implementation bends the given matrix in a fixed fashion. All it does is to increment the negative
#' eigenvalues between 0 and the smallest positive eigenvalue. If the smallest positive eigenvalue is very small, then
#' the corrected eigenvalues will be about 100 times smaller which corresponds to the hard-coded factor in the computation
#' of the weight `wr`.
#'
#' For some applications it might be interesting to be able to specify a ratio between the largest and the smallest
#' eigenvalue. With that it should be possible to scale all eigenvalues to be within that ratio boundary.
#'
#' ## Ratio of Eigenvalues
#' Instead of just changing the negative eigenvalues by mapping them between zero and the smallest eigenvalue,
#' we would like to correct eigenvalues such that all of them are poitive and such that the ratio between the
#' largest and the smallest eigenvalue is below a certain threshold.
#'
#' In principle this approach can be implemented similarly to the one that corrects the negative eigenvalues,
#' presented so far. The only adaptation that one must do is to replace the limit of the eigenvalues that must
#' be changed from 0 to that number determined by the maximum ration. Everything else should stay the same,
#' modulo a special treatment of any negative eigenvalues. An easy solution to that might be to do it in two
#' steps, first make all eigenvalues positive using an approach used in `makePD2()`. Then we can use a second
#' step to come up with a correction that ensures a maximum ration between eigenvalues. On the other hand it
#' should be possible to provide a one step solution.
#'
#' The function is called `make_pd_rat_ev()`. As arguments the function takes the input matrix and a maximum
#' ratio between largest and smallest eigenvalue.
#+ def-make_pd_rat_ev
make_pd_rat_ev <- function(A, pn_max_ratio){
  # get eigenvalue/eigenvector decomposition
  D <- eigen(A)
  # assign eigenvectors and eigen values to separate variables
  vec_eval <- D$values
  mat_evec <- D$vectors
  # number of negative eigenvalues
  nneg <- sum(vec_eval < 0)
  # correction based on ratio of max and minimum of the absolute eigenvalues
  max_ev <- max(vec_eval)
  n_abs_rat <- max_ev / min(abs(vec_eval))
  # if the ratio is ok and no negative eigenvalues, then we can return the input matrix
  if (n_abs_rat < pn_max_ratio && nneg < 1)
    return(A)
  # in case we have to correct them, we correct all eigenvalues that are negative or
  # outside the ratio boundary. Find eigenvalue that must be corrected due to the ratio boundary
  vec_idx_ev_out <- which(max_ev/vec_eval > pn_max_ratio)
  # add those which are negative
  if (nneg > 0) vec_idx_ev_out <- unique(c(vec_idx_ev_out, which(vec_eval < 0)))
  # the boundary for the smallest ev is
  n_min_ev <- max_ev / pn_max_ratio
  # the range between the last eigenvalue that does not need correction (n_last_ev_not_corrected)
  # and the smallest eigenvalue of the input matrix is projected to the range between the
  # n_last_ev_not_corrected and the minimum eigenvalue after correction based on max desired ratio
  n_last_ev_not_corrected <- vec_eval[vec_idx_ev_out[1]-1]
  n_out_range <- max(n_last_ev_not_corrected - vec_eval[vec_idx_ev_out])
  # the corrected range between the last eigenvalue inside the range and the ratio boundary is
  n_corr_range <- n_last_ev_not_corrected - n_min_ev
  n_range_rat <- n_corr_range / n_out_range
  # correct according to the range ratios
  vec_eval[vec_idx_ev_out] <- n_last_ev_not_corrected - (n_last_ev_not_corrected - vec_eval[vec_idx_ev_out]) * n_range_rat
  # return reconstructed matrix
  return(mat_evec %*% diag(vec_eval) %*% t(mat_evec))

}

# test
(G_rat_100 <- make_pd_rat_ev(G, pn_max_ratio = 100))
eigen(G_rat_100)


(G_rat_200 <- make_pd_rat_ev(G, pn_max_ratio = 200))
eigen(G_rat_200)



