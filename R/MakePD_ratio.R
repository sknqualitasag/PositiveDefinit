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
