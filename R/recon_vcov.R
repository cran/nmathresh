#' Reconstruct likelihood covariance matrix
#'
#' Reconstruct the contrast-level likelihood covariance matrix from prior and
#' posterior covariance matrices. The resulting likelihood covariance matrix can
#' then be used to perform a contrast-level threshold analysis with the function
#' \code{nma_thresh}.
#'
#' @param post Posterior covariance matrix.
#' @param prior.prec Prior precision. Defaults to .0001 which is a common flat
#'   prior for NMA. Not used if \code{prior.vcov} is specified.
#' @param prior.vcov Prior covariance matrix. Defaults to a diagonal matrix of
#'   the same size as post, with elements \code{1/prior.prec}.
#' @param X Contrast design matrix. If omitted a complete network is assumed.
#' @param verbose Print intermediate matrices? Defaults to \code{FALSE}.
#'
#' @details Full details of the calculation are given by Phillippo \emph{et al.}
#'   (2017). Briefly, the aim is to recover the contrast-level
#'   likelihood covariance matrix \eqn{V} that would have led to the posterior
#'   covariance matrix \eqn{\Sigma} being obtained from a fixed effects NMA,
#'   with design matrix \eqn{X} and prior covariance matrix \eqn{\Sigma_d} for a
#'   normal prior on the basic treatment parameters. This is possible in this
#'   case via the equation (resulting from conjugacy):
#'
#'   \deqn{\Sigma^{-1} = X^TV^{-1}X + \Sigma_d^{-1}.}
#'
#'   When the treatment network is complete (i.e. fully connected), this
#'   equation may be rearranged exactly.
#'
#'   When the treatment network is incomplete (i.e. not all treatments are
#'   directly compared), this equation may be solved through the use of
#'   non-negative least squares (NNLS).
#'
#' @return A matrix; the reconstructed likelihood covariance matrix.
#' @import nnls
#' @seealso \code{\link{nma_thresh}}.
#' @export
#'
#' @examples
#' # Please see the vignette "Examples" for worked examples including use of
#' # this function, including more information on the brief code below.
#'
#' vignette("Examples", package = "nmathresh")
#'
#' ### Contrast level thresholds for Thrombolytic treatments NMA
#' K <- 6   # Number of treatments
#'
#' # Contrast design matrix is
#' X <- matrix(ncol = K-1, byrow = TRUE,
#'             c(1, 0, 0, 0, 0,
#'               0, 1, 0, 0, 0,
#'               0, 0, 1, 0, 0,
#'               0, 0, 0, 1, 0,
#'               0, -1, 1, 0, 0,
#'               0, -1, 0, 1, 0,
#'               0, -1, 0, 0, 1))
#'
#' # Reconstruct hypothetical likelihood covariance matrix using NNLS
#' lik.cov <- recon_vcov(Thrombo.post.cov, prior.prec = .0001, X = X)
#'

recon_vcov <- function(post, prior.prec=.0001,
                       prior.vcov=diag(1 / prior.prec, dim(post)[1]),
                       X=NULL, verbose=FALSE) {


  # Prior and posterior precision matrices
  tryCatch(prior.P <- solve(prior.vcov),
           error = function(err) {
             err$message <- paste0(
               "Cannot invert prior.vcov, should be a (square, positive definite) covariance matrix.\n",
               err$message)
             stop(err)
           })
  tryCatch(post.P <- solve(post),
           error = function(err) {
             err$message <- paste0(
               "Cannot invert post, should be a (square, positive definite) covariance matrix.\n",
               err$message)
             stop(err)
           })

  # X^T V^-1 X is equal to the difference of these precision matrices
  D <- post.P - prior.P

  #### COMPLETE NETWORK CASE, SOLVE EXACTLY ####
  if (is.null(X)) {
    # Match terms row by row to get likelihood precisions
    # Off-diagonal elements are simply switched sign of D elements
    M1 <- -D + diag(diag(D))

    # Diagonal elements are D diagonal minus sum of other row elements
    M2 <- diag(diag(D) - rowSums(M1))

    # M is the matrix of p elements
    M <- M1 + M2

    # Now read off the elements

    # Diagonal elements are the p_1b
    p1 <- diag(M)

    # Upper triangular rows are
    #    p_23 ......... p_2K
    #         p_34 .... p_3K
    #              .........
    #                   p_(K-1)K

    # R reads down columns, so use lower.tri instead (symmetric matrix) to get
    # the remaining parameters in the right order.
    p2 <- M[lower.tri(M)]

    # The likelihood precision matrix is then
    lik.cov <- diag(1 / c(p1, p2))
    return(lik.cov)
  }

  #### INCOMPLETE NETWORK CASE, FIND OPTIMUM SOLUTION WITH NNLS ####
  # Find NNLS solution to the equation Rp=s
  else{

    # Basic parameter check
    if (dim(post)[1] != dim(X)[2]) {
      stop("Number of parameters in design matrix does not match posterior covariance matrix.")
    }

    # Number of treatments
    K <- dim(post)[1] + 1

    # Construct s vector from difference in precision matrices
    s <- D[lower.tri(D, diag=TRUE)]

    # Construct R matrix from the design matrix X
    R <- matrix(0, nrow=K*(K-1)/2, ncol=K*(K-1)/2)

    # Loop through rows of X
    for (i in 1:nrow(X)){
      a.ind <- ifelse(any(X[i,] == -1), which(X[i,] == -1), 0)
      b.ind <- ifelse(any(X[i,] == 1), which(X[i,] == 1), 0)

      # Edges ab
      if (a.ind>0 && b.ind>0) {
        R[1 + K*(K-1)/2 - (K - a.ind + 1)*(K - a.ind)/2,
          1 + K*(K-1)/2 - (K - a.ind)*(K - a.ind - 1)/2 + b.ind - a.ind - 1] <- 1
        R[1 + K*(K-1)/2 - (K - b.ind + 1)*(K - b.ind)/2,
          1 + K*(K-1)/2 - (K - a.ind)*(K - a.ind - 1)/2 + b.ind - a.ind - 1] <- 1

        R[1 + K*(K-1)/2 - (K - a.ind + 1)*(K - a.ind)/2 + b.ind - a.ind,
          1 + K*(K-1)/2 - (K - a.ind)*(K - a.ind - 1)/2 + b.ind - a.ind - 1] <- -1

      }
      # Edges 1b
      else if (a.ind == 0 && b.ind > 0) {
        R[1 + K*(K-1)/2 - (K - b.ind + 1)*(K - b.ind)/2, b.ind] <- 1
      }
    }

    # Drop parameters with no edges (ie remove columns with all zeros)
    R <- R[, !apply(R, 2, identical, rep(0, K*(K-1)/2))]

    if (verbose) {
      cat("Design matrix R for NNLS system Rp=q is:\n")
      print(R)
    }

    # Solve using NNLS to obtain likelihood precisions
    sol <- nnls::nnls(R, s)
    message(
      "Likelihood precisions found using NNLS.\n",
      "Residual Sum of Squares: ", format(sol$deviance, digits=6), "\n",
      "             --------------------"
    )

    # RSS from elements which are fixed (i.e. trying to equate posterior precision with zero)
    ind.R0 <- ifelse(any(apply(R, 1, function(x) all(x == 0))),
                         which(apply(R, 1, function(x) all(x == 0))),
                         NA)
    rss0 <- ifelse(is.na(ind.R0), 0, sum((R %*% sol$x - s)[ind.R0]^2))
    message(
		  "              RSS fixed: ", format(rss0, digits=6)
		)

    # RSS from elements which are actually fitted
    rss1 <- ifelse(is.na(ind.R0), sum((R %*% sol$x - s)^2), sum((R %*% sol$x - s)[-ind.R0]^2))
    message(
      "             RSS fitted: ", format(rss1, digits=6), "\n",
      "             --------------------"
    )

    lik.cov <- diag(1 / sol$x)

    # Check whether any variances calculated are infinite
    inftf <- is.infinite(diag(lik.cov))
    if (any(inftf)) {
      warning("Returned some infinite variances. These will be not be included in KL calculation.")
    }

    # Kullback-Leibler Divergence of estimated posterior (arising from fitted
    # likelihood with infinite variances removed) from the "true" posterior

    # The posterior cov from the fitted likelihood cov
    post.cov.fit <- solve(solve(prior.vcov) + t(X[!inftf,]) %*%
                                              solve(lik.cov[!inftf,!inftf]) %*%
                                              X[!inftf,])

    # KL divergence of this from "true" post cov
    KL.div <- 1/2 * (sum(diag(solve(post.cov.fit) %*% post)) - dim(post)[1] +
                     log(det(post.cov.fit)/det(post)))
    message("Kullback-Leibler Divergence of fitted from 'true' posterior: ",
            KL.div)

    return(lik.cov)
  }
}
