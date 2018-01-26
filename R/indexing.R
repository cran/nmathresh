#' Convert contrast indexing
#'
#' Functions for converting between \eqn{d_{ab}} indexing of contrasts (useful
#' notationally) and \code{d[i]} indexing used by R.
#'
#' @param a Vector of treatment codes \eqn{a}.
#' @param b Vector of treatment codes \eqn{b}.
#' @param i Vector of indices \code{i}.
#' @param K Total number of treatments.
#'
#' @return \code{d_ab2i} returns a vector of indices \code{i}. \code{d_i2ab}
#'   returns a data frame of indices \code{a} and \code{b}.
#' @export
#' @aliases d_i2ab
#'
#' @examples
#' d_ab2i(c(1,1,1, 2,2, 3), c(2,3,4, 3,4, 4), K=4)
#' d_i2ab(1:6, K=4)
#'
#' @note By convention, \eqn{1 \le a < b \le K}. If this is not the case, an
#'   error will be thrown. For a given number of treatments \eqn{K}, the total
#'   number of possible contrasts \eqn{d_{ab}} is \eqn{K(K-1)/2}, and hence
#'   \eqn{i \le K}. Again, if this is not the case, an error will be thrown.
#'

#' @describeIn d_ab2i Convert \code{d[i]} type indices to \eqn{d_{ab}} type indices.
#'
d_ab2i <- function(a, b, K) {
  stopifnot(b>a, a<K, b<=K, a>=1, b>=2, length(a)==length(b), length(K)==1, K==trunc(K))

  #sum(1:(K-1)) - sum(1:(K-a[i])) + (b-a)
  sum(1:(K-1)) - sapply((sapply(K-a, seq, from=1)), sum) + (b-a)
}

#' @describeIn d_ab2i Convert \eqn{d_{ab}} type indices to \code{d[i]} type indices.
#' @export
#'
d_i2ab <- function(i, K) {
  if (any(i > K*(K-1)/2)) {
    stop("Index i is greater than the total number of contrasts K(K-1)/2")
  }
  if (length(K) != 1 || K != trunc(K)) stop("K should be a single integer.")

  # Generate "triangular" numbers
  tri <- lower.tri(matrix(nrow=K-1, ncol=K-1), diag=T) %*% (K-1):1

  # Get a,b
  a <- b <- rep(NA, length(i))
  for (j in 1:length(i)) {
    a[j] <- which(i[j] <= tri)[1]
    b[j] <- i[j] - c(0, tri)[a[j]] + a[j]
  }

  return(data.frame(a=a, b=b))
}
