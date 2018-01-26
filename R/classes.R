#' The thresh class
#'
#' The function \code{nmathresh} returns S3 objects of class \code{thresh}.
#'
#' @rdname thresh-class
#' @name thresh-class
#'
#' @details Objects of class \code{thresh} have the following components:
#'   \describe{
#'   \item{\code{thresholds}}{A data frame with columns \code{lo} and \code{hi}
#'   for the lower and upper thresholds, and \code{lo.newkstar} and
#'   \code{hi.newkstar} for the new optimal (or rank-\code{trt.rank}) treatments
#'   at each of the thresholds.}
#'   \item{\code{U}}{The threshold solutions matrix. One column for each data
#'   point \eqn{m}, one row for each contrast \eqn{d_{ab}} (in ascending order).
#'   The elements \eqn{U_{ab,m}} describe the amount of adjustment to data point
#'   \eqn{y_m} required to reverse the relative ranking of treatments \eqn{a}
#'   and \eqn{b}. This matrix is particularly useful for deriving thresholds for
#'   more complex decisions (e.g. bias-adjustment thresholds for a new treatment
#'   entering the top two, for any change in rank of the top three, etc.)}
#'   \item{\code{Ukstar}}{The threshold solutions matrix limited to contrasts
#'   involving \eqn{k^*}. In other words, the rows of \code{U} corresponding to
#'   contrasts of the form \eqn{d_{ak^*}} or \eqn{d_{k^*a}}. Elements
#'   \eqn{U_{ak^*,m}} of this matrix describe the amount of adjustment to data
#'   point \eqn{y_m} required to make treatment \eqn{a} optimal (or
#'   rank-\code{trt.rank}) over \eqn{k^*}.}
#'   \item{\code{H}}{The influence matrix of the data on the basic treatment
#'   parameters. One column for each data point \eqn{m}, one row for each basic
#'   treatment parameter \eqn{d_k}. Elements \eqn{H_{k,m}} describe the
#'   influence of data point \eqn{y_m} on parameter \eqn{d_k}. This matrix can
#'   be used to derive more complex thresholds (e.g. 2D thresholds
#'   for simultaneous adjustments to two data points, or thresholds for common
#'   adjustments to a group of data points).}
#'   \item{\code{kstar}}{The base-case optimal (or rank-\code{trt.rank})
#'   treatment \eqn{k^*}.}
#'   \item{\code{call}}{A list containing all the arguments defined in the
#'   original call to \code{nma_thresh}.}
#'   }
#'
#' @seealso \code{\link{nma_thresh}}
#'
#' @exportClass thresh
#'
NULL

#' @export
print.thresh <- function(x, n = 6L, ...){

  cat("A thresh object. For help, see ?'thresh-class'.\n\n")

  if (x$call$trt.rank == 1) {
    cat("Base-case optimal treatment is k* = ", paste0(x$kstar, collapse = ", "), ".\n\n", sep = "")
  } else {
    cat("Base-case rank ", x$call$trt.rank, " treatment is k* = ",
        x$kstar, ".\n\n", sep = "")
  }

  print(utils::head(x$thresholds, n), ...)

  N <- nrow(x$thresholds)
  if (N > n & n > 0) cat("... ", N - n, " further rows omitted ...", sep = "")
  else if (n < 0 & -n < N) cat("... ", -n, " further rows omitted ...", sep = "")

}
