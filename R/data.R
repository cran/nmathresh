#' Posterior summary from Social Anxiety NMA
#'
#' A \code{summary.mcmc} object of the type produced by the \code{coda} package,
#' containing the requisite posterior summary information on the variables
#' \code{d} (basic treatment effect parameters), \code{delta} (shrunken random
#' effects estimates for each study), and \code{diff} (contrasts of treatment
#' effect parameters).
#'
#' @format A \code{summary.mcmc} object. The key components for our use are:
#'   \describe{ \item{statistics}{Matrix containing the posterior summary
#'   statistics of the variables \code{d}, \code{delta}, and \code{diff}, with
#'   columns for \code{Mean}, \code{SD}, \code{Naive SE}, and \code{Time-series
#'   SE} (also known as the Monte-Carlo standard error)} \item{quantiles}{Matrix
#'   containing the posterior 2.5\%, 25\%, 50\%, 75\%, and 97.5\% quantiles of
#'   the variables \code{d}, \code{delta}, and \code{diff}}}
#'
#' @seealso \code{\link[coda]{summary.mcmc}}, \code{\link{SocAnx.post.cov}}
#' @source Generated from WinBUGS output, using the WinBUGS code from
#'   Mayo-Wilson et al. (2014). See also \code{vignette("Examples", package =
#'   "nmathresh")}.
#'
#' @references Mayo-Wilson E, Dias S, Mavranezouli I, Kew K, Clark DM, Ades AE,
#'   et al. Psychological and pharmacological interventions for social anxiety
#'   disorder in adults: a systematic review and network meta-analysis. Lancet
#'   Psychiatry 2014;1:368-76.
#'   http://dx.doi.org/10.1016/S2215-0366(14)70329-3
#'

"SocAnx.post.summary"

#' Posterior covariance matrix from Social Anxiety NMA
#'
#' The posterior covariance matrix of the variables \code{d} (basic treatment
#' effect parameters) and \code{delta} (shrunken random effects estimates for
#' each study).
#'
#' @seealso \code{\link{SocAnx.post.summary}}
#' @source Generated from WinBUGS output, using the WinBUGS code from
#'   Mayo-Wilson et al. (2014). See also \code{vignette("Examples", package =
#'   "nmathresh")}.
#'
#' @references Mayo-Wilson E, Dias S, Mavranezouli I, Kew K, Clark DM, Ades AE,
#'   et al. Psychological and pharmacological interventions for social anxiety
#'   disorder in adults: a systematic review and network meta-analysis. Lancet
#'   Psychiatry 2014;1:368-76.
#'   http://dx.doi.org/10.1016/S2215-0366(14)70329-3
#'

"SocAnx.post.cov"

#' Posterior summary from Thrombolytics NMA
#'
#' A \code{summary.mcmc} object of the type produced by the \code{coda} package,
#' containing the requisite posterior summary information on the variables
#' \code{dd}, the contrasts of the treatment effect parameters.
#'
#' @format A \code{summary.mcmc} object. The key components for our use are:
#'   \describe{ \item{statistics}{Matrix containing the posterior summary
#'   statistics, with columns for \code{Mean}, \code{SD}, \code{Naive SE}, and
#'   \code{Time-series SE} (also known as the Monte-Carlo standard error)}
#'   \item{quantiles}{Matrix containing the posterior 2.5\%, 25\%, 50\%, 75\%,
#'   and 97.5\% quantiles}}
#'
#' @seealso \code{\link[coda]{summary.mcmc}}, \code{\link{Thrombo.post.cov}}
#' @source Generated from WinBUGS output, using the WinBUGS code from Caldwell
#'   et al. (2005). See also \code{vignette("Examples", package = "nmathresh")}.
#'
#' @references Caldwell DM, Ades AE, Higgins JPT. Simultaneous comparison of
#'   multiple treatments: combining direct and indirect evidence. Brit Med J
#'   2005;331:897-900. http://dx.doi.org/10.1136/bmj.331.7521.897
#'

"Thrombo.post.summary"

#' Posterior covariance matrix from Thrombolytics NMA
#'
#' The posterior covariance matrix of the basic treatment effect parameters.
#'
#' @seealso \code{\link{Thrombo.post.summary}}
#' @source Generated from WinBUGS output, using the WinBUGS code from Caldwell
#'   et al. (2005). See also \code{vignette("Examples", package = "nmathresh")}.
#'
#' @references Caldwell DM, Ades AE, Higgins JPT. Simultaneous comparison of
#'   multiple treatments: combining direct and indirect evidence. Brit Med J
#'   2005;331:897-900. http://dx.doi.org/10.1136/bmj.331.7521.897
#'

"Thrombo.post.cov"