#' Calculate thresholds and invariant intervals
#'
#' This function calculates decision-invariant bias-adjustment thresholds and
#' intervals for Bayesian network meta-analysis, as described by Phillippo
#' \emph{et al.} (2018). Thresholds are derived from the joint
#' posterior, and reflect the amount of change to a data point before the
#' treatment decision changes. Calculation is achieved using fast matrix
#' operations.
#'
#' @param mean.dk Posterior means of basic treatment parameters \eqn{d_k}.
#' @param lhood Likelihood (data) covariance matrix.
#' @param post Posterior covariance matrix (see details).
#' @param nmatype Character string, giving the type of NMA performed. One of
#'   "fixed" (fixed effects, the default) or "random" (random effects). May be
#'   abbreviated.
#' @param X [FE models only] Design matrix for basic treatment parameters.
#' @param mu.design [RE models only] Design matrix for any extra parameters.
#'   Defaults to \code{NULL} (no extra parameters).
#' @param delta.design [RE models only] Design matrix for delta, defaults to the
#'   \eqn{N \times N}{N x N} identity matrix.
#' @param opt.max Should the optimal decision be the maximal treatment effect
#'   (\code{TRUE}, default) or the minimum (\code{FALSE}).
#' @param trt.rank Rank of the treatment to derive thresholds for. Defaults to
#'   1, thresholds for the optimum treatment.
#' @param trt.code Treatment codings of the reference treatment and in the
#'   parameter vector \eqn{d_k}. Use if treatments re-labelled or re-ordered.
#'   Default is equivalent to \code{1:K}.
#' @param trt.sub Only look at thresholds in this subset of treatments in
#'   \code{trt.code}, e.g. if some are excluded from the ranking. Default is
#'   equivalent to \code{1:K}.
#' @param mcid Minimal clinically important difference for the decision (when
#'   \code{mcid.type = 'decision'}) or for changing the decision (when
#'   \code{mcid.type = 'change'}). Defaults to 0, use the maximal efficacy
#'   decision rule.
#' @param mcid.type Default \code{'decision'}, the decision rule is based on
#'   MCID (see details). Otherwise \code{'change'}, use the maximum efficacy
#'   rule, but only consider changing the decision when the alternative
#'   treatment becomes more effective than the base case by \code{mcid} or more.
#'
#' @details This function provides bias-adjustment threshold analysis for both
#'   fixed and random effects NMA models, as described by Phillippo \emph{et
#'   al.} (2018). Parameters \code{mean.dk}, \code{lhood}, and
#'   \code{post} are always required, however there are differences in the
#'   specification of \code{post} and other required parameters and between the
#'   fixed and random effects cases:
#'
#'   \describe{ \item{\strong{FE models}}{The design matrix \code{X} for basic
#'   treatment parameters is required. The posterior covariance matrix specified
#'   in \code{post} should only refer to the basic treatment parameters.}
#'   \item{\strong{RE models}}{The design matrix \code{mu.design} for additional
#'   parameters (e.g. covariates) is required, as is the design matrix
#'   \code{delta.design} for random effects terms. The posterior covariance
#'   matrix specified in \code{post} should refer to the basic treatment
#'   parameters, RE terms, and additional parameters \emph{in that order}; i.e.
#'   \code{post} should be the posterior covariance matrix of the vector
#'   \eqn{(d^T, \delta^T, \mu^T)^T}.} }
#'
#' @section Model details: \strong{The FE NMA model}
#'
#'   The fixed effects NMA model is assumed to be of the form \describe{
#'   \item{Prior:}{\eqn{d \sim \mathrm{N} ( d_0, \Sigma_d )}{d ~ N(d_0,
#'   \Sigma_d)}}
#'   \item{Likelihood:}{\eqn{y|d \sim \mathrm{N} ( \delta, V )}{y|d ~ N(\delta,
#'   V)}} \item{FE model:}{\eqn{\delta = Xd + M\mu}} }
#'
#'   The additional parameters \eqn{\mu} may be given any sensible prior; they
#'   do not affect the threshold analysis in any way.
#'
#'   \strong{The RE NMA model}
#'
#'   The random effects NMA model is assumed to be of the form \describe{
#'   \item{Priors:}{\eqn{ d \sim \mathrm{N} ( d_0, \Sigma_d ), \quad \mu \sim
#'   \mathrm{N} ( \mu_0, \Sigma_\mu )}{d ~ N(d_0, \Sigma_d), \mu ~ N(\mu_0,
#'   \Sigma_\mu)}}
#'   \item{Likelihood:}{\eqn{y|d,\mu,\tau^2 \sim \mathrm{N} ( L\delta + M\mu, V
#'   )}{y|d,\mu,\tau^2 ~ N(L\delta + M\mu, V)}}
#'   \item{RE model:}{\eqn{\delta \sim \mathrm{N} ( Xd, \tau^2 )}{\delta ~ N(Xd,
#'   \tau^2)}} }
#'
#'   The between-study heterogeneity parameter \eqn{\tau^2} may be given any
#'   sensible prior. In the RE case, the threshold derivations make the
#'   approximation that \eqn{\tau^2} is fixed and known.
#'
#' @section Decision rules:
#'
#'   The default decision rule is maximal efficacy; the optimal treatment is
#'   \eqn{ k^* = \mathrm{argmax}_k \mathrm{E}(d_{k})}{k* = argmax(E(d_k))}.
#'
#'   When \eqn{\epsilon} = \code{mcid} is greater than zero and
#'   \code{mcid.type = 'decision'}, the decision rule is no longer for a single
#'   best treatment, but is based on minimal clinically important difference. A
#'   treatment is in the optimal set if \eqn{\mathrm{E}(d_k) \ge
#'   \epsilon}{E(d_k) \ge \epsilon} and \eqn{\max_a \mathrm{E}(d_a) -
#'   \mathrm{E}(d_k) \le \epsilon}{max E(d_a) - E(d_k) \le \epsilon}.
#'
#'   When \code{mcid.type = 'change'}, the maximal efficacy rule is used, but
#'   thresholds are found for when a new treatment is better than the base-case
#'   optimal by at least \code{mcid}.
#'
#' @return An object of class \code{thresh}.
#' @seealso \code{\link{recon_vcov}}, \code{\link{thresh_forest}},
#'   \code{\link{thresh-class}}.
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
#' # Thresholds are then
#' thresh <- nma_thresh(mean.dk = Thrombo.post.summary$statistics[1:(K-1), "Mean"],
#'                      lhood = lik.cov,
#'                      post = Thrombo.post.cov,
#'                      nmatype = "fixed",
#'                      X = X,
#'                      opt.max = FALSE)
#'

nma_thresh <- function(mean.dk, lhood, post,
                       nmatype="fixed",
                       X=NULL,
                       mu.design=NULL, delta.design=NULL,
                       opt.max=TRUE, trt.rank=1, trt.code=NULL, trt.sub=NULL,
                       mcid=0, mcid.type='decision') {


## Basic parameter checks --------------------------------------------------

  # Fixed or random effects
  tryCatch(isFE <- match.arg(tolower(nmatype), c("fixed","random")) == "fixed",
           error = function(err) {
             stop('nmatype should be one of "fixed" or "random"')
           })

  # Warn if mu.design or delta.design given in FE case, or X given in RE case
  if (isFE && (!is.null(mu.design) || !is.null(delta.design))) {
    warning('nmatype = "fixed", so arguments mu.design and delta.design are ignored.')
  }
  if (!isFE && !is.null(X)) {
    warning('nmatype = "random", so argument X is ignored.')
  }

  # Get number of data points N
  if (dim(lhood)[1] == dim(lhood)[2]) {
    N <- dim(lhood)[1]
    message("Likelihood for N = ", N, " data points.")
  } else stop("Likelihood covariance matrix lhood should be square.")

  # Get number of treatments
  K <- length(mean.dk) + 1
  message("Number of treatments is K = ", K, ".")

  # Check X matrix for FE model
  if (isFE) {
    if (is.null(X)) stop("Design matrix X must be provided for FE models.")
    else if (dim(X)[1] != N || dim(X)[2] != K-1) {
      stop("Design matrix X should be N x (K-1).")
    }
  }

  # Get number of extra parameters
  if (isFE || is.null(mu.design)) {
    m <- 0
  } else if (nrow(mu.design) != N) {
    stop("Number of rows in mu.design does not equal N.")
  } else {
    m <- ifelse(is.null(mu.design),0,ncol(mu.design))
    message("Number of extra parameters in m.design is m = ",m,".")
  }

  # Get number of delta parameters
  if (isFE) {
    n.delta <- 0
  } else if (is.null(delta.design)) {
    delta.design <- diag(nrow = N)
  } else if (nrow(delta.design) != N) {
    stop("Number of rows in delta.design does not equal N.")
  }
  if (!isFE) {
    n.delta <- ncol(delta.design)
    message("Number of delta parameters is n.delta = ",n.delta,".")
  }


  # Check posterior covariance matrix is n.delta+m+K-1 square
  if (nrow(post) != ncol(post)) {
    stop("Posterior covariance matrix should be square.")
  } else if (nrow(post) != n.delta + m + K-1) {
    if (isFE) stop("Posterior covariance matrix should be K-1 square.")
    else stop("Posterior covariance matrix should be n.delta+m+K-1 square.")
  }

  # Treatment rank
  if (length(trt.rank) > 1 | trt.rank != round(trt.rank)) {
    stop("trt.rank should be a single integer.")
  } else if (trt.rank < 1 | trt.rank > K) {
    stop("trt.rank should be between 1 and K (number of trts).")
  }

  # Note about recoded treatments
  if (is.null(trt.code)) {
    trt.code <- 1:K
  }
  else if (length(trt.code) != K) stop("trt.code should be of length K.")
  else {
    message("Using recoded treatments. Reference treatment is ", trt.code[1],
        ". Parameter vector is:\n",
        "\t", paste0("d[", trt.code[-1], "]", collapse=", ")
        )
  }

  # Treatment subset
  if (is.null(trt.sub)){
    trt.sub <- trt.code
  } else if (length(trt.sub)>K) stop("Length of trt.sub should be <= K.")
  else {
    message("Deriving thresholds on a subset of treatments:")
    message("\t", paste(trt.sub, collapse=", "))
  }

  trt.sub.internal <- which(trt.code %in% trt.sub)

  # Error if trt.rank > length(trt.sub)
  if (trt.rank > length(trt.sub)) {
    stop("trt.rank is larger than the length of trt.sub")
  }

  # mcid should be a single non-negative numeric value
  if (!is.numeric(mcid) | length(mcid) != 1 | mcid < 0) {
    stop("mcid should be a single non-negative numeric value")
  }

  # Check mcid.type
  mcid.type <- match.arg(mcid.type, c('decision', 'change'))

  # Can't use mcid decision rule and treatment ranks
  if (mcid > 0 & mcid.type == 'decision' & trt.rank > 1) {
    stop("Can't use mcid decision rule and trt.rank at the same time.")
  }


## Pre-processing ----------------------------------------------------------

  # Get contrast details
  d_ab <- d_i2ab(1:(K*(K-1)/2), K)

  # Create contrast "design" matrix
  D <- matrix(0, nrow=K*(K-1)/2, ncol=K-1)
  D[cbind(1:nrow(D), d_ab$a - 1)] <- -1
  D[cbind(1:nrow(D), d_ab$b - 1)] <- 1

  # Create vector of contrasts d_ab
  contr <- as.vector(D %*% mean.dk)


## Derive influence matrix -------------------------------------------------

  ## FE models

  # If the likelihood covariance matrix was reconstructed to use in a
  # contrast-level analysis, e.g. using recon_vcov, there might be infinite
  # variances. We'll check that lhood is diagonal, and then handle these
  # correctly.

  if (isFE) {

    if (!Matrix::isDiagonal(lhood)) {
      inflmat <- post %*% crossprod(X, solve(lhood))
    } else {
      inflmat <- post %*% crossprod(X, diag(1/diag(lhood)))
    }

  } else {

  ## RE models
    Bstar <- post[1:(K-1), K:(K-1+n.delta)] # Posterior cov submatrix B*
    if (m > 0) {
      Cstar <- post[1:(K-1), (K+n.delta):(K-1+n.delta+m)] # Posterior cov submatrix C*
      inflmat <-
        (tcrossprod(Bstar, delta.design) + tcrossprod(Cstar, mu.design)) %*% solve(lhood)
    } else {
      inflmat <- Bstar %*% crossprod(delta.design, solve(lhood))
    }

  }

  # Add row names to H matrix (inflmat)
  rownames(inflmat) <- paste0("d[", trt.code[2:K], "]")


## Derive solution matrix U -------------------------------------------------

  if (mcid > 0 & mcid.type == 'change') {
    threshmat <- sweep(1 / (D %*% inflmat), 1, -contr - sign(contr)*mcid, "*")

    ## -- For mcid.type = "change" --
    # For mcid > 0, if a contrast is negative, we want a new decision when the
    # contrast is > +mcid. If a contrast is positive, we want a new decision
    # when the contrast is < -mcid. In other words, the contrast has to be
    # overturned by an extra mcid.new amount.
  } else {
    threshmat <- sweep(1 / (D %*% inflmat), 1, -contr + sign(contr)*mcid, "*")

    ## -- For mcid.type = "decision" --
    # And also for standard maximal efficacy rule, when mcid = 0 anyway.
    # For mcid > 0, if a contrast is negative, we want to know when the
    # contrast is = -mcid. If a contrast is positive, we want to know when
    # the contrast is = +mcid
  }

  # Add row names to U matrix (threshmat)
  rownames(threshmat) <- paste0("d[", trt.code[d_ab$a], ",", trt.code[d_ab$b], "]")

  # Now we only need to look at contrasts involving the optimal treatment k*
  # Updated to handle trt.rank, to pick out other ranked treatments than the
  # optimal treatment k* in first place.
  # Updated to handle trt.sub, only look for k* in a subset of treatments.

  mean.dk.subNA <- mean.dk
  mean.dk.subNA[!(1:(K - 1) %in% (trt.sub.internal - 1))] <- NA

  if (opt.max) {
    if (mcid > 0 & mcid.type == 'decision') {
      kstar <- which(mean.dk.subNA >= mcid & max(mean.dk.subNA, na.rm = TRUE) - mean.dk.subNA <= mcid) + 1
    } else {
      kstar <- order(c(0, mean.dk.subNA), decreasing = TRUE)[trt.rank]
    }
  } else if (!opt.max) {
    if (mcid > 0 & mcid.type == 'decision') {
      kstar <- which(mean.dk.subNA <= -mcid & min(mean.dk.subNA, na.rm = TRUE) - mean.dk.subNA >= -mcid) + 1
    } else {
      kstar <- order(c(0, mean.dk.subNA), decreasing = FALSE)[trt.rank]
    }
  }

  if (mcid > 0 & mcid.type == 'decision') {
    message("Current optimal treatment set is k* = ", paste(trt.code[kstar], collapse = ", "), ".")
  } else if (trt.rank == 1) {
    message("Current optimal treatment is k* = ", trt.code[kstar], ".")
  } else {
    message("Current rank ", trt.rank, " treatment is k = ", trt.code[kstar], ".")
  }

  # Which rows of U correspond to treatments in kstar?
  # We also only want contrasts involving treatments in trt.sub.
  if (mcid > 0 & mcid.type == 'decision') {
    # If using MCID decision rule, we do want thresholds for contrasts between treatments in kstar
    contr.kstar <- which((d_ab$a %in% kstar & d_ab$b %in% trt.sub.internal) |
                           (d_ab$b %in% kstar & d_ab$a %in% trt.sub.internal))
  } else {
    # Otherwise we don't care about switches within kstar, so only one of a or b can be in kstar (xor).
    contr.kstar <- which(xor(d_ab$a %in% kstar, d_ab$b %in% kstar) &
                           d_ab$a %in% trt.sub.internal &
                           d_ab$b %in% trt.sub.internal)
  }

  # So we look in the corresponding rows of the threshold matrix
  threshmat.kstar <- threshmat[contr.kstar, , drop = FALSE]


## Derive thresholds -------------------------------------------------------

  # Split threshmat and inflmat into a list of columns to mapply over
  threshmat.list <- lapply(seq_len(ncol(threshmat.kstar)), function(i) threshmat.kstar[,i])
  inflmat.list <- lapply(seq_len(ncol(inflmat)), function(i) inflmat[,i])

  # Get thresholds for each data point
  thresholds <- as.data.frame(
    do.call(rbind,
            mapply(get.int,
                   x = threshmat.list,
                   inflmat = inflmat.list,
                   MoreArgs = list(
                     kstar = kstar,
                     trt.code = trt.code,
                     contrs = d_ab[contr.kstar,],
                     mcid = mcid.type == "decision" & mcid > 0,
                     mean.dk = mean.dk,
                     opt.max = opt.max
                     ),
                   SIMPLIFY = FALSE)))


## Return thresh object ----------------------------------------------------
  return(structure(
    list(thresholds = thresholds,
         U = threshmat,
         Ukstar = threshmat.kstar,
         H = inflmat,
         kstar = trt.code[kstar],
         call = list(
           mean.dk = mean.dk,
           lhood = lhood,
           post = post,
           nmatype = nmatype,
           X = X,
           mu.design = mu.design,
           delta.design = delta.design,
           opt.max = opt.max,
           trt.rank = trt.rank,
           trt.code = trt.code,
           trt.sub = trt.sub,
           mcid = mcid,
           mcid.type = mcid.type
           )
         ),
    class="thresh")
    )

}



## Function get.int to return thresholds from U ----------------------------

# Return the positive and negative thresholds for each observation
# Define a function to do this which we can apply over the columns of U (observations)

#' Get thresholds from U matrix
#'
#' Return the positive and negative thresholds for an observation, given a
#' vector of possible threshold solutions. This function is intended for
#' internal use, and is called by \code{nma_thresh} automatically.
#'
#' @param x Column of \eqn{U} matrix, containing all possible threshold
#'   solutions for a data point.
#' @param kstar Base-case optimal treatment.
#' @param trt.code Vector of (possibly recoded) treatments. See
#'   \code{nma_thresh} parameter of the same name.
#' @param contrs Details of contrasts corresponding to rows in \code{x},
#'   as rows of the data.frame output by \code{d_i2ab}.
#' @param mcid Use MCID decision rule? Default \code{FALSE}.
#' @param mean.dk Posterior means of basic treatment parameters, required when
#'   \code{mcid} is \code{TRUE}.
#' @param inflmat Column of influence matrix \eqn{H} for the data point,
#'   required when \code{mcid} is \code{TRUE}.
#' @param opt.max Is the maximum treatment effect optimal? See
#'   \code{nma_thresh} parameter of same name. Required when \code{mcid} is
#'   \code{TRUE}.
#'
#' @return Data frame of thresholds and new optimal treatments with columns
#'   \code{lo}, \code{lo.newkstar}, \code{hi}, and \code{hi.newkstar}.
#' @export
#'
get.int <- function(x, kstar, trt.code, contrs, mcid = FALSE,
                    mean.dk = NULL, inflmat = NULL, opt.max = NULL) {

  # Basic parameter checks
  if (mcid == TRUE & (is.null(mean.dk) | is.null(opt.max) | is.null(inflmat))) {
    stop("Provide mean.dk, inflmat, and opt.max when mcid = TRUE")
  }

  # Using standard decision rule
  if (mcid == FALSE & length(kstar) == 1) {
    # If both thresholds are infinite
    if (all(is.infinite(x))) {
      hi <- Inf
      lo <- -Inf
      hi.newkstar <- lo.newkstar <- NA_character_

    # If lower threshold is infinite
    } else if (all(x[!is.infinite(x)] > 0)) {
      hi <- min(x[!is.infinite(x)])
      i.hi <- which(x == hi)
      hi.newkstar <- trt.code[contrs[i.hi, contrs[i.hi, ] != kstar]]
      lo <- -Inf
      lo.newkstar <- NA_character_

    # If upper threshold is infinite
    } else if (all(x[!is.infinite(x)] < 0)) {
      hi <- Inf
      hi.newkstar <- NA_character_
      lo <- max(x[!is.infinite(x)])
      i.lo <- which(x == lo)
      lo.newkstar <- trt.code[contrs[i.lo, contrs[i.lo, ] != kstar]]

    # If neither threshold is infinite
    } else {
      hi <- min(x[x > 0 & !is.infinite(x)])
      i.hi <- which(x == hi)
      hi.newkstar <- trt.code[contrs[i.hi, contrs[i.hi, ] != kstar]]
      lo <- max(x[x < 0 & !is.infinite(x)])
      i.lo <- which(x == lo)
      lo.newkstar <- trt.code[contrs[i.lo, contrs[i.lo, ] != kstar]]
    }
  } else if (mcid == TRUE) {
    # Using MCID decision rule
    # Thresholds fit one of two criteria:
    #   1. Either d_1k = mcid, or
    #   2. d_ak* = mcid, where k* is the maximally effective treatment at the threshold
    # For d_ab contrasts there is no threshold if b is not the most effective.
    # We shall step along the solution vector to check these criteria.

    # Vector of treatments (possibly subset)
    dk <- sort(unique(c(contrs$a, contrs$b)))

    # First, negative thresholds
    i.neg <- order(x[x < 0 & !is.infinite(x)], decreasing = TRUE)
    if (length(i.neg) == 0) {
      # No negative threshold
      lo <- -Inf
      lo.newkstar <- NA_character_

    } else {
      x.neg <- x[x < 0 & !is.infinite(x)][i.neg]
      ab.neg <- contrs[x < 0 & !is.infinite(x),][i.neg,]

      for (i in 1:length(i.neg)) {
        if (ab.neg[i, "a"] == 1) {
          # Treatment has dropped below MCID
          lo <- x.neg[i]
          lo.newkstar <- trt.code[setdiff(kstar, ab.neg[i, "b"])]
          break
        } else {
          # Evaluate treatment effects at threshold
          d.thr <- c(0, mean.dk + x.neg[i] * inflmat)[dk]

          # Which is the best?
          if (opt.max) kstar.thr <- dk[d.thr == max(d.thr)]
          else kstar.thr <- dk[d.thr == min(d.thr)]

          # Does the contrast involve the best treatment?
          if (any(ab.neg[i,] == kstar.thr)) {
            # If both treatments were in kstar, one has dropped out
            if (all(ab.neg[i,] %in% kstar)) {
              lo <- x.neg[i]
              lo.newkstar <- trt.code[setdiff(kstar, setdiff(ab.neg[i,], kstar.thr))]
              break
            } else {
              # Otherwise, another has joined kstar
              lo <- x.neg[i]
              lo.newkstar <- trt.code[sort(c(kstar, setdiff(ab.neg[i,], kstar.thr), recursive = TRUE))]
              break
            }
          }
          # Otherwise continue, no threshold
        }
        # If we have reached the end of the loop, no threshold found
        if (i == length(i.neg)) {
          lo <- -Inf
          lo.newkstar <- NA_character_
        }
      }
    }

    # Now positive thresholds
    i.pos <- order(x[x > 0 & !is.infinite(x)], decreasing = FALSE)
    if (length(i.pos) == 0) {
      # No positive threshold
      hi <- Inf
      hi.newkstar <- NA_character_

    } else {
      x.pos <- x[x > 0 & !is.infinite(x)][i.pos]
      ab.pos <- contrs[x > 0 & !is.infinite(x),][i.pos,]

      for (i in 1:length(i.pos)) {
        if (ab.pos[i, "a"] == 1) {
          # Treatment has dropped below MCID
          hi <- x.pos[i]
          hi.newkstar <- trt.code[setdiff(kstar, ab.pos[i, "b"])]
          break
        } else {
          # Evaluate treatment effects at threshold
          d.thr <- c(0, mean.dk + x.pos[i] * inflmat)[dk]

          # Which is the best?
          if (opt.max) kstar.thr <- dk[d.thr == max(d.thr)]
          else kstar.thr <- dk[d.thr == min(d.thr)]

          # Does the contrast involve the best treatment?
          if (any(ab.pos[i,] == kstar.thr)) {
            # If both treatments were in kstar, one has dropped out
            if (all(ab.pos[i,] %in% kstar)) {
              hi <- x.pos[i]
              hi.newkstar <- trt.code[setdiff(kstar, setdiff(ab.pos[i,], kstar.thr))]
              break
            } else {
              # Otherwise, another has joined kstar
              hi <- x.pos[i]
              hi.newkstar <- trt.code[sort(c(kstar, setdiff(ab.pos[i,], kstar.thr), recursive = TRUE))]
              break
            }
          }
          # Otherwise continue, no threshold
        }
        # If we have reached the end of the loop, no threshold found
        if (i == length(i.pos)) {
          hi <- Inf
          hi.newkstar <- NA_character_
        }
      }
    }
  }

  thresholds <- data.frame(lo = lo,
                           lo.newkstar = paste0(lo.newkstar, collapse = ", "),
                           hi = hi,
                           hi.newkstar = paste0(hi.newkstar, collapse = ", "),
                           stringsAsFactors = FALSE)
  rownames(thresholds) <- NULL

  return(thresholds)
}
