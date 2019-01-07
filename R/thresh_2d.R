#' Producing two-dimensional invariant regions
#'
#' This function produces two-dimensional threshold lines and invariant regions,
#' as shown by Phillippo \emph{et al.} (2018).
#'
#' @param thresh A \code{thresh} object, as produced by
#'   \code{\link{nma_thresh}}.
#' @param idx Integer specifying the index (with respect to
#'   \code{thresh$thresholds}) of the first data point to consider adjusting.
#'   Will be shown on the x axis.
#' @param idy Integer specifying the index (with respect to
#'   \code{thresh$thresholds}) of the second data point to consider adjusting.
#'   Will be shown on the y axis.
#' @param xlab Character string giving the label for the x axis.
#' @param ylab Character string giving the label for the y axis.
#' @param xlim Numeric vector of length 2, giving the x axis limits.
#' @param ylim Numeric vector of length 2, giving the y axis limits.
#' @param breaks Numeric vector giving position of tick marks on the x and y
#'   axes. Calculated automatically by default.
#' @param xbreaks Numeric vector giving position of tick marks on the x axis.
#'   Equal to \code{breaks} by default, if set this overrides any value given to
#'   \code{breaks}.
#' @param ybreaks Numeric vector giving position of tick marks on the y axis.
#'   Equal to \code{breaks} by default, if set this overrides any value given to
#'   \code{breaks}.
#' @param fill Fill colour for invariant region. Defaults to a nice shade of
#'   blue \code{rgb(.72, .80, .93, .7)}.
#' @param lwd Line width for threshold lines. Default 1.
#' @param fontsize Font size for labels. Default 12.
#'
#' @import ggplot2
#'
#' @return A \code{ggplot} object containing the 2D threshold plot, which is
#'   returned invisibly and plotted (unless assigned).
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
#' # Produce an invariant region for simultaneous adjustments to both arms of Study 1
#' thresh_2d(thresh, 1, 2,
#'           xlab = "Adjustment in Study 1 LOR: 3 vs. 1",
#'           ylab = "Adjustment in Study 1 LOR: 4 vs. 1",
#'           xlim = c(-1.5, 0.5), ylim = c(-2, 14),
#'           ybreaks = seq(-2, 14, 2))
#'
thresh_2d <- function(thresh, idx, idy,
                      xlab = paste("Adjustment to data point", idx), ylab = paste("Adjustment to data point", idy),
                      xlim = NULL, ylim = NULL,
                      breaks = waiver(), xbreaks = breaks, ybreaks = breaks,
                      fill = rgb(.72, .80, .93, .7),
                      lwd = 1,
                      fontsize = 12){

  # Cannot handle mcid decisions / vector kstar yet
  if (length(thresh$kstar) > 1 ||
      (!is.null(thresh$call$mcid) &&
       thresh$call$mcid > 0 &&
       thresh$call$mcid.type == 'decision')) {
    stop("Decision rules with multiple optimal treatments not yet supported.")
  }

  # Number of treatments
  K <- nrow(thresh$Ukstar) + 1

  # Derive intercept and gradient of threshold lines using Ukstar matrix
  linedat <- data.frame(intercept = thresh$Ukstar[, idy],
                        gradient = - thresh$Ukstar[, idy] / thresh$Ukstar[, idx])


  # Intercepts array [line1, line2, x/y]
  l1s <- rep(1:(K-1), each = K-1)
  l2s <- rep(1:(K-1), times = K-1)

  intarray <- array(NA, dim = c(K-1 , K-1, 2))

  intarray[,,1] <- (linedat[l2s, "intercept"] - linedat[l1s, "intercept"]) /
                      (linedat[l1s, "gradient"] - linedat[l2s, "gradient"])

  diag(intarray[,,1]) <- NA   # intercept of line with itself

  intarray[,,2] <- linedat[, "gradient"] * intarray[,,1] + linedat[, "intercept"]

  # Intercept array to data frame [x,y]
  uptri <- upper.tri(intarray[,,1])
  intdat <- data.frame(x = intarray[,,1][uptri],
                       y = intarray[,,2][uptri],
                       l1 = matrix(rep(1:(K-1), times = K-1), nrow = K-1)[uptri],
                       l2 = matrix(rep(1:(K-1), each = K-1), nrow = K-1)[uptri])

  # Calculate invariant region
  # For each threshold line:
  #  1. Check which side of the line the origin lies
  #  2. Exclude all intercept points which lie the opposite side
  #  3. Repeat 1-2 for each line

  # Function to determine inclusion/exclusion of points
  IRvertex <- function(Mx,    # x value of test point
                       My,    # y value of test point
                       xint,  # x intercepts of all threshold lines
                       yint,  # y intercepts of all threshold lines
                       eps=1e-12){  # for testing equality

    # Quick parameter checks
    stopifnot(length(xint)==length(yint))
    stopifnot(length(Mx)==1)
    stopifnot(length(My)==1)

    # Is test point the same side of every line as the origin?
    # Ignore lines which a point lies on (within eps)
    all(ifelse(abs((xint - Mx)*(yint - My) - Mx*My) <= eps,
               TRUE,
               sign((xint - Mx)*(yint - My) - Mx*My) == sign(xint * yint)))
  }


  # For each intercept point, check inclusion/exclusion
  inIR <- mapply(IRvertex, Mx = intdat$x, My = intdat$y,
                 MoreArgs = list(xint = thresh$Ukstar[, idx],
                                 yint = thresh$Ukstar[, idy]))


  # Invariant region data
  IRdat <- intdat[inIR, ]

  # Check that the invariant region is closed.
  # If it isn't, add some lines out of view to close it.
  # If no lines are parallel, then there is at most one open side, and we can:
  #   - Check if any lines enter only one vertex
  #   - If so, add a new edge between the two lines, closing the IR

  IRlines <- unique(c(IRdat$l1, IRdat$l2))  # Which lines bound the IR?

  # Set some values for xlim/ylim if not given, rather than leaving it up to
  # ggplot. We also use these values at this stage, so we need them before
  # ggplot is called anyway. Make the view window large enough to contain the
  # vertices of the IR (and the x/y intercepts, in the open case), times a small
  # expansion factor.
  if (is.null(xlim)) {
    xlim <- range(IRdat$x, thresh$Ukstar[IRlines, idx], na.rm = TRUE, finite = TRUE) * 1.1
    if (xlim[1] >= 0) xlim[1] <- -xlim[2]
    if (xlim[2] <= 0) xlim[2] <- -xlim[1]
  }
  if (is.null(ylim)) {
    ylim <- range(IRdat$y, thresh$Ukstar[IRlines, idy], na.rm = TRUE, finite = TRUE) * 1.1
    if (ylim[1] >= 0) ylim[1] <- -ylim[2]
    if (ylim[2] <= 0) ylim[2] <- -ylim[1]
  }

  # Check for parallel boundaries
  if (anyDuplicated(IRdat[IRlines, "gradient"])) {
    warning("Parallel boundary lines of the invariant region, may break shading/labelling for now.")
  }
  # No parallel lines, so see if there are any entering only one vertex
  else if (any(tabulate(c(IRdat$l1, IRdat$l2)) == 1)) {

    # List the edges with an open end
    oedge <- which(tabulate(c(IRdat$l1, IRdat$l2)) == 1)

    # Check which end is open between the two
    # Note: osign is the half of the x-axis where the two edges intersect, so
    # the open end is on the other side.
    osign <- sign(intdat[intdat$l1 %in% oedge & intdat$l2 %in% oedge, "x"])


    # Add a vertex on each open line way past the view window
    # Note: can't use xlim as if NULL is calculated by ggplot
    for (i in 1:2){
      IRdat <- rbind(IRdat,
                     data.frame(x = xlim[ifelse(osign == 1, 1, 2)]*1.2,
                                y = linedat[oedge[i], "gradient"]*xlim[ifelse(osign == 1, 1, 2)]*1.2 + linedat[oedge[i], "intercept"],
                                l1 = oedge[i], l2 = -999))
    }

    # Provide a note to the user that an open invariant region was found
    message("NOTE: The resulting invariant region is open at one side. There is no threshold boundary in this direction.")
  }

  # Sort vertices into clockwise order
  IRdat$angle <- atan2(IRdat$y, IRdat$x)
  IRdat <- IRdat[order(IRdat$angle),]

  # Create line labels data frame
  labdat <- linedat[IRlines, ]
  labdat$l <- IRlines

  # Labels
  labdat$lab <- paste0("paste(tilde(k),'* = ',", labdat$l + (labdat$l >= thresh$kstar), ")")

  # Get min and max x values of visible invariant region bounded lines
  for (i in 1:length(IRlines)){
    labdat[i, "x1"] <- max(xlim[1],
                           min(IRdat[IRdat$l1 == labdat[i, "l"] | IRdat$l2 == labdat[i, "l"], "x"]),
                           min((ylim - labdat[i, "intercept"])/labdat[i, "gradient"]))

    labdat[i, "x2"] <- min(xlim[2],
                           max(IRdat[IRdat$l1 == labdat[i, "l"] | IRdat$l2 == labdat[i, "l"], "x"]),
                           max((ylim - labdat[i, "intercept"])/labdat[i, "gradient"]))
  }

  labdat$x <- rowMeans(labdat[, c("x1", "x2")])
  labdat$y <- labdat$gradient * labdat$x + labdat$intercept

  # Derive alignment of each label to try to avoid overlap with boundary of IR
  labdat$vjust <- ifelse(labdat$intercept < 0, 1, 0)
  labdat$hjust <- ifelse(-labdat$intercept/labdat$gradient < 0, 1, 0)

  # Construct plot
  ggplot() +

    # Axes
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +

    # Threshold lines
    geom_abline(aes_string(intercept = "intercept", slope = "gradient"),
                data = linedat, colour = "grey60") +

    # Invariant region
    geom_polygon(aes_string(x = "x", y = "y"), data = IRdat,
                 fill = fill, colour = "black") +

    # Line labels
    geom_label(aes_string(x = "x", y = "y", label = "lab",
                          vjust = "vjust", hjust = "hjust"),
               data = labdat, parse = TRUE, na.rm = TRUE) +

    # Axis labels
    xlab(xlab) +
    ylab(ylab) +

    # Axis setup
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_x_continuous(breaks = xbreaks) +
    scale_y_continuous(breaks = ybreaks) +

    # Theme
    theme_bw() +
    theme(panel.grid = element_blank())
}
