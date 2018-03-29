#' Producing threshold forest plots
#'
#' This function produces threshold forest plots, overlaying the
#' decision-invariant intervals on the data points and their confidence/credible
#' intervals, as shown by Phillippo \emph{et al.} (2017).
#'
#' @param thresh A \code{thresh} object, as produced by
#'   \code{\link{nma_thresh}}.
#' @param y Data points. Either a column of \code{data}, or a numeric vector.
#' @param CI.lo Confidence/credible interval lower limits. Either a column of
#'   \code{data}, or a numeric vector.
#' @param CI.hi Confidence/credible interval upper limits. Either a column of
#'   \code{data}, or a numeric vector.
#' @param label Row labels (for each data point). Either a column of
#'   \code{data}, or a character vector.
#' @param orderby Variable(s) to order the table rows by. Either a column or
#'   columns of \code{data}, or a vector. By default, the rows are not
#'   reordered. Further arguments and/or multiple ordering columns may be passed
#'   to the function \code{order} by instead providing a \code{list} containing
#'   the arguments to \code{order}.
#' @param data A data frame containing the data points \code{y},
#'   confidence/credible intervals (\code{CI.lo}, \code{CI.hi}), and row labels
#'   \code{labels}. If \code{data} is not provided, the above variables will be
#'   searched for in the calling environment.
#' @param CI.title Title for CI column, default "95\% Confidence Interval".
#' @param label.title Character string giving the heading for the row labels
#'   column.
#' @param y.title Character string giving the heading for the data points
#'   column, default "Mean".
#' @param II.title Title for invariant interval column, default "Invariant
#'   Interval".
#' @param xlab Character string giving the label for the \eqn{x}-axis.
#' @param xlim Numeric vector (length 2) of lower and upper limits for the
#'   \eqn{x}-axis. If not set, tries to choose a sensible default.
#' @param sigfig Number of significant digits to display in the table. Default
#'   3.
#' @param digits Number of decimal places to display in the table. Overrides
#'   \code{sigfig} if set.
#' @param refline \eqn{x} intercept of vertical reference line, if desired.
#' @param clinsig Set the clinical significance level. Rows are marked with a
#'   dagger if they have one or more thresholds less than this value. Not set by
#'   default.
#' @param cutoff A single numeric value or numeric vector pair. Thresholds
#'   larger in magnitude than this value, or lying outside this interval, will
#'   be cut off and marked as NT (no threshold). Not set by default.
#' @param II.colw Colour for "wide" invariant intervals.
#' @param II.cols Colour for "short" invariant intervals.
#' @param II.lwd Line width of invariant intervals. Default 8.
#' @param CI.lwd Line width of confidence/credible intervals. Default 1.
#' @param pointsize Point size for forest plot means. Default 4.
#' @param fontsize Base font size. Default 12.
#' @param xbreaks Position of tick marks on the \eqn{x}-axis as a numeric
#'   vector.
#' @param add.columns Data frame (or matrix, vector) of additional columns to
#'   add to table.
#' @param add.columns.title Optional titles for the additional columns,
#'   otherwise use names from provided data.
#' @param add.columns.after Which column to add the new columns after? Default
#'   adds the columns to the far right.
#' @param add.columns.hjust Vector of horizontal justifications for the new
#'   columns, from \code{0} (left) to \code{1} (right). Default centres every
#'   column.
#' @param add.columns.uline Underline the header of the new columns? Default
#'   \code{TRUE}.
#' @param calcdim Logical, calculate suggested output dimensions for saving to
#'   pdf? Calculates output size when \code{TRUE}; saves time when \code{FALSE}.
#' @param display Logical, display the plot? Defaults to \code{TRUE}.
#'
#' @return Displays the forest plot on the current plot device (if \code{display
#'   = TRUE}). Also returns invisibly the underlying \code{gtable} object, which
#'   can be further manipulated.
#' @export
#'
#' @import gtable grid gridExtra grDevices
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
#' # Get treatment codes for the contrasts with data
#' d.a <- d.b <- vector(length = nrow(X))
#' for (i in 1:nrow(X)){
#'   d.a[i] <- ifelse(any(X[i, ] == -1), which(X[i, ] == -1), 0) + 1
#'   d.b[i] <- ifelse(any(X[i, ] == 1), which(X[i, ] == 1), 0) + 1
#' }
#'
#' # Transform from d_ab style contrast references into d[i] style from the full set of contrasts
#' # for easy indexing in R
#' d.i <- d_ab2i(d.a, d.b, K = 6)
#'
#' # Create plot data
#' plotdat <- data.frame(lab = paste0(d.b, " vs. ", d.a),
#'                       contr.mean = Thrombo.post.summary$statistics[d.i, "Mean"],
#'                       CI2.5 = Thrombo.post.summary$quantiles[d.i, "2.5%"],
#'                       CI97.5 = Thrombo.post.summary$quantiles[d.i, "97.5%"])
#'
#' # Plot
#' thresh_forest(thresh, contr.mean, CI2.5, CI97.5, label = lab, data = plotdat,
#'               label.title = "Contrast", xlab = "Log Odds Ratio", CI.title = "95% Credible Interval",
#'               xlim = c(-.3, .3), refline = 0, digits = 2)
#'

thresh_forest <- function(thresh,
                         y, CI.lo, CI.hi, label, orderby = NULL, data = NULL,
                         CI.title = "95% Confidence Interval",
                         label.title = "", y.title = "Mean",
                         II.title = "Invariant Interval", xlab = "",
                         xlim = NULL, sigfig = 3, digits = NULL,
                         refline = NULL, clinsig = NULL, cutoff = NULL,
                         II.colw = rgb(.72, .80, .93),
                         II.cols = rgb(.93, .72, .80),
                         II.lwd = 8, CI.lwd = 1,
                         pointsize = 4, fontsize = 12,
                         xbreaks = NULL,
                         add.columns = NULL,
                         add.columns.title = NULL,
                         add.columns.after = -1,
                         add.columns.hjust = 0.5,
                         add.columns.uline = TRUE,
                         calcdim = display, display = TRUE){

  # Evaluate data arguments
  y <- eval(substitute(y), data, parent.frame())
  CI.lo <- eval(substitute(CI.lo), data, parent.frame())
  CI.hi <- eval(substitute(CI.hi), data, parent.frame())
  label <- eval(substitute(label), data, parent.frame())

  # Number of data points
  N <- nrow(thresh$thresholds)

  # Set up orderby
  orderby <- eval(substitute(orderby), data, parent.frame())
  if (is.null(orderby)) {
    orderby <- list(1:N)
  } else if (!is.list(orderby)) {
    orderby <- list(orderby)
  }
  row_order <- do.call("order", orderby)

  # Set up additional columns
  if (!is.null(add.columns)) {
    add.columns <- as.data.frame(add.columns, stringsAsFactors = FALSE)

    if (is.null(add.columns.title)) {
      add.columns.title <- colnames(add.columns)
    }
    if (length(add.columns.title) != ncol(add.columns)) {
      stop("Mismatch number of additional columns and titles")
    }

    add.columns.hjust <- rep_len(add.columns.hjust, ncol(add.columns))
  }

  # Check number of rows
  stopifnot(length(y) == N,
            length(CI.lo) == N,
            length(CI.hi) == N,
            length(label) == N)

  if (!is.null(add.columns) && nrow(add.columns) != N) stop("Mismatch number of rows for add.columns")

  # Set up plot data
  pd <- cbind(data.frame(y = y, CI.lo = CI.lo, CI.hi = CI.hi, label = label),
              thresh$thresholds)

  # Use significant figures or decimal places?
  sf <- is.null(digits)

  # Set up cutoff
  if (is.null(cutoff)) cutoff <- c(-Inf, Inf)
  else if (length(cutoff) <= 1) cutoff <- c(-abs(cutoff), abs(cutoff))

  # Set up xlim and xbreaks
  if (is.null(xbreaks)) {
    if(is.null(xlim)) xbreaks <- pretty(c(min(CI.lo, na.rm = TRUE), max(CI.hi, na.rm = TRUE)))
    else xbreaks <- pretty(xlim)
  }

  if (is.null(xlim)) {
    xlim <- range(xbreaks)
  }

  # Invariant interval
  pd$II.lo <- pd$y + pd$lo
  pd$II.hi <- pd$y + pd$hi

  # If y and lo/hi are infinite, reset II limit to Inf (not NaN)
  pd$II.lo[is.infinite(pd$y) & is.infinite(pd$lo)] <- -Inf
  pd$II.hi[is.infinite(pd$y) & is.infinite(pd$hi)] <- Inf


  # Function to print with significant digits or decimal places
  printsig <- function(num, dig = ifelse(sf, sigfig, digits), cutyn = FALSE){
    if (cutyn && (num >= min(Inf, cutoff[2]) || num <= max(-Inf, cutoff[1]))) "NT"
    else if (sf) gsub("\\.$", "", formatC(signif(num, digits = dig), digits = dig, format = "fg", flag = "#"))
    else formatC(num, dig, format = "f")
  }


  # Present CIs
  for (i in 1:N) {
    pd[i, "CI.txt"] <- paste0("(", ifelse(is.na(pd[i, "CI.lo"]), "\u2013", printsig(pd[i, "CI.lo"])), ", ",
                              ifelse(is.na(pd[i, "CI.hi"]), "\u2013", printsig(pd[i, "CI.hi"])), ")")
  }

	# Format means
	pd$y.txt <- printsig(pd$y)

  # Present invariant intervals
  for (i in 1:N) {
    pd[i, "II.txt"] <- paste0("(", printsig(pd[i, "II.lo"], cutyn = TRUE), ", ",
                                  printsig(pd[i, "II.hi"], cutyn = TRUE), ")")
  }

  # If no thresholds found (or beyond cutoff), set newkstar to "-"
  pd$lo.newkstar[is.na(pd$lo.newkstar) | pd$II.lo <= cutoff[1]] <- "\u2013"
  pd$hi.newkstar[is.na(pd$hi.newkstar) | pd$II.hi >= cutoff[2]] <- "\u2013"

  # Get details of short (statistically insignificant) intervals
  pd$is.short <- (!is.na(pd$CI.lo) & pd$II.lo > pd$CI.lo) | (!is.na(pd$CI.hi) & pd$II.hi < pd$CI.hi)
  pd$lab.ff <- ifelse(pd$is.short, "bold", "plain")

  # Get details of clinically insignificant intervals
  if (!is.null(clinsig)) {
    pd$is.clinshort <- pmin(abs(pd$lo), abs(pd$hi)) < clinsig
  } else {
    pd$is.clinshort <- FALSE
  }

  # Mark these rows with daggers
  pd$lab.clinshort <- rep("", N)
  pd$lab.clinshort[pd$is.clinshort] <- "\206"

  # Reorder rows
  pd <- pd[row_order, ]

  # Save the number of rows again, in case altered by orderby (e.g. with
  # na.last = NA.)
  Nrows <- nrow(pd)

  # Make expression II.title bold
  if (is.expression(II.title)) II.title <- eval(bquote(expression(bold(.(II.title[[1]])))))

  # Arrange table
  g_tab <- tableGrob(
    d = pd[, c("lab.clinshort", "label", "y.txt", "CI.txt",
            "lo.newkstar", "II.txt", "hi.newkstar")],
    rows = NULL,
    cols = c("", label.title, y.title, CI.title, "", "", ""),
    theme = gridExtra::ttheme_minimal(
      base_size = fontsize,
      core = list(fg_params = list(
        hjust = rep(c(0, 0, 0, .5, .5, .5, .5), each = Nrows),
        x = rep(c(.6, 0,.5, .5, .5, .5, .5), each = Nrows),
        fontface = c(rep("plain", Nrows), pd$lab.ff, rep("plain", Nrows*5)),
        vjust = rep(c(1, .5, .5, .5, .5, .5, .5), each = Nrows),
        y = rep(c(1, .5, .5, .5, .5, .5, .5), each = Nrows)
        )),
      colhead = list(fg_params = list(
        hjust = c(0, 0, .5, .5, .5, .5, .5),
        vjust = c(0, 0, 0, 0, 0, 0, 0),
        x = c(0, 0,.5, .5, .5, .5, .5),
        y = c(.25, .25, .25, .25, .25, .25, .25)
        ))
      )
    )

  # Add II header separately, so that it spans the newkstar columns as well
  II.title.grob <- textGrob(label = II.title,
                            y = 0.25, hjust = 0.5, vjust = 0,
                            gp = g_tab$grobs[[4]]$gp)  # Copy gp from CI column header
  g_tab <- gtable_add_grob(g_tab,
                           grobs = II.title.grob,
                           t = 1, b = 1, l = 5, r = 7)

  # If II header is really long, add extra width to either side of newkstar columns
  extrawidth <- max(unit(1, "grobwidth", II.title.grob) -
                      sum(g_tab$widths[[5]], g_tab$widths[[6]], g_tab$widths[[7]]),
                    unit(0, "mm")) * 0.5

  g_tab$widths[[5]] <- g_tab$widths[[5]] + extrawidth
  g_tab$widths[[7]] <- g_tab$widths[[7]] + extrawidth

  # Align the mean column at the decimal point
  # Get string widths to left of dp
  y.strwd <- lapply(gsub("(-?[0-9]+[.]?)[0-9]+", "\\1", pd$y.txt), function(x) unit(1, "strwidth", x))

  # Find grobs for mean column, update x location
  meancol <- which(g_tab$layout$l == 3 & g_tab$layout$name == "core-fg")
  for (i in 1:length(meancol)) {
    g_tab$grobs[[meancol[i]]]$x <- unit(0.5, "npc") - y.strwd[[i]]
  }

  # Set number of table columns
  Ntabcols <- 7

  # Add in rows of forest plot to the table
  g_all <- gtable_add_cols(g_tab, unit(1, "null"))

  g_for <- mapply(forestgrob,
                  pd$y, pd$CI.lo, pd$CI.hi, pd$II.lo, pd$II.hi,
                  MoreArgs = list(xlim, CI.lwd, II.cols, II.colw, II.lwd, pointsize),
                  SIMPLIFY = FALSE)

  g_all <- gtable_add_grob(g_all, g_for, t = 2:(Nrows+1), l = Ntabcols+1, z = 100)

  # Add reference line
    if (!is.null(refline)) {
    g_all <- gtable_add_grob(g_all, linesGrob(x = x2c(rep(refline, 2), xlim),
                                              gp = gpar(lty = 2, col = "grey40")),
                             t = 2, b = Nrows+1, l = Ntabcols+1, r = Ntabcols+1, z = 99)
  }

  # Add x-xaxis
  g_all <- gtable_add_rows(g_all, heights = unit(2, "lines"))
  g_all <- gtable_add_grob(g_all, gTree(children = gList(xaxisGrob(x2c(xbreaks, xlim),
                                                                   label = xbreaks,
                                                                   gp = gpar(fontsize = fontsize - 1))),
                                        vp = viewport(y = 1, just = "bottom", height = unit(2, "lines"))
                                        ),
                           t = Nrows+2, b = Nrows+2, l = Ntabcols+1, r = Ntabcols+1, z = 102, clip = "off")

  # Add x-xaxis label
  g_all <- gtable_add_rows(g_all, heights = unit(2, "lines"))
  g_all <- gtable_add_grob(g_all, textGrob(xlab, vjust = 0.5, gp = gpar(fontsize = fontsize - 1)),
                           t = Nrows+3, b = Nrows+3, l = Ntabcols+1, r = Ntabcols+1, z = 101, clip = "off")

  # Add in additional columns (if any)
  if (!is.null(add.columns)) {
    # Reorder if necessary
    add.columns <- add.columns[row_order, , drop = FALSE]

    # Format numeric columns
    add.columns[] <- lapply(add.columns, function(x){if (is.numeric(x)) printsig(x) else x})

    # Create table grob
    g_add <- tableGrob(add.columns,
                       rows = NULL, cols = add.columns.title,
                       theme = gridExtra::ttheme_minimal(
                         base_size = fontsize,
                         core = list(fg_params = list(
                           hjust = rep(add.columns.hjust, each = Nrows),
                           x = rep(add.columns.hjust*0.9 + 0.05, each = Nrows)
                         )),
                         colhead = list(fg_params = list(
                           hjust = add.columns.hjust,
                           x = add.columns.hjust*0.9 + 0.05
                         ))
                       ))

    # Add blank rows to account for axes
    g_add <- gtable_add_rows(g_add, heights = g_all$heights[-(1:nrow(g_add))], pos = -1)

    # Add zero width column to fix underline bug with only one add.column
    if (ncol(add.columns) == 1 && (add.columns.after == -1 || add.columns.after > Ntabcols)) {
      g_add <- gtable_add_cols(g_add, widths = unit(0, "npc"), pos = -1)
    }

    # Add extra columns at desired position
    if (add.columns.after == -1) {
      g_all <- cbind(g_all, g_add, size = "first", z = c(0, 1))
    } else {
      g_all <- cbind(g_all[, 1:add.columns.after], g_add, g_all[, -(1:add.columns.after)],
                     size = "first", z = c(0, 1, 2))
    }

    # Update Ntabcols, if necessary
    if (add.columns.after != -1 && add.columns.after <= Ntabcols) Ntabcols <- Ntabcols + ncol(add.columns)
  }

  # Add legend manually (constructed as an inset table)
  leg <- tableGrob(matrix(c("   ", paste0(y.title, "   "),
                            "   ", paste0(CI.title, "   "),
                            "   ", "Invariant Interval"), nrow = 1),
                   rows = NULL, cols = NULL, theme = ttheme_minimal(base_size = fontsize - 1))

  leg <- gtable_add_grob(leg, circleGrob(x = .5, y = .5, r = unit(pointsize, "pt")),
                         t = 1, l = 1, z = -200)
  leg <- gtable_add_grob(leg, linesGrob(x = c(0, 1), y = c(.5, .5),
                                        gp = gpar(lwd = CI.lwd)),
                         t = 1, l = 3, z = -201)
  leg <- gtable_add_grob(leg, linesGrob(x = c(0, 1), y = c(.5, .5),
                                        gp = gpar(lwd = unit(II.lwd, "pt"), col = II.colw)),
                         t = 1, l = 5, z = -202)

  g_all <- gtable_add_grob(g_all, leg, t = Nrows+3, b = Nrows+3, l = 2, r = Ntabcols, z = 98)

  # Header underline
  rightul <- (!is.null(add.columns)) &&
    add.columns.uline &&
    (add.columns.after == -1 || add.columns.after > Ntabcols)

  ulgrobs <- replicate(ifelse(rightul, 2, 1),
                       segmentsGrob(
                         x0 = unit(0, "npc"), y0 = unit(0, "npc"),
                         x1 = unit(1, "npc"), y1 = unit(0, "npc"),
                         gp = gpar(lwd = 1)),
                       simplify = FALSE)

  g_all <- gtable_add_grob(g_all,
                           grobs = ulgrobs,
                           t = 1, b = 1,
                           l = if (rightul) c(1, Ntabcols + 2) else 1,
                           r = if (rightul) c(Ntabcols, -1) else Ntabcols,
                           z = Inf)

  # g_all <- gtable_add_grob(g_all,
  #                          grobs = segmentsGrob(
  #                            x0 = unit(0, "npc"), y0 = unit(0, "npc"),
  #                            x1 = unit(1, "npc"), y1 = unit(0, "npc"),
  #                            gp = gpar(lwd = 1)
  #                            ),
  #                          t = 1, b = 1, l = 1, r = Ntabcols, z = -2)
  #
  # if (add.columns.uline & (add.columns.after == -1 | add.columns.after > Ntabcols)) {
  #   g_all <- gtable_add_grob(g_all,
  #                            grobs = segmentsGrob(
  #                              x0 = unit(0, "npc"), y0 = unit(0, "npc"),
  #                              x1 = unit(1, "npc"), y1 = unit(0, "npc"),
  #                              gp = gpar(lwd = 1)
  #                              ),
  #                            t = 1, b = 1, l = Ntabcols + 2, r = ncol(g_all), z = -1)
  # }

  # Add padding in between table and plot
  g_all <- gtable_add_cols(g_all, unit(1, "lines"), Ntabcols)

  if (!is.null(add.columns) && add.columns.after == -1 || add.columns.after > Ntabcols) {
    g_all <- gtable_add_cols(g_all, unit(1, "lines"), Ntabcols + 2)
  }

  # Pad edges of table
  g_all <- gtable_add_padding(g_all, unit(c(.5, 1, .5, .5), "lines"))

  # Plot
  if (display) {
    grid.newpage()
    grid.draw(g_all)
  }

  # Print sizes for easy saving
  if (calcdim) {
    message("Suggested output dimensions (inches)")
    message("   width:", convertWidth(sum(g_all$widths), "in", valueOnly = TRUE)*2)
    message("   height:", convertHeight(sum(g_all$heights), "in", valueOnly = TRUE))
  }

  # Invisibly return the gtable object
  invisible(g_all)
}


### Internal functions ###

# Function to convert x value to coordinate in (0, 1)
x2c <- function(x, xlim){
  if (is.null(x)) return(NULL)
  (x - xlim[1])/(xlim[2] - xlim[1])
}

# Function to create grid grob of mean, CI, and threshold interval
forestgrob <- function(y,
                       CI.lo, CI.hi, II.lo, II.hi,
                       xlim, CI.lwd, II.cols, II.colw, II.lwd, pointsize){

  # If y is infinite/NA then nothing is defined, return empty gTree
  if (is.infinite(y) || is.na(y)) return(grobTree(NULL))

  # Truncate CIs, IIs if necessary
  if (is.na(CI.lo)) {
    tCI.lo <- NA
    istCI.lo <- FALSE
  } else {
    tCI.lo <- max(CI.lo, xlim[1])
    istCI.lo <- CI.lo < xlim[1]
  }

  if (is.na(CI.hi)) {
    tCI.hi <- NA
    istCI.hi <- FALSE
  } else {
    tCI.hi <- min(CI.hi, xlim[2])
    istCI.hi <- CI.hi > xlim[2]
  }

  tII.lo <- max(II.lo, xlim[1])
  istII.lo <- II.lo < xlim[1]

  tII.hi <- min(II.hi, xlim[2])
  istII.hi <- II.hi > xlim[2]

  # Half invariant interval line width for drawing
  hlwd <- II.lwd / 2

  # Construct grob
  grobTree(
    # Invariant interval, lower
    {if (istII.lo) {
      polygonGrob(x = unit(x2c(c(y, tII.lo, tII.lo, tII.lo, y), xlim),
                           units = "npc") +
                    unit(c(0, hlwd, 0, hlwd, 0), units = "pt"),
                  y = unit(rep(0.5, 5), units = "npc") +
                    unit(c(hlwd, hlwd, 0, -hlwd, -hlwd), units = "pt"),
                  gp = gpar(lwd = NA, fill = ifelse(!is.na(CI.lo) && II.lo >= CI.lo, II.cols, II.colw)))
    } else {
      polygonGrob(x = x2c(c(y, tII.lo, tII.lo, y), xlim),
                  y = unit(rep(0.5, 4), units = "npc") +
                    unit(c(hlwd, hlwd, -hlwd, -hlwd), units = "pt"),
                  gp = gpar(lwd = NA, fill = ifelse(!is.na(CI.lo) && II.lo >= CI.lo, II.cols, II.colw)))
    }},

    # Invariant interval, upper
    {if (istII.hi) {
      polygonGrob(x = unit(x2c(c(y, tII.hi, tII.hi, tII.hi, y), xlim),
                           units = "npc") +
                    unit(c(0, -hlwd, 0, -hlwd, 0), units = "pt"),
                  y = unit(rep(0.5, 5), units = "npc") +
                    unit(c(hlwd, hlwd, 0, -hlwd, -hlwd), units = "pt"),
                  gp = gpar(lwd = NA, fill = ifelse(!is.na(CI.hi) && II.hi <= CI.hi, II.cols, II.colw)))
    } else {
      polygonGrob(x = x2c(c(y, tII.hi, tII.hi, y), xlim),
                  y = unit(rep(0.5, 4), units = "npc") +
                    unit(c(hlwd, hlwd, -hlwd, -hlwd), units = "pt"),
                  gp = gpar(lwd = NA, fill = ifelse(!is.na(CI.hi) && II.hi <= CI.hi, II.cols, II.colw)))
    }},

    # CI
    linesGrob(x = x2c(c(tCI.lo, tCI.hi), xlim), y = c(.5, .5),
              arrow = {
                if (!istCI.lo && !istCI.hi) {
                  NULL
                } else if (istCI.lo && istCI.hi) {
                  arrow(ends = "both", type = "closed", length = unit(6, "pt"))
                } else if (istCI.lo) {
                  arrow(ends = "first", type = "closed", length = unit(6, "pt"))
                } else {
                  arrow(ends = "last", type = "closed", length = unit(6, "pt"))
                }
              },
              gp = gpar(lwd = CI.lwd, fill = "black")),

    # Mean
    circleGrob(x = x2c(y, xlim), y = .5, r = unit(pointsize, "pt"),
               gp = gpar(col = "black", fill = "white")),
    vp = viewport(clip = "off")
  )

}
