#' Draw pedigree parts
#'
#' Draw a plot object made by [partsPed()]. This function only draws the given
#' tables and does not inspect the pedigree.
#'
#' @param parts Output from [partsPed()].
#' @param scale Output from [scalePed()].
#' @param new A logical. If TRUE, start a new plot.
#' @param ... Not used.
#'
#' @return Invisibly returns `parts`.
#' @export
#'
#' @examples
#' x = data.frame(id = 1:3, fid = c(0, 0, 1), mid = c(0, 0, 2), sex = c(1, 2, 1))
#' res = plotPed(x, draw = FALSE)
#' drawPed(res$parts, res$scale)
drawPed = function(parts, scale, new = TRUE, ...) {
  if(new) {
    # Set either line or physical margins before initialising user coordinates
    if(identical(scale$marginType, "mai"))
      graphics::par(mai = scale$margins, xpd = TRUE)
    else
      graphics::par(mar = scale$margins, xpd = TRUE)
    graphics::plot.new()
    graphics::plot.window(xlim = scale$usr[1:2], ylim = scale$usr[3:4],
                          xaxs = "i", yaxs = "i")
  }
  else {
    graphics::par(usr = scale$usr, xpd = TRUE)
  }

  segs = parts$segments
  if(nrow(segs)) {
    graphics::segments(segs$x0, segs$y0, segs$x1, segs$y1,
                       col = segs$col, lty = segs$lty, lwd = segs$lwd)
  }

  curves = parts$curves
  if(nrow(curves)) {
    # Draw each duplicate curve separately
    n = tabulate(curves$curve)
    ends = cumsum(n) # curve rows are contiguous
    starts = ends - n + 1L
    for(i in seq_along(n)) {
      j = starts[i]
      ii = seq.int(j, ends[i])
      graphics::lines.default(curves$x[ii], curves$y[ii], col = curves$col[j],
                              lty = curves$lty[j], lwd = curves$lwd[j])
    }
  }

  syms = parts$symbols
  if(nrow(syms)) {
    for(i in seq_len(nrow(syms))) {
      poly = .shapePoly(syms$shape[i])
      graphics::polygon(syms$x[i] + poly$x * syms$width[i],
                        syms$y[i] + poly$y * syms$height[i],
                        border = syms$col[i], col = syms$fill[i],
                        lty = syms$lty[i], lwd = syms$lwd[i])
    }

    # Draw hatches only for selected symbols
    hatch = which(syms$hatched & !is.na(syms$hatchCol))
    for(i in hatch) {
      poly = .shapePoly(syms$shape[i])
      x = syms$x[i] + poly$x * syms$width[i]
      y = syms$y[i] + poly$y * syms$height[i]
      graphics::polygon(x, y, border = NA, col = syms$hatchCol[i],
                        density = parts$hatchDensity, angle = parts$hatchAngle,
                        lty = 1, lwd = parts$hatchLwd)

      # Restore the border independently of the hatch width
      graphics::polygon(x, y, border = syms$col[i], col = NA,
                        lty = syms$lty[i], lwd = syms$lwd[i])
    }

    # Draw carrier dots at symbol centres
    carrier = which(syms$carrier)
    if(length(carrier))
      graphics::points(syms$x[carrier],
                       syms$y[carrier] + syms$height[carrier] / 2,
                       pch = 16, cex = scale$cex, col = syms$col[carrier])

    # Draw deceased slashes over symbols
    deceased = which(syms$deceased)
    if(length(deceased)) {
      x = syms$x[deceased]
      y = syms$y[deceased]
      w = syms$width[deceased]
      h = syms$height[deceased]
      graphics::segments(x - 0.6 * w, y - 0.1 * h,
                         x + 0.6 * w, y + 1.1 * h,
                         col = syms$col[deceased])
    }

    # Draw proband arrows from the lower left
    proband = which(syms$proband)
    if(length(proband)) {
      x = syms$x[proband]
      y = syms$y[proband]
      w = syms$width[proband]
      h = syms$height[proband]
      graphics::arrows(x - 2.2 * w, y - 0.9 * h, x - w, y,
                       length = 0.15, lwd = 1.2,
                       col = syms$col[proband])
    }

  }

  # ID labels
  txt = parts$texts
  if(nrow(txt)) {
    for(i in seq_len(nrow(txt))) {
      graphics::text.default(txt$x[i], txt$y[i], labels = txt$label[i],
                             cex = txt$cex[i], col = txt$col[i],
                             adj = c(txt$adjx[i], txt$adjy[i]), xpd = NA)
    }
  }

  # Text annotation
  ann = parts$textAnnot
  if(nrow(ann)) {
    # Group annotations sharing non-vectorised arguments
    group = as.integer(interaction(ann$adjx, ann$adjy, ann$family, ann$srt,
                                   drop = TRUE))
    for(g in unique.default(group)) {
      i = which(group == g)
      j = i[1L] # read shared arguments from the first row in the group
      graphics::text.default(
        ann$x[i], ann$y[i], labels = ann$label[i], cex = ann$cex[i],
        adj = c(ann$adjx[j], ann$adjy[j]), col = ann$col[i], font = ann$font[i],
        family = ann$family[j], srt = ann$srt[j], xpd = NA)
    }
  }

  # Title
  tit = parts$title
  if(!is.null(tit))
    do.call(graphics::mtext, c(list(text = tit$text, side = 3), tit$args))

  invisible(parts)
}

.shapePoly = function(shape) {
  switch(shape,
    diamond = list(x = c(0, -0.55, 0, 0.55), y = c(-0.05, 0.5, 1.05, 0.5)),
    square = list(x = c(-0.5, -0.5, 0.5, 0.5), y = c(0, 1, 1, 0)),
    circle = {
      a = seq.int(0, 2 * pi, length.out = 50L)
      list(x = 0.5 * cos(a), y = 0.5 * sin(a) + 0.5)
    },
    triangle = list(x = c(0, -0.6, 0.6), y = c(1, 0.4, 0.4)),
    list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5)))
}
