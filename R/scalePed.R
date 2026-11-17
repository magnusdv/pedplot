#' Scale a pedigree plot
#'
#' Resolve symbol size, label dimensions and plot limits for an aligned
#' pedigree.
#'
#' `scalePed()` converts an alignment from [alignPed()] into sizes and limits
#' used when drawing the pedigree. Lengths are computed internally in
#' inches, matching base graphics text and device measurements, while user input
#' may be given in other units.
#'
#' Automatic scaling first measures labels and duplication curves, then chooses
#' a symbol size fitting the available device region after margins. Symbol boxes,
#' labels below symbols, side labels and duplication curves determine the plot
#' limits. Extra `padding` may be added to give space for additional external
#' annotations.
#'
#'
#' @param align Output from [alignPed()].
#' @param labs Labels printed below symbols.
#' @param textAnnot A named list of text annotations around or inside symbols.
#' @param proband IDs of proband individuals, used to reserve arrow space.
#' @param title Optional title printed in the top margin.
#' @param titleArgs A list of arguments passed to [graphics::mtext()] for the
#'    title.
#' @param cex Text size multiplier for labels.
#' @param symbolsize Symbol size. With `symbolsizeUnit = "auto"` this is a
#'   multiplier of the automatic size. With `"in"`, `"cm"` or `"mm"` it is an
#'   absolute size.
#' @param symbolsizeUnit Unit for `symbolsize`: `"auto"` (default), `"in"`,
#'   `"cm"` or `"mm"`.
#' @param margins Plot margins, recycled to bottom, left, top, right.
#' @param marginsUnit Unit for `margins`: `"lines"` (default), `"in"`, `"cm"`
#'   or `"mm"`. Line margins use `par(mar)`, physical margins use `par(mai)`.
#' @param padding Extra inner plot space for annotations, recycled to bottom,
#'   left, top, right.
#' @param paddingUnit Unit for `padding`: `"in"` (default), `"cm"`, `"mm"` or
#'   `"sym"` (= size of a pedigree symbol).
#' @param xlim,ylim Optional numeric vectors of length 2 giving plot limits in
#'   user coordinates.
#' @param ... Further arguments, currently ignored.
#'
#' @return A list with scaling values, including:
#'
#' * `boxw`, `boxh`: Symbol width and height in user coordinates.
#' * `boxIn`: Symbol size in inches.
#' * `labh`, `labhIn`: Single-line label height in user coordinates and inches.
#' * `legh`: Length of short connector legs in user coordinates.
#' * `cex`: Text size multiplier.
#' * `margins`: Plot margins.
#' * `marginType`: Graphics margin type, either `"mar"` or `"mai"`.
#' * `paddingIn`: Annotation padding in inches, named bottom, left, top, right.
#' * `hscale`, `vscale`: Inches per user-coordinate unit.
#' * `usr`: Plot limits, named left, right, bottom, top.
#'
#' @examples
#' x = data.frame(id = 1:3, fid = c(0, 0, 1),
#'                mid = c(0, 0, 2), sex = c(1, 2, 1))
#' align = alignPed(x)
#' scalePed(align)
#' scalePed(align, symbolsize = 4, symbolsizeUnit = "mm")
#' scalePed(align, padding = c(1, 0, 0, 0), paddingUnit = "sym")
#'
#' @export
#' @importFrom graphics par strheight strwidth
scalePed = function(align, labs = align$input$id, cex = 1, symbolsize = 1,
                    symbolsizeUnit = c("auto", "in", "cm", "mm"),
                    margins = 1, marginsUnit = c("lines", "in", "cm", "mm"),
                    padding = 0, paddingUnit = c("in", "cm", "mm", "sym"),
                    xlim = NULL, ylim = NULL, textAnnot = NULL,
                    proband = NULL, title = NULL, titleArgs = list(), ...) {

  symbolsizeUnit = match.arg(symbolsizeUnit)
  marginsUnit = match.arg(marginsUnit)
  paddingUnit = match.arg(paddingUnit)

  if(length(symbolsize) != 1L)
    stop2("`symbolsize` must have length 1")

  nodes = align$nodes
  if(!nrow(nodes))
    stop2("No visible individuals")

  textUnder = .prepLabs(align$input, labs)
  labsNode = textUnder[nodes$idx]
  labSize = .textSizeIn(labsNode, cex)
  labhIn = labSize$lineh
  textAnnot = .prepTextAnnot(align, textAnnot, cex)

  # Single-parent labels are drawn beside symbols rather than below them
  matings = align$matings
  singleParent = matings$singleParent
  sideNode = if(!nrow(matings)) integer(0) else
    unique(c(matings$faNode[singleParent], matings$moNode[singleParent]))

  hasLab = nzchar(labsNode)
  side = hasLab & nodes$node %in% sideNode[sideNode > 0L]
  below = hasLab & !side

  curve = .curveRange(align)
  xrange = range(nodes$x)
  yrange = range(nodes$y)
  if(length(curve$x)) {
    xrange[1L] = min(xrange[1L], curve$x[1L])
    xrange[2L] = max(xrange[2L], curve$x[2L])
    yrange[1L] = min(yrange[1L], curve$y[1L])
    yrange[2L] = max(yrange[2L], curve$y[2L])
  }

  eps = sqrt(.Machine$double.eps)
  xspan = xrange[2L] - xrange[1L]
  yspan = yrange[2L] - yrange[1L]
  if(xspan < eps) {
    xrange = mean(xrange) + c(-0.5, 0.5)
    xspan = 1
  }
  if(yspan < eps) {
    yrange = mean(yrange) + c(-0.5, 0.5)
    yspan = 1
  }

  margin = .marginSpec(margins, marginsUnit)

  # Titles use margin space, not the inner plot region
  if(length(title) && nzchar(title)) {
    titleArgs = .setTitleArgs(titleArgs, cex)
    if(margin$type == "mar")
      margin$value[3L] = max(margin$value[3L], titleArgs$line + titleArgs$cex)
    else
      margin$value[3L] = max(margin$value[3L],
                             (titleArgs$line + titleArgs$cex) * graphics::par("csi"))
  }

  oldpar = graphics::par(c(margin$type, "xpd"))
  on.exit(graphics::par(oldpar), add = TRUE)
  if(margin$type == "mar")
    graphics::par(mar = margin$value, xpd = TRUE)
  else
    graphics::par(mai = margin$value, xpd = TRUE)
  psize = graphics::par("pin")
  if(any(psize <= 0))
    stop2("Margins leave no room for the plot")

  padCoef = .paddingCoef(padding, paddingUnit)
  fixedPad = padCoef$fixed
  symPad = padCoef$sym
  labGapIn = 0.7 * labhIn
  bottomGapIn = 0.3 * labhIn
  levelGapIn = labhIn

  ymin = min(nodes$y)
  nLev = max(nodes$y) - ymin + 1

  if(symbolsizeUnit == "auto") {
    between = below & nodes$y > ymin
    bottom = below & nodes$y == ymin
    maxBetween = if(any(between)) max(labSize$height[between]) else 0
    maxBottom = if(any(bottom)) max(labSize$height[bottom]) else 0
    bottomNeed = if(maxBottom > 0) labGapIn + maxBottom + bottomGapIn else 0
    betweenNeed = if(nLev > 1L) labGapIn + maxBetween + levelGapIn else 0

    vpadIn = fixedPad[1L] + fixedPad[3L]
    vpadSym = symPad[1L] + symPad[3L]
    ht1 = (psize[2L] - vpadIn - bottomNeed - betweenNeed * (nLev - 1L)) /
      (nLev + vpadSym)
    ht2 = (psize[2L] - vpadIn) / (nLev + (nLev - 1L) / 2 + vpadSym)

    wd1 = graphics::strwidth("ABC", units = "inches", cex = cex) * 2.5 / 3
    wd2 = (psize[1L] - fixedPad[2L] - fixedPad[4L]) * 0.8 /
      (0.8 + xspan + 0.8 * (symPad[2L] + symPad[4L]))
    boxIn = symbolsize * min(ht1, ht2, wd1, wd2)
  }
  else {
    boxIn = .unitIn(symbolsize, symbolsizeUnit)
  }

  if(!is.finite(boxIn) || boxIn <= 0)
    stop2("Cannot fit the graph; increase plot region or reduce size")
  boxIn = unname(boxIn)

  paddingIn = fixedPad + symPad * boxIn
  names(paddingIn) = c("bottom", "left", "top", "right")

  # Reserve room for lower-left proband arrows
  probandNode = nodes$id %in% proband

  # Measure annotation boxes in physical units
  annWidth = annHeight = annX = annY = numeric()
  if(nrow(textAnnot)) {
    # Group labels sharing non-vectorised measurement arguments
    group = as.integer(interaction(textAnnot$cex, textAnnot$font,
                                   textAnnot$family, drop = TRUE))
    width = height = numeric(nrow(textAnnot))
    for(g in unique.default(group)) {
      i = which(group == g)

      # Read shared arguments from the first row in the group
      j = i[1L]
      width[i] = graphics::strwidth(
        textAnnot$label[i], units = "inches", cex = textAnnot$cex[j],
        font = textAnnot$font[j], family = textAnnot$family[j])
      height[i] = graphics::strheight(
        textAnnot$label[i], units = "inches", cex = textAnnot$cex[j],
        font = textAnnot$font[j], family = textAnnot$family[j])
    }

    # Rotate the label centre and its bounding box
    angle = textAnnot$srt * pi / 180
    cosa = cos(angle)
    sina = sin(angle)
    midX = (0.5 - textAnnot$adjx) * width
    midY = (0.5 - textAnnot$adjy) * height
    annWidth = abs(cosa) * width + abs(sina) * height
    annHeight = abs(sina) * width + abs(cosa) * height

    # Convert grid positions to physical offsets from node tops
    annGap = textAnnot$offset * labhIn * textAnnot$cex / cex
    annX = textAnnot$hpos * (boxIn / 2 + annGap) +
      cosa * midX - sina * midY
    annY = -(1 - textAnnot$vpos) * boxIn / 2 +
      (textAnnot$hpos == 0) * textAnnot$vpos * annGap +
      sina * midX + cosa * midY
  }

  # Physical extents around node anchors
  xpos = c(nodes$x, nodes$x[below], nodes$x[side], textAnnot$x)

  # Add horizontal overhang from symbols, labels and annotations
  leftIn = c((0.5 + 1.7 * probandNode) * boxIn, labSize$width[below] / 2,
             rep.int(0, sum(side)),
             pmax.int(0, annWidth / 2 - annX))
  rightIn = c(rep.int(boxIn / 2, nrow(nodes)), labSize$width[below] / 2,
              0.75 * boxIn + labSize$width[side],
              pmax.int(0, annX + annWidth / 2))

  ypos = c(nodes$y, nodes$y[below], nodes$y[side], textAnnot$y)

  # Add vertical overhang from symbols, labels and annotations
  lowerIn = c((1 + 0.9 * probandNode) * boxIn,
              boxIn + labGapIn + labSize$height[below] + bottomGapIn,
              boxIn / 2 + labSize$height[side] / 2,
              pmax.int(0, annHeight / 2 - annY))
  upperIn = c(rep.int(0, nrow(nodes)), rep.int(0, sum(below)),
              pmax.int(0, labSize$height[side] / 2 - boxIn / 2),
              pmax.int(0, annY + annHeight / 2))

  if(is.null(xlim)) {
    xs = .axisScale(xpos, xrange, leftIn, rightIn, psize[1L],
                    paddingIn[2L], paddingIn[4L])
    left = xs$limits[1L]
    right = xs$limits[2L]
    hscale = xs$scale
  }
  else {
    left = xlim[1L]
    right = xlim[2L]
    if(right - left <= 0)
      stop2("`xlim` must be increasing")
    hscale = psize[1L] / (right - left)
  }

  if(is.null(ylim)) {
    ys = .axisScale(ypos, yrange, lowerIn, upperIn, psize[2L],
                    paddingIn[1L], paddingIn[3L])
    bottom = ys$limits[1L]
    top = ys$limits[2L]
    vscale = ys$scale
  }
  else {
    bottom = min(ylim)
    top = max(ylim)
    if(top - bottom <= 0)
      stop2("`ylim` must span a positive range")
    vscale = psize[2L] / (top - bottom)
  }

  usr = c(left = left, right = right, bottom = bottom, top = top)

  aboveBottom = below & nodes$y > ymin
  maxBetween = if(any(aboveBottom)) max(labSize$height[aboveBottom]) else 0

  genNeed = boxIn + labGapIn + maxBetween + levelGapIn
  if(nLev > 1L && vscale < genNeed - eps)
    stop2("Labels leave no room between generations; reduce `cex` or symbol size")

  boxw = boxIn / hscale
  boxh = boxIn / vscale
  labh = labhIn / vscale
  legh = min(0.25, boxh * 1.5)

  list(boxw = boxw, boxh = boxh, boxIn = boxIn, labh = labh, labhIn = labhIn,
       legh = legh, cex = cex, margins = margin$value, marginType = margin$type,
       paddingIn = paddingIn, hscale = hscale, vscale = vscale, usr = usr)
}

# Convert a scalar or four-vector to bottom, left, top, right
.recycle4 = function(x, what) {
  if(length(x) == 1L)
    x = rep(x, 4L)
  if(length(x) != 4L)
    stop2("`", what, "` must have length 1 or 4")
  as.numeric(x)
}

# Convert supported physical units to inches
.unitIn = function(x, unit) {
  switch(unit,
         `in` = x,
         cm = x / 2.54,
         mm = x / 25.4,
         stop2("Unknown unit: ", unit))
}

# Resolve margins to par type and values
.marginSpec = function(margins, unit) {
  value = .recycle4(margins, "margins")
  if(unit == "lines")
    return(list(type = "mar", value = value))
  list(type = "mai", value = .unitIn(value, unit))
}

# Split padding into fixed-inch and symbol-size components
.paddingCoef = function(padding, unit) {
  val = .recycle4(padding, "padding")
  if(unit == "sym")
    list(fixed = numeric(4L), sym = val)
  else
    list(fixed = .unitIn(val, unit), sym = numeric(4L))
}

# Measure possibly multiline text in inches
.textSizeIn = function(txt, cex) {
  txt = as.character(txt)
  lineh = graphics::strheight("M", units = "inches", cex = cex)
  width = height = numeric(length(txt))
  nz = nzchar(txt)
  if(!any(nz))
    return(list(width = width, height = height, lineh = lineh))

  # Measure unique labels; split only labels with line breaks
  u = unique(txt[nz])
  uw = graphics::strwidth(u, units = "inches", cex = cex)
  multi = grepl("\n", u, fixed = TRUE)
  if(any(multi)) {
    z = strsplit(u[multi], "\n", fixed = TRUE)
    nline = lengths(z)
    zz = unlist(z, use.names = FALSE)
    w = graphics::strwidth(zz, units = "inches", cex = cex)
    end = cumsum(nline)
    start = end - nline + 1L
    uw[multi] = vapply(seq_along(z), function(i) max(w[start[i]:end[i]]),
                        numeric(1))
  }

  uh = graphics::strheight(u, units = "inches", cex = cex)
  m = match(txt[nz], u)
  width[nz] = uw[m]
  height[nz] = uh[m]
  list(width = width, height = height, lineh = lineh)
}

# Compute the user-coordinate range of duplicate curves
.curveRange = function(align) {
  dup = align$edges$type == "duplicate"
  if(!any(dup))
    return(list(x = numeric(), y = numeric()))

  nodes = align$nodes
  nodeRow = integer(max(nodes$node))
  nodeRow[nodes$node] = seq_len(nrow(nodes))
  fr = nodeRow[align$edges$from[dup]]
  to = nodeRow[align$edges$to[dup]]
  keep = fr > 0L & to > 0L
  if(!any(keep))
    return(list(x = numeric(), y = numeric()))

  x0 = nodes$x[fr[keep]]
  x1 = nodes$x[to[keep]]
  y0 = nodes$y[fr[keep]]
  y1 = nodes$y[to[keep]]

  # Use the analytical range of the parabolic duplicate curve
  dy = y1 - y0
  t = pmin(1, pmax(0, (dy + 2) / 4))
  ymax = y0 + dy * t + 2 * t - 2 * t^2

  list(x = c(min(x0, x1), max(x0, x1)),
       y = c(min(y0, y1), max(ymax)))
}

# Solve one axis scale from user positions and physical extents
.axisScale = function(pos, base, lowerExt, upperExt, psize, lowerPad, upperPad) {
  avail = psize - lowerPad - upperPad
  if(avail <= 0)
    stop2("Padding leaves no room for the graph")

  span = base[2L] - base[1L]
  if(span <= 0)
    stop2("Plot range must be positive")

  distLower = pos - base[1L]
  distUpper = pos - base[2L]

  # Iterate because physical extents interact with user-coordinate scaling
  scale = avail / span
  for(i in seq_len(8L)) {
    lowerIn = max(0, lowerExt - scale * distLower)
    upperIn = max(0, upperExt + scale * distUpper)
    scale = (avail - lowerIn - upperIn) / span
    if(!is.finite(scale) || scale <= 0)
      stop2("Cannot fit the graph; increase plot region or reduce size")
  }

  # Recompute extents at the final scale before returning limits
  lowerIn = max(0, lowerExt - scale * distLower)
  upperIn = max(0, upperExt + scale * distUpper)

  from = base[1L] - (lowerIn + lowerPad) / scale
  to = base[2L] + (upperIn + upperPad) / scale
  list(scale = unname(scale),
       limits = unname(c(from, to)))
}
