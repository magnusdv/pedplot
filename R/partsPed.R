#' Build drawing parts
#'
#' Convert alignment and scale into simple drawing tables.
#'
#' @param align Output from [alignPed()].
#' @param scale Output from [scalePed()].
#' @param labs Labels printed below symbols.
#' @param textAnnot A named list of text annotations at `top`, `left`, `right`,
#'   `bottom`, the four corners or `inside`. Each entry is a label specification,
#'   or a list with labels first and graphical arguments after.
#' @param title Optional title printed in the top margin.
#' @param titleArgs A list of arguments passed to [graphics::mtext()] for the
#'   title.
#' @param col,fill,lty,lwd Graphical values recycled over individuals. Named
#'   values apply to every corresponding ID.
#' @param hatched IDs of individuals whose symbols are hatched.
#' @param hatchCol Hatch colours, defaulting to the corresponding border colours.
#' @param hatchDensity,hatchAngle Hatch density and angle.
#' @param hatchLwd Hatch line width.
#' @param aff IDs of affected individuals, filled black unless overridden by
#'   `fill`.
#' @param carrier IDs of carrier individuals.
#' @param deceased IDs of deceased individuals.
#' @param miscarriage IDs of miscarriage individuals.
#' @param proband IDs of proband individuals.
#' @param ... Not used.
#'
#' @return A list with symbols, segments, curves and texts.
#' @export
#'
#' @examples
#' x = data.frame(id = 1:3, fid = c(0, 0, 1), mid = c(0, 0, 2), sex = c(1, 2, 1))
#' align = alignPed(x)
#' scale = scalePed(align)
#' parts = partsPed(align, scale)
#' names(parts)
partsPed = function(align, scale, labs = align$input$id,
                    col = 1, fill = NA, lty = 1, lwd = 1,
                    hatched = NULL, hatchCol = col, hatchDensity = 25, hatchAngle = 45,
                    hatchLwd = 1, aff = NULL, carrier = NULL, deceased = NULL,
                    miscarriage = NULL, proband = NULL,
                    textAnnot = NULL, title = NULL, titleArgs = list(), ...) {
  input = align$input
  nodes = align$nodes

  # Resolve symbol values once before indexing duplicate nodes
  col = .prepSymbolArg(col, input, 1)
  lty = .prepSymbolArg(lty, input, 1)
  lwd = .prepSymbolArg(lwd, input, 1)
  hatchCol = .prepSymbolArg(hatchCol, input, col)
  hatched = input$id %in% hatched
  carrier = input$id %in% carrier
  deceased = input$id %in% deceased
  proband = input$id %in% proband

  # Apply affected fill unless overridden explicitly
  fillDefault = rep.int(NA_integer_, nrow(input))
  fillDefault[input$id %in% aff] = 1L
  fill = if(missing(fill)) fillDefault else
    .prepSymbolArg(fill, input, fillDefault)

  # Symbol shape
  shape = .shapeSex(nodes$sex)
  shape[nodes$id %in% miscarriage] = "triangle"

  # Convert node-top y to the lower symbol edge used by drawPed
  symbols = list2DF(list(
    node = nodes$node,
    idx = nodes$idx,
    id = nodes$id,
    x = nodes$x,
    y = nodes$y - scale$boxh,
    width = rep(scale$boxw, nrow(nodes)),
    height = rep(scale$boxh, nrow(nodes)),
    sex = nodes$sex,
    shape = shape,
    col = col[nodes$idx],
    fill = fill[nodes$idx],
    lty = lty[nodes$idx],
    lwd = lwd[nodes$idx],
    hatched = hatched[nodes$idx],
    hatchCol = hatchCol[nodes$idx],
    carrier = carrier[nodes$idx],
    deceased = deceased[nodes$idx],
    proband = proband[nodes$idx]))

  titlePart = NULL
  if(length(title) && nzchar(title)) {
    titlePart = list(text = title[1L],
                     args = .setTitleArgs(titleArgs, cex = scale$cex))
  }

  # Move annotation anchors from node tops to their symbol positions
  textAnnot = .prepTextAnnot(align, textAnnot, scale$cex)
  if(nrow(textAnnot)) {
    gapIn = textAnnot$offset * scale$labhIn * textAnnot$cex / scale$cex

    # Side positions move horizontally; centred top and bottom move vertically
    textAnnot$x = textAnnot$x + textAnnot$hpos * (scale$boxw / 2 + gapIn / scale$hscale)
    textAnnot$y = textAnnot$y - (1 - textAnnot$vpos) * scale$boxh / 2 +
      (textAnnot$hpos == 0) * textAnnot$vpos * gapIn / scale$vscale
  }
  textAnnot = textAnnot[c("node", "idx", "id", "x", "y", "label", "cex",
                          "col", "font", "family", "srt", "adjx", "adjy")]

  list(symbols = symbols,
       segments = .partsSegments(align, scale),
       curves = .partsCurves(align),
       texts = .partsTexts(align, scale, labs),
       textAnnot = textAnnot,
       title = titlePart,
       hatchDensity = hatchDensity[1L],
       hatchAngle = hatchAngle[1L],
       hatchLwd = hatchLwd[1L])
}

# Expand graphical values by input position or ID
.prepSymbolArg = function(x, input, default) {
  n = nrow(input)
  nms = names(x)
  if(is.null(nms))
    return(rep_len(x, n))

  ans = rep_len(default, n)
  m = match(input$id, nms, nomatch = 0L)
  ans[m > 0L] = x[m]
  ans
}


# Build all straight connector segments
.partsSegments = function(align, scale) {
  nodes = align$nodes
  edges = align$edges
  mats = align$matings
  boxw = scale$boxw
  boxh = scale$boxh
  legh = scale$legh
  branch = 0.6
  pconnect = 0.5

  # Preallocate an upper bound and trim at the end
  nMax = 2L * sum(edges$type == "spouse") + nrow(edges) + 4L * nrow(mats)
  type = character(nMax)
  x0 = y0 = x1 = y1 = numeric(nMax)
  col = lty = lwd = numeric(nMax)
  nSeg = 0L

  nodeRow = integer(max(nodes$node))
  nodeRow[nodes$node] = seq_len(nrow(nodes))

  # Horizontal spouse connectors
  sp = edges[edges$type == "spouse", ]
  if(nrow(sp)) {
    fr = nodeRow[sp$from]
    to = nodeRow[sp$to]
    keep = fr > 0L & to > 0L

    if(any(keep)) {
      fr = fr[keep]
      to = to[keep]
      sp = sp[keep, ]
      left = ifelse(nodes$x[fr] <= nodes$x[to], fr, to)
      right = ifelse(nodes$x[fr] <= nodes$x[to], to, fr)
      yy = nodes$y[left] - boxh / 2
      len = length(left)
      ii = nSeg + seq_len(len)

      type[ii] = "spouse"
      x0[ii] = nodes$x[left] + boxw / 2
      y0[ii] = yy
      x1[ii] = nodes$x[right] - boxw / 2
      y1[ii] = yy
      col[ii] = lty[ii] = lwd[ii] = 1
      nSeg = nSeg + len

      # Consanguineous matings are drawn as a second spouse line
      double = sp$code == "2"
      if(any(double)) {
        dd = which(double)
        len = length(dd)
        ii = nSeg + seq_len(len)
        yy2 = yy[dd] - boxh / 10

        type[ii] = "spouse2"
        x0[ii] = nodes$x[left[dd]] + boxw / 2
        y0[ii] = yy2
        x1[ii] = nodes$x[right[dd]] - boxw / 2
        y1[ii] = yy2
        col[ii] = lty[ii] = lwd[ii] = 1
        nSeg = nSeg + len
      }
    }
  }

  child = edges[edges$type == "child", ]
  if(nrow(child) && nrow(mats)) {
    matRows = match(unique(child$mating), mats$mating)

    # Loop over mating groups and add child bars plus parent stems
    for(mm in matRows) {
      m = mats$mating[mm]
      kidEdges = child[child$mating == m, ]
      ch = nodeRow[kidEdges$to]
      ch = ch[ch > 0L]
      if(!length(ch))
        next

      target = nodes$x[ch]
      yTop = max(nodes$y[ch]) + legh
      len = length(ch)
      ii = nSeg + seq_len(len)

      # Each child gets one short vertical line up to the sibship bar
      type[ii] = "child"
      x0[ii] = nodes$x[ch]
      y0[ii] = nodes$y[ch]
      x1[ii] = target
      y1[ii] = yTop
      col[ii] = lty[ii] = lwd[ii] = 1
      nSeg = nSeg + len

      if(len > 1L) {
        nSeg = nSeg + 1L
        type[nSeg] = "sib"
        x0[nSeg] = min(target)
        x1[nSeg] = max(target)
        y0[nSeg] = y1[nSeg] = yTop
        col[nSeg] = lty[nSeg] = lwd[nSeg] = 1
      }

      parNodes = c(mats$faNode[mm], mats$moNode[mm])
      parNodes = parNodes[parNodes > 0L]
      if(!length(parNodes))
        next

      par = nodeRow[parNodes]
      par = par[par > 0L]
      if(!length(par))
        next

      px = mean(nodes$x[par])
      py = if(length(par) == 1L) nodes$y[par] - boxh else mean(nodes$y[par] - boxh / 2)
      trng = range(target)

      # Keep the parent stem away from the ends of wide sibship bars
      xA = if(trng[2L] - trng[1L] < 2 * pconnect)
        mean(trng)
      else
        pmax.int(trng[1L] + pconnect, pmin.int(trng[2L] - pconnect, px))

      xB = px
      yA = yTop
      yB = py
      yDel = ((yB - yA) * branch) / 2
      ii = nSeg + seq_len(3L)

      # Stem: vertical, slanted, vertical
      type[ii] = "parent"
      x0[ii] = c(xA, xA, xB)
      y0[ii] = c(yA, yA + yDel, yB - yDel)
      x1[ii] = c(xA, xB, xB)
      y1[ii] = c(yA + yDel, yB - yDel, yB)
      col[ii] = lty[ii] = lwd[ii] = 1
      nSeg = nSeg + 3L
    }
  }

  nsq = seq_len(nSeg)
  list2DF(list(type = type[nsq],
                x0 = x0[nsq],
                y0 = y0[nsq],
                x1 = x1[nsq],
                y1 = y1[nsq],
                col = col[nsq],
                lty = lty[nsq],
                lwd = lwd[nsq]))
}


# Build dotted duplicate curves as sampled points
.partsCurves = function(align) {
  edges = align$edges[align$edges$type == "duplicate", ]
  nodes = align$nodes
  if(!nrow(edges))
    return(list2DF(list(type = character(), curve = integer(), x = numeric(),
                      y = numeric(), id = character(), col = numeric(),
                      lty = numeric(), lwd = numeric())))

  nodeRow = integer(max(nodes$node))
  nodeRow[nodes$node] = seq_len(nrow(nodes))
  fr = nodeRow[edges$from]
  to = nodeRow[edges$to]
  keep = fr > 0L & to > 0L
  if(!any(keep))
    return(list2DF(list(type = character(), curve = integer(), x = numeric(),
                      y = numeric(), id = character(), col = numeric(),
                      lty = numeric(), lwd = numeric())))

  fr = fr[keep]
  to = to[keep]
  nPair = length(fr)
  nPts = 15L
  nTot = nPair * nPts
  curve = rep(seq_len(nPair), each = nPts)
  x = y = numeric(nTot)
  id = rep(nodes$id[fr], each = nPts)
  bend = 0.5 - (seq(-7, 7))^2 / 98

  # Loop over duplicate pairs and write one point block per curve
  for(i in seq_len(nPair)) {
    ii = (i - 1L) * nPts + seq_len(nPts)
    x[ii] = seq.int(nodes$x[fr[i]], nodes$x[to[i]], length.out = nPts)
    y[ii] = seq.int(nodes$y[fr[i]], nodes$y[to[i]], length.out = nPts) + bend
  }

  list2DF(list(
    type = rep("duplicate", nTot),
    curve = curve, x = x, y = y,
    id = id, col = rep(1, nTot),
    lty = rep(2, nTot),
    lwd = rep(1, nTot)))
}

# Build below and side label positions from visible nodes
.partsTexts = function(align, scale, labs) {
  nodes = align$nodes
  labs = .prepLabs(align$input, labs)

  useBelow = which(nzchar(labs[nodes$idx]))
  if(!length(useBelow))
    return(list2DF(list(type = character(), node = integer(), idx = integer(),
                       id = character(), x = numeric(), y = numeric(),
                       label = character(), cex = numeric(), col = numeric(),
                       adjx = numeric(), adjy = numeric())))

  hiddenSpouse = .singleParentNodes(align)
  side = nodes$node[useBelow] %in% hiddenSpouse
  n = length(useBelow)
  type = rep.int("label", n)
  x = nodes$x[useBelow]
  y = nodes$y[useBelow] - scale$boxh - scale$labh * 0.7
  adjx = rep.int(0.5, n)
  adjy = rep.int(1, n)

  # Side labels leave room for child connectors from single parents
  if(any(side)) {
    x[side] = nodes$x[useBelow[side]] + 0.75 * scale$boxw
    y[side] = nodes$y[useBelow[side]] - scale$boxh / 2
    adjx[side] = 0
    adjy[side] = 0.5
    type[side] = "labelSide"
  }

  list2DF(list(
    type = type, node = nodes$node[useBelow], idx = nodes$idx[useBelow],
    id = nodes$id[useBelow], x = x, y = y,
    label = labs[nodes$idx[useBelow]], cex = rep(scale$cex, n),
    col = rep(1, n), adjx = adjx, adjy = adjy))
}


# Prepare text annotations in visible-node order
.prepTextAnnot = function(align, textAnnot, cex) {
  # Encode left/centre/right and bottom/middle/top as -1/0/1
  positions = c("topleft", "top", "topright", "left", "inside", "right",
                "bottomleft", "bottom", "bottomright")
  hpos = rep.int(-1:1, 3L)
  vpos = c(1, 1, 1, 0, 0, 0, -1, -1, -1)

  # Empty annotation table
  empty = list2DF(list(
    node = integer(), idx = integer(), id = character(), x = numeric(),
    y = numeric(), label = character(), cex = numeric(), col = character(),
    font = numeric(), family = character(), srt = numeric(), adjx = numeric(),
    adjy = numeric(), offset = numeric(), hpos = numeric(), vpos = numeric()))
  if(!length(textAnnot))
    return(empty)

  # Match list names to the supported grid positions
  if(is.null(names(textAnnot)))
    stop2("`textAnnot` must be named")
  pos = match(names(textAnnot), positions)
  if(anyNA(pos))
    stop2("Unknown `textAnnot` position: ", names(textAnnot)[is.na(pos)][1L])

  # Use visible nodes so plot duplicates receive the same annotation
  nodes = align$nodes

  # Build textAnnot table by rows
  rows = lapply(seq_along(textAnnot), function(i) {
    a = textAnnot[[i]]
    if(is.atomic(a) || is.function(a))
      a = list(a)

    # Prepare labels and discard empty rows
    label = .prepLabs(align$input, a[[1L]])[nodes$idx]
    use = which(nzchar(label))
    if(!length(use))
      return(NULL)

    # Recycle graphical values in input order before selecting nodes
    p = pos[i]
    n = length(use)
    pick = function(x, default)
      rep_len(x %||% default, nrow(align$input))[nodes$idx[use]]

    # Convert the grid position to text adjustment at its anchor
    hp = hpos[p]
    vp = vpos[p]
    adj0 = c((1 - hp) / 2, if(hp == 0) (1 - vp) / 2 else 0.5)
    adj = c(a$adj %||% adj0, 0.5)[1:2]

    # Centre missing adjustment components
    adj[is.na(adj)] = 0.5

    # Store node coordinates. NB: partsPed moves the anchors after scaling
    list2DF(list(
      node = nodes$node[use], idx = nodes$idx[use], id = nodes$id[use],
      x = nodes$x[use], y = nodes$y[use], label = label[use],
      cex = pick(a$cex, cex), col = pick(a$col, 1),
      font = pick(a$font, 1), family = as.character(pick(a$family, "")),
      srt = pick(a$srt, 0), adjx = rep(adj[1L], n),
      adjy = rep(adj[2L], n), offset = pick(a$offset, 0.5),
      hpos = rep(hp, n), vpos = rep(vp, n)))
  })

  # Bind position tables once
  rows = rows[lengths(rows) > 0L]
  if(length(rows)) .myRbind(rows) else empty
}


# Find visible parent nodes in single-parent matings
.singleParentNodes = function(align) {
  mats = align$matings
  if(!nrow(mats) || !"singleParent" %in% names(mats))
    return(integer())

  single = mats$singleParent
  if(!any(single))
    return(integer())

  # Hidden founder nodes are already removed, so the non-zero parent is visible
  ans = c(mats$faNode[single], mats$moNode[single])
  unique.default(ans[ans > 0L])
}
