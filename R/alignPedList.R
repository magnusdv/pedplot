# Align and arrange pedigree components before merging their tables
.alignPedList = function(x, hideFounders, gap = 1,
                         nrow = NULL, ncol = NULL, widths = NULL, heights = NULL,
                         ...) {
  if(!length(x))
    stop2("No pedigrees")

  inputs = lapply(x, .getInput)
  nInput = vapply(inputs, NROW, integer(1L))
  if(any(!nInput))
    stop2("No individuals")

  input = .myRbind(inputs)
  input$idx = seq_len(NROW(input))

  # Apply hidden-founder requests separately within each component
  hide = unique.default(as.character(hideFounders))
  if(length(hide)) {
    unknown = hide[!hide %in% input$id]
    if(length(unknown))
      message("Ignoring unknown `hideFounders`: ", toString(unknown))
    hide = hide[hide %in% input$id]

    inputs = lapply(inputs, function(z) {
      h = .myintersect(hide, z$id)
      if(length(h)) .checkHidden(z, h) else z
    })
  }

  input$hidden = unlist(lapply(inputs, `[[`, "hidden"), use.names = FALSE)
  aligns = lapply(inputs, .alignBlocks)
  nComp = length(aligns)

  # Default to 1 row
  if(is.null(nrow) && is.null(ncol)) {
    nrow = 1L
    ncol = nComp
  }
  else if(is.null(nrow))
    nrow = ceiling(nComp / ncol)
  else if(is.null(ncol))
    ncol = ceiling(nComp / nrow)
  nrow = as.integer(nrow)
  ncol = as.integer(ncol)
  if(nrow < 1L || ncol < 1L || nrow * ncol < nComp)
    stop2("Grid is too small")

  comp = seq_len(nComp)
  gridRow = (comp - 1L) %/% ncol + 1L
  gridCol = (comp - 1L) %% ncol + 1L

  # Find the space required by each component
  xr = vapply(aligns, function(a) range.default(a$nodes$x), numeric(2L))
  yr = vapply(aligns, function(a) range.default(a$nodes$y), numeric(2L))
  compWidth = xr[2L, ] - xr[1L, ]
  compHeight = yr[2L, ] - yr[1L, ]
  colWidth = vapply(seq_len(ncol), function(j)
    max(compWidth[gridCol == j], 0), numeric(1L))
  rowHeight = vapply(seq_len(nrow), function(j)
    max(compHeight[gridRow == j], 0), numeric(1L))

  # Enlarge cells to obtain the requested relative dimensions
  if(length(widths)) {
    widths = rep_len(widths, ncol)
    if(any(widths <= 0))
      stop2("`widths` must be positive")
    colWidth = widths * max(pmax.int(colWidth, 1) / widths)
  }
  if(length(heights)) {
    heights = rep_len(heights, nrow)
    if(any(heights <= 0))
      stop2("`heights` must be positive")
    rowHeight = heights * max(pmax.int(rowHeight, 1) / heights)
  }

  # Use separate horizontal and vertical gaps when supplied
  gap = rep_len(gap, 2L)

  # Centre components in their cells and rows from top to bottom
  colLeft = c(0, cumsum(colWidth[-ncol] + gap[1L]))
  rowTop = -c(0, cumsum(rowHeight[-nrow] + gap[2L]))

  xShift = colLeft[gridCol] + colWidth[gridCol] / 2 -
    (xr[1L, ] + xr[2L, ]) / 2
  yShift = rowTop[gridRow] - rowHeight[gridRow] / 2 -
    (yr[1L, ] + yr[2L, ]) / 2

  # Centre horizontally and place the completed alignment above zero
  xShift = xShift - mean(range.default(xr[1L, ] + xShift,
                                       xr[2L, ] + xShift))
  yShift = yShift + 1 - min(yr[1L, ] + yShift, yr[2L, ] + yShift)

  # Offset all local row references before binding each table once
  nNode = vapply(aligns, function(a) NROW(a$nodes), integer(1L))
  nMating = vapply(aligns, function(a) NROW(a$matings), integer(1L))
  idxOff = c(0L, cumsum(nInput[-length(nInput)]))
  nodeOff = c(0L, cumsum(nNode[-length(nNode)]))
  matingOff = c(0L, cumsum(nMating[-length(nMating)]))

  aligns = lapply(seq_along(aligns), function(i) {
    a = aligns[[i]]
    a$nodes$node = a$nodes$node + nodeOff[i]
    a$nodes$idx = a$nodes$idx + idxOff[i]
    a$nodes$x = a$nodes$x + xShift[i]
    a$nodes$y = a$nodes$y + yShift[i]

    a$matings$mating = a$matings$mating + matingOff[i]
    a$matings$faIdx = a$matings$faIdx + idxOff[i] * (a$matings$faIdx > 0L)
    a$matings$moIdx = a$matings$moIdx + idxOff[i] * (a$matings$moIdx > 0L)
    a$matings$faNode = a$matings$faNode + nodeOff[i] * (a$matings$faNode > 0L)
    a$matings$moNode = a$matings$moNode + nodeOff[i] * (a$matings$moNode > 0L)

    # Child edges store mating numbers in `from`; other edges store node numbers
    child = a$edges$type == "child"
    a$edges$from[child] = a$edges$from[child] + matingOff[i]
    a$edges$from[!child] = a$edges$from[!child] + nodeOff[i]
    a$edges$to = a$edges$to + nodeOff[i]
    a$edges$mating = a$edges$mating + matingOff[i]
    a$edges$child = a$edges$child + idxOff[i]
    a
  })

  nodes = .myRbind(lapply(aligns, `[[`, "nodes"))
  matings = .myRbind(lapply(aligns, `[[`, "matings"))
  edges = .myRbind(lapply(aligns, `[[`, "edges"))
  edges$edge = seq_len(NROW(edges))

  ids = input$id
  align = list(
    input = input, nodes = nodes, matings = matings, edges = edges,
    meta = list(hidden = ids[input$hidden],
                duplicated = unique.default(
                  ids[nodes$idx[duplicated.default(nodes$idx)]]),
                maxlev = max(nodes$y)))

  class(align) = "pedalign"
  align
}
