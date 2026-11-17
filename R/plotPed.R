#' Plot pedigrees
#'
#' Plot a pedigree or a list of pedigrees.
#'
#' `plotPed()` is the main plotting function in **pedplot**. The input may be a
#' pedigree data frame, a `pedtools::ped` object, or a list of either.
#'
#' The workflow is approximately as follows:
#'
#' ```
#' align = alignPed(x, ...)
#' scale = scalePed(align, labs = labs, ...)
#' parts = partsPed(align, scale, labs = labs, ...)
#'
#' if(draw)
#'   drawPed(parts, scale)
#' ```
#'
#' Most arguments are passed through `...`. Arguments controlling pedigree
#' alignment are handled by [alignPed()], sizes and plot limits by [scalePed()],
#' and graphical content by [partsPed()].
#'
#' Graphical vectors may be supplied in input order or as named vectors, in
#' which case the names are interpreted as individual IDs. The `textAnnot`
#' argument gives text around or inside symbols; see the examples below.
#'
#' With `draw = FALSE`, the alignment, scaling and drawing parts are returned
#' without producing a plot. The prepared plot may subsequently be drawn with
#' `drawPed(res$parts, res$scale)`. The `scale` object is required separately
#' because `parts` contains drawing coordinates, while `scale` contains the
#' margins and plot limits needed to interpret those coordinates.
#'
#' @param x A pedigree data frame with columns `id`, `fid`, `mid` and `sex`, a
#'   `pedtools::ped` object, or a list of either.
#' @param labs Labels printed below pedigree symbols. By default, individual IDs
#'   are used. A named vector applies labels by ID, while an unnamed vector is
#'   interpreted in input order.
#' @param draw A logical indicating whether the plot should be drawn.
#' @inheritDotParams alignPed hideFounders
#' @inheritDotParams scalePed cex symbolsize symbolsizeUnit margins marginsUnit padding paddingUnit xlim ylim
#' @inheritDotParams partsPed col fill lty lwd hatched hatchCol hatchDensity hatchAngle hatchLwd aff carrier deceased miscarriage proband textAnnot title titleArgs
#'
#' @returns Invisibly, an object containing the alignment (`align`), scaling
#'   information (`scale`), drawing tables (`parts`) and matched call (`call`).
#'
#' @seealso [alignPed()], [scalePed()] and [partsPed()] for the individual
#'   preparation steps, and [drawPed()] for drawing prepared parts.
#'
#' @examples
#' # Trio pedigree
#' x = data.frame(id = 1:3, fid = c(0, 0, 1),
#'                mid = c(0, 0, 2), sex = c(1, 2, 1))
#' plotPed(x)
#'
#' # Different colours, line types and widths
#' plotPed(x, col = 2, fill = "grey90", lwd = 3)
#' plotPed(x, col = 2:4, lty = 1:3, cex = 1.5, lwd = 3)
#'
#' plotPed(x, cex = 2, fill = c("2" = "gold"))
#' plotPed(x, cex = 2, fill = c("2" = "gold"), col = c("2" = "gold"))
#' plotPed(x, cex = 2, fill = c("2" = "gold"), hatched = "2", hatchDensity = 20)
#'
#' # Hatching
#' plotPed(x, hatched = 3, cex = 1.5)
#' plotPed(x, hatched = 3, lwd = 2, cex = 1.5)
#' plotPed(x, hatched = 1:3, hatchCol = 2:4, hatchAngle = 0, hatchLwd = 2)
#'
#' # Standard pedigree annotations
#' plotPed(x, carrier = 1:2, deceased = 2, aff = 3)
#' plotPed(x, miscarriage = 3, deceased = 3)
#' plotPed(x, aff = 3, proband = 3)
#'
#' # Text annotations
#' plotPed(x, textAnnot = list(top = c("1" = "Father"),
#'                             inside = c("3" = "P")))
#' plotPed(x, textAnnot = list(
#'   left = list(c("3" = "Proband"), col = 2, font = 2)))
#'
#' # Multiple pedigrees
#' plotPed(list(x, x, x))
#'
#' # Control layout and gap between multiple pedigrees
#' plotPed(list(x,x,x), nrow = 2, gap = 0.5)
#'
#'
#' # Prepare first and draw later
#' res = plotPed(x, draw = FALSE)
#' drawPed(res$parts, res$scale)
#'
#' # Manual edit of the plot data
#' res$parts$symbols$shape[1] = "diamond"
#' res$parts$symbols$fill[1] = "yellow"
#' drawPed(res$parts, res$scale)
#'
#' @export
plotPed = function(x = NULL, labs = NULL, draw = TRUE, ...) {
  call = match.call()

  align = alignPed(x, ...)

  labs = labs %||% align$input$id

  scale = scalePed(align, labs = labs, ...)
  parts = partsPed(align, scale, labs = labs, ...)

  res = list(align = align, scale = scale, parts = parts, call = call)
  class(res) = "pedplot"

  if(draw)
    drawPed(parts, scale)

  invisible(res)
}


