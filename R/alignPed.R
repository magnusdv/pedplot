#' Align a pedigree
#'
#' Build plot positions for a pedigree. Input data must use the pedsuite column
#' names `id`, `fid`, `mid` and `sex`.
#'
#' @param x A data frame with columns `id`, `fid`, `mid` and `sex`.
#' @param hideFounders Founder IDs to hide from the plot.
#' @param ... Further arguments.
#'
#' @return A list with input, nodes, matings, edges and metadata.
#' @export
#'
#' @examples
#' x = data.frame(id = 1:4, fid = c(0, 0, 1, 1),
#'                mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
#' align = alignPed(x)
#' align$nodes
alignPed = function(x = NULL, hideFounders = NULL, ...) {

  if(is.list(x) && !is.data.frame(x) && !inherits(x, "ped"))
    return(.alignPedList(x, hideFounders, ...))

  input = x |> .getInput() |> .checkHidden(hideFounders)
  .alignBlocks(input)
}


# Check pedigree input and return standardised data frame
.getInput = function(x) {
  # Convert to data frame (e.g. from `ped`)
  if(!is.data.frame(x))
    x = as.data.frame(x)

  miss = .mysetdiff(c("id", "fid", "mid", "sex"), names(x))
  if(length(miss))
    stop2("Missing required column: ", miss)

  x$id = as.character(x$id)
  if(anyDuplicated.default(x$id))
    stop2("Duplicated ID: ", x$id[duplicated(x$id)])

  x$idx = seq_len(nrow(x))
  x$hidden = FALSE

  x[c("idx", "id", "fid", "mid", "sex", "hidden")]
}
