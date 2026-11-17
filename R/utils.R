`%||%` = function(x, y) if(is.null(x)) y else x

stop2 = function(...) {
    a = lapply(list(...), toString)
    a = append(a, list(call. = FALSE))
    do.call(stop, a)
}

.mysetdiff = function (x, y)
  unique.default(x[match(x, y, 0L) == 0L])

.myintersect = function (x, y)
  y[match(x, y, 0L)]

.parentIdx = function(par, id, what) {
  idx = match(par, id, nomatch = 0L)
  bad = par != "0" & idx == 0L
  if(any(bad))
    stop2("Unknown ", what, " ID: ", par[bad][1])
  idx
}

.sexName = function(sex) {
  ans = rep("unknown", length(sex))
  ans[sex == 1L] = "male"
  ans[sex == 2L] = "female"
  ans
}

.shapeSex = function(sex) {
  ans = rep("diamond", length(sex))
  ans[sex == 1L] = "square"
  ans[sex == 2L] = "circle"
  ans[sex == 3L] = "triangle"
  ans
}

.checkHidden = function(input, hideFounders) {
  if(is.null(hideFounders))
    return(input)

  hid = unique.default(as.character(hideFounders))
  idx = match(hid, input$id)

  bad = is.na(idx)
  if(any(bad))
    message("Ignoring unknown `hideFounders`: ", toString(hid[bad]))
  idx = idx[!bad]
  if(!length(idx))
    return(input)

  founder = input$fid == 0 & input$mid == 0
  bad = !founder[idx]
  if(any(bad))
    message("Ignoring non-founder `hideFounders`: ",
            toString(input$id[idx[bad]]))
  idx = idx[!bad]
  if(!length(idx))
    return(input)

  # Restrict hiding to simple monogamous founder matings
  fa = .parentIdx(input$fid, input$id, "father")
  mo = .parentIdx(input$mid, input$id, "mother")
  i = which(fa > 0L | mo > 0L)
  i = i[!duplicated.default(fa[i] * (length(fa) + 1) + mo[i])]
  pair = cbind(fa[i], mo[i])
  nMate = tabulate(pair[pair > 0L], nbins = nrow(input))

  messTmp = "Ignoring `hideFounders = %s`: %s"

  # Process requests in order so valid founders are still hidden
  for(i in idx) {
    mat = pair[, 1L] == i | pair[, 2L] == i
    near = unique.default(c(pair[mat, ]))
    spouse = near[near > 0L & near != i]

    if(any(input$hidden[spouse])) {
      message(sprintf(messTmp, input$id[i], "both parents cannot be hidden"))
      next
    }
    if(any(nMate[near[near > 0L]] > 1L)) {
      message(sprintf(messTmp, input$id[i], "non-monogamous founder or spouse"))
      next
    }
    if(sum(!input$hidden) == 1L) {
      message(sprintf(messTmp, input$id[i], "cannot hide all individuals"))
      next
    }
    input$hidden[i] = TRUE
  }

  input
}


.prepLabs = function(input, labs) {
  if(is.null(labs))
    return(rep("", nrow(input)))

  if(is.function(labs))
    labs = labs(input)

  if(length(labs) == 1L && identical(labs, "num"))
    return(as.character(input$idx))

  mode(labs) = "character"
  nms = names(labs)

  if(is.null(nms) && length(labs) == nrow(input))
    return(labs)

  ans = rep("", nrow(input))
  mtch = match(input$id, nms %||% labs, nomatch = 0L)
  ok = mtch > 0L
  ans[ok] = labs[mtch[ok]]

  # Hide hidden (works with NULL and empty)
  ans[input$hidden] = ""

  ans
}

# Prepare graphical arguments for title
.setTitleArgs = function(titleArgs, cex) {
  if(!is.list(titleArgs))
    stop2("`titleArgs` must be a list")

  args = utils::modifyList(list(cex = max(1.2 * cex, 1), font = 2, col = 1, line = 0.5),
                    titleArgs)
  args[c("text", "side")] = NULL
  args
}


.myRbind = function(dfs) {
  nms = colnames(dfs[[1L]])
  lst = lapply(nms, function(nm) unlist(lapply(dfs, `[[`, nm), use.names = FALSE))
  names(lst) = nms
  list2DF(lst)
}
