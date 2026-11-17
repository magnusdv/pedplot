test_that("parts contain finite drawing data", {
  x = data.frame(
    id = c(1:8),
    fid = c(0, 0, 1, 0, 1, 0, 3, 5),
    mid = c(0, 0, 2, 0, 2, 0, 4, 6),
    sex = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L))
  a = alignPed(x)
  s = scalePed(a)
  p = partsPed(a, s)

  expect_true(all(is.finite(p$symbols$width)))
  expect_true(all(is.finite(p$symbols$height)))
  if(nrow(p$segments))
    expect_true(all(is.finite(unlist(p$segments[c("x0", "y0", "x1", "y1")]))))
})


test_that("hidden-founder parts are visible-only", {
  x = data.frame(id = 1:4, fid = c(0, 0, 1, 1),
                 mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
  a = alignPed(x, hideFounders = "1")
  s = scalePed(a)
  p = partsPed(a, s)

  expect_lt(s$usr[3L], s$usr[4L])
  expect_false("1" %in% p$symbols$id)
  expect_false(any(p$segments$type == "spouse"))
  expect_true(any(p$texts$type == "labelSide"))
  expect_true(all(is.finite(unlist(p$symbols[c("x", "y", "width", "height")]))))
  if(nrow(p$segments))
    expect_true(all(is.finite(unlist(p$segments[c("x0", "y0", "x1", "y1")]))))
})

