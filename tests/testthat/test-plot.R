test_that("plot wrapper can return without drawing", {
  x = data.frame(id = 1:3, fid = c(0, 0, 1), mid = c(0, 0, 2), sex = c(1, 2, 1))
  res = plotPed(x, draw = FALSE)
  expect_s3_class(res, "pedplot")
  expect_named(res, c("align", "scale", "parts", "call"))
})

test_that("plot titles are margin parts", {
  x = data.frame(id = 1:3, fid = c(0, 0, 1),
                 mid = c(0, 0, 2), sex = c(1L, 2L, 1L))

  r = plotPed(x, draw = FALSE)
  expect_null(r$parts$title)

  r = plotPed(x, title = "Pedigree A", draw = FALSE)
  expect_equal(r$parts$title$text, "Pedigree A")
  expect_equal(r$parts$title$args$font, 2)
  expect_equal(r$parts$title$args$line, 0.5)

  r = plotPed(x, title = "Pedigree A",
              titleArgs = list(cex = 0.8, col = 2, line = 1),
              draw = FALSE)
  expect_equal(r$parts$title$args$cex, 0.8)
  expect_equal(r$parts$title$args$col, 2)
  expect_equal(r$parts$title$args$line, 1)
})
