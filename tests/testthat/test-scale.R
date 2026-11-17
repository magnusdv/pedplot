test_that("scaling returns explicit sizes and conversion factors", {
  x = data.frame(id = 1:3, fid = c(0, 0, 1), mid = c(0, 0, 2), sex = c(1, 2, 1))
  a = alignPed(x)
  s = scalePed(a)

  expect_true(is.finite(s$boxw))
  expect_true(is.finite(s$boxh))
  expect_true(is.finite(s$boxIn))
  expect_null(s$boxUser)
  expect_null(s$paddingUser)
  expect_equal(s$marginType, "mar")
})

test_that("fixed symbol sizes use physical units", {
  x = data.frame(id = 1:3, fid = c(0, 0, 1), mid = c(0, 0, 2), sex = c(1, 2, 1))
  a = alignPed(x)

  s = scalePed(a, symbolsize = 4, symbolsizeUnit = "mm")
   expect_equal(s$boxIn, 4 / 25.4)

  s = scalePed(a, symbolsize = 0.5, symbolsizeUnit = "cm")
  expect_equal(s$boxIn, 0.5 / 2.54)
})

test_that("physical margins and padding are resolved in inches", {
  x = data.frame(id = 1:4, fid = c(0, 0, 1, 1), mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
  a = alignPed(x)

  s1 = scalePed(a, margins = c(1, 2, 1, 2), marginsUnit = "cm")
  expect_equal(s1$marginType, "mai")
  expect_equal(s1$margins, c(1, 2, 1, 2) / 2.54)

  s2 = scalePed(a, symbolsize = 4, symbolsizeUnit = "mm",
               padding = c(1, 2, 3, 4), paddingUnit = "sym")
  expect_equal(unname(s2$paddingIn), c(1, 2, 3, 4) * 4 / 25.4)

})

test_that("scaling keeps non-degenerate limits", {
  x = data.frame(id = 1:7, fid = c(0, 0, 1, 0, 3, 0, 5),
                 mid = c(0, 0, 2, 0, 4, 0, 6), sex = c(1L, 2L, 1L, 2L, 1L, 2L, 1L))
  r = plotPed(x, hideFounders = c(2, 4, 6), draw = FALSE)

  expect_gt(diff(r$scale$usr[1:2]), 0)
  expect_gt(diff(r$scale$usr[3:4]), 0)
})

test_that("long and multiline labels increase horizontal limits", {
  x = data.frame(id = 1:5, fid = c(0, 0, 1, 1, 1), mid = c(0, 0, 2, 2, 2), sex = c(1, 2, 1, 2, 1))
  a = alignPed(x)
  s0 = scalePed(a)
  labs = c("founder 1", "founder 2", "a very very long label", "short", "line 1\nline 2")
  s1 = scalePed(a, labs = labs)

  expect_gt(diff(s1$usr[1:2]), diff(s0$usr[1:2]))
})

test_that("duplicate curves contribute to top limits", {
  x = data.frame(id = 1:8,
                 fid = c(0, 0, 1, 1, 1, 4, 5, 6),
                 mid = c(0, 0, 2, 2, 2, 2, 2, 2),
                 sex = c(1L, 2L, rep(1L, 6L)))
  a = alignPed(x)
  s = scalePed(a)

  expect_true(any(a$edges$type == "duplicate"))
  expect_gt(s$usr[4L], max(a$nodes$y))
})

test_that("multiline bottom labels fit inside plot region", {
  x = data.frame(id = 1:5, fid = c(0, 0, 1, 1, 1),
                 mid = c(0, 0, 2, 2, 2), sex = c(1, 2, 1, 2, 1))
  labs = c("founder 1", "founder 2", "a very long label",
           "short", "very long first line\nline 2\nline 3\nline 4")
  r = plotPed(x, labs = labs, draw = FALSE, margins = 0)

  txt = r$parts$texts
  last = txt$label == labs[5L]
  h = graphics::strheight(txt$label[last], units = "inches",
                        cex = txt$cex[last]) / r$scale$vscale
  expect_gte(txt$y[last] - h,
             r$scale$usr[3L] - sqrt(.Machine$double.eps))
})


test_that("titles increase top margin only when needed", {
  x = data.frame(id = 1:3, fid = c(0, 0, 1),
                 mid = c(0, 0, 2), sex = c(1L, 2L, 1L))
  a = alignPed(x)

  s0 = scalePed(a, margins = 1)
  s1 = scalePed(a, margins = 1, title = "Title")
  expect_gt(s1$margins[3L], s0$margins[3L])
  expect_equal(s1$marginType, "mar")

  s2 = scalePed(a, margins = c(1, 1, 5, 1), title = "Title")
  expect_equal(s2$margins[3L], 5)

  s3 = scalePed(a, margins = 1, marginsUnit = "mm", title = "Title")
  expect_equal(s3$marginType, "mai")
  expect_gt(s3$margins[3L], 1 / 25.4)
})

