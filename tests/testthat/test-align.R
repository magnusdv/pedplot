test_that("singletons are aligned", {
  x = data.frame(id = 1, fid = 0, mid = 0, sex = 1L)
  a = alignPed(x)
  expect_s3_class(a, "pedalign")
  expect_equal(nrow(a$nodes), 1L)
  expect_equal(a$nodes$x, 0)
})

test_that("nuclear family is aligned", {
  x = data.frame(id = 1:4, fid = c(0, 0, 1, 1),
                 mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
  a = alignPed(x)
  expect_equal(nrow(a$input), 4L)
  expect_true(nrow(a$matings) >= 1L)
  expect_true(all(a$nodes$idx %in% seq_len(4)))
})

test_that("hidden founders are removed from nodes", {
  x = data.frame(id = 1:4, fid = c(0, 0, 1, 1),
                 mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
  a = alignPed(x, hideFounders = 1)
  m = a$matings

  expect_true(a$input$hidden[a$input$id == "1"])
  expect_false("1" %in% a$nodes$id)
  expect_equal(nrow(m), 1L)
  expect_equal(m$faIdx, 1L)
  expect_equal(m$moIdx, 2L)
  expect_true(m$faHidden)
  expect_false(m$moHidden)
  expect_true(m$singleParent)
  expect_equal(m$faNode, 0L)
  expect_true(m$moNode > 0L)
})

test_that("hidden-founder single parents are centred above children", {
  x = data.frame(id = 1:4, fid = c(0, 0, 1, 1),
                 mid = c(0, 0, 2, 2), sex = c(1, 2, 1, 2))
  a = alignPed(x, hideFounders = 1)

  xparent = a$nodes$x[a$nodes$id == 2] # xpos of remanining parent
  xch = a$nodes$x[a$nodes$id %in% 3:4] # xpos of children

  expect_equal(xparent, mean(xch))
})

test_that("hidden-founder shifts propagate to ancestors", {
  x = data.frame(id = 1:5, fid = c(0, 0, 1, 0, 3),
                 mid = c(0, 0, 2, 0, 4), sex = c(1, 2, 1, 2, 1))
  a = alignPed(x, hideFounders = "4")
  m = a$matings[a$matings$faIdx == 1L & a$matings$moIdx == 2L, ]
  parx = mean(a$nodes$x[a$nodes$node %in% c(m$faNode, m$moNode)])
  childx = a$nodes$x[a$nodes$idx == 3L]

  expect_equal(parx, childx)
})


test_that("node copies are explicit", {
  x = data.frame(id = 1:7, fid = c(0, 0, 1, 1, 0, 3, 4),
                 mid = c(0, 0, 2, 2, 0, 5, 5), sex = c(1, 2, 1, 2, 2, 1, 2))
  a = alignPed(x)
  for(i in unique(a$nodes$idx)) {
    jj = a$nodes$idx == i
    expect_equal(a$nodes$copy[jj], seq_len(sum(jj)))
  }
})

test_that("hidden-founder branches move inward under ancestors", {
  skip("not implemented yet")
  x = data.frame(id = 1:8,
                 fid = c(0, 0, 1, 0, 1, 0, 3, 5),
                 mid = c(0, 0, 2, 0, 2, 0, 4, 6),
                 sex = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L))

  a = alignPed(x, hideFounders = "4")
  xx = a$nodes$x
  names(xx) = a$nodes$id

  expect_equal(unname(xx["3"]), unname(xx["1"]))
  expect_equal(unname(xx["7"]), unname(xx["3"]))

  a = alignPed(x, hideFounders = c("4", "6"))
  xx = a$nodes$x
  names(xx) = a$nodes$id

  expect_equal(unname(xx["3"]), unname(xx["1"]))
  expect_equal(unname(xx["5"]), unname(xx["2"]))
  expect_equal(unname(xx["7"]), unname(xx["3"]))
  expect_equal(unname(xx["8"]), unname(xx["5"]))
})

test_that("stacked hidden-founder spouses remain vertical", {
  x = data.frame(id = 1:7,
                 fid = c(0, 0, 1, 0, 3, 0, 5),
                 mid = c(0, 0, 2, 0, 4, 0, 6),
                 sex = c(1L, 2L, 1L, 2L, 1L, 2L, 1L))

  a = alignPed(x, hideFounders = c("2", "4", "6"))
  xx = a$nodes$x
  names(xx) = a$nodes$id

  expect_equal(unname(xx["1"]), unname(xx["3"]))
  expect_equal(unname(xx["3"]), unname(xx["5"]))
  expect_equal(unname(xx["5"]), unname(xx["7"]))
})

test_that("alignment with hidden founders: half cousins", {
  x = data.frame(id = 1:9,
                 fid = c(0, 0, 0, 2, 0, 2, 0, 4, 6),
                 mid = c(0, 0, 0, 1, 0, 3, 0, 5, 7),
                 sex = c(2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 1L))

  a = alignPed(x, hideFounders = c(5,7))
  xx = a$nodes$x
  names(xx) = a$nodes$id

  expect_equal((xx[["1"]]+xx[["2"]])/2, xx[["4"]])
  expect_equal((xx[["1"]]+xx[["2"]])/2, xx[["8"]])
  expect_equal((xx[["3"]]+xx[["2"]])/2, xx[["6"]])
  expect_equal((xx[["3"]]+xx[["2"]])/2, xx[["9"]])
})

