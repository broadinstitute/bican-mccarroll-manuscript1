test_that("default works", {
  expect_equal(hello(), 'Hello, world!')
})
test_that("hello function works with mom argument", {
  expect_equal(hello(mom = TRUE), "Hi, mom!")
  expect_equal(hello(mom = FALSE), "Hello, world!")
})

test_that("hello function works with invalid mom argument", {
  expect_error(hello(mom = "not a logical"), "Argument 'mom' must be TRUE or FALSE.")
})
