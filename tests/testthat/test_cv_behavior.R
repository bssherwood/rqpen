
set.seed(42)

x <- matrix(rnorm(800),nrow=100)

y <- 1 + x[,1] - 3*x[,5] + rnorm(100)

suppressWarnings({
  cv_model <- cv.rq.pen(x,y)
})

last_model <- cv_model$models[[length(cv_model$models)]]
penultimate_model <- cv_model$models[[length(cv_model$models) - 1]]

testthat::test_that("We match prior cv.rq.pen behavior.", {
  expect_equal(length(cv_model$models), expected = 89)
  expect_equal(sum(abs(last_model$coefficients) > 0), 1)
  expect_gt(sum(abs(penultimate_model$coefficients) > 0), 1)
})
