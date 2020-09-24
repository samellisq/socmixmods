library(socmixmods)
#############################################################
## Generate some data
#############################################################

sim.data = socmixmod_simulate_data()

################################################
# Tests: fit_binnm_mixture_model
###################################################


context("Checking fit_binom_mixture_model")

test_that("Check incorrect inputs are rejected", {
  expect_error(fit_binom_mixture_model(Den = as.character(sim.data$d), Num = x, K = 3))
  x.error = sim.data$x
  x.error[1] = NA
  expect_error(fit_binom_mixture_model(Den = d, Num = x.error, K = 3))
  error.x = sim.data$x
  error.x = sim.data$x + sim.data$d
  expect_error(fit_binom_mixture_model(Den = d, Num = x.error, K = 3))
  expect_error(fit_binom_mixture_model(Den = d, Num = x, K = 3.1))
  expect_error(fit_binom_mixture_model(Den = d, Num = x, K = -3))

})


test_that("Check output is of the correct type",{
  mod = fit_binom_mixture_model(Den = sim.data$d,
                                Num = sim.data$n,
                                K = 3,
                                edge.ids = sim.data$edge.ids)
  expect_match(class(mod), "socmixmod_model")

})

context("Checking binom_assoc_mixt")

test_that("binom_assoc_mixt incorrect inputs are rejected", {
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = as.character(sim.data$x),criterion="ICL"))
  x.error = sim.data$x
  x.error[1] = NA
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = sim.data$x.error,criterion="ICL"))
  error.x = sim.data$x
  error.x = sim.data$x + sim.data$d
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = x.error,criterion="ICL"))
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = sim.data$x, minK = 5, maxK = 2, criterion="ICL"))
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = sim.data$x,criterion="ILC"))
  expect_error(binom_assoc_mixt(Den = sim.data$d,Num = sim.data$x,criterion="ICL", run.all.K = TURE))
})

test_that("Check output is of the correct type",{
  mod2 = binom_assoc_mixt(Den = sim.data$d,
                          Num = sim.data$n,
                          edge.ids = sim.data$edge.ids,
                          minK = 1,
                          maxK = 3,
                          criterion = "ICL",
                          run.all.K = TRUE,
                          verbose = FALSE
                          )
  expect_match(class(mod2), "socmixmod_fittedmodels")
  mod2.2 = mod2$all.models[[2]]
  expect_match(class(mod2.2), "socmixmod_model")

})



