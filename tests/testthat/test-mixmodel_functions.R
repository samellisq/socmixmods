library(socmixmods)
#############################################################
## Generate some data
#############################################################

K = 3 # number of types (K/J)
N = 20 # number of individuals
N.dyad = (N*(N-1))/2 #number of dyads
mean.d = 80 #sampling effort
bad.par = T
while(bad.par){
  mu = runif(K,0,1)
  a = runif(K,0,1)
  a = a/sum(a)
  if(min(dist(mu))>=0.1 & min(a)>=0.1/K) bad.par = F
}
rho = runif(K,0,0.015) #overdispersion
b1 = mu*(1/rho - 1) #shape parameters from means and overdispersion
b2 = ((mu-1)*(rho-1))/rho
k = sample(K, size = N.dyad, rep = T, prob = a) #assign classes
p = rbeta(n=N.dyad,shape1=b1[k],shape2=b2[k]) #assign association probabilities
d = rpois(N.dyad,mean.d) #assign denominators
x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators
names(d) = as.roman(seq(1,length(d),1)) #Add edge names
names(x) = as.roman(seq(1,length(x),1)) # Add edge names

################################################
# Tests
###################################################


context("Check that incorrect inputs are rejected")

test_that("fit_binom_mixture_model incorrect inputs are rejected", {
  expect_error(fit_binom_mixture_model(Den = as.character(d), Num = x, J = 3))
  x.error = x
  x.error[1] = NA
  expect_error(fit_binom_mixture_model(Den = d, Num = x.error, J = 3))
  error.x = x
  error.x = x + d
  expect_error(fit_binom_mixture_model(Den = d, Num = x.error, J = 3))
  expect_error(fit_binom_mixture_model(Den = d, Num = x, J = 3.1))
  expect_error(fit_binom_mixture_model(Den = d, Num = x, J = -3))

})

test_that("binom_assoc_mixt incorrect inputs are rejected", {
  expect_error(binom_assoc_mixt(Den = d,Num = as.character(x),criterion="ICL"))
  x.error = x
  x.error[1] = NA
  expect_error(binom_assoc_mixt(Den = d,Num = x.error,criterion="ICL"))
  error.x = x
  error.x = x + d
  expect_error(binom_assoc_mixt(Den = d,Num = x.error,criterion="ICL"))
  expect_error(binom_assoc_mixt(Den = d,Num = x, minJ = 5, maxJ = 2, criterion="ICL"))
  expect_error(binom_assoc_mixt(Den = d,Num = x,criterion="ILC"))
  expect_error(binom_assoc_mixt(Den = d,Num = x,criterion="ICL", run.all.J = TURE))
})


context("context 2")
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
