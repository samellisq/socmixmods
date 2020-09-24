# # # # Run model
# dat = socmixmod_simulate_data(K=3)
# fitted.model =
#   fit_binom_mixture_model(Den = dat$d, Num = dat$n, K = 1, edge.ids = dat$edge.ids)
# #
# #
# X = fitted.model
#
# summary(X)
#
# by_edge(X)
#
# plot(X)
#
# # #######
# #
# fitted.models = binom_assoc_mixt(Den = dat$d,
#                                  Num = dat$n,
#                                  edge.ids = dat$edge.ids,
#                                  minK = 1,
#                                  maxK = 5)
# Y= fitted.models
# summary(Y)
# A = get_bestmodel(Y)
# summary(A)
# by_edge(A)
#
# B = Y$all.models[[2]]
# class(B)
# summary(B)
