#' Simulate Data
#'
#'Function to simulate data with a known number of components (K) to use with
#'socmixmod functions. See Weiss et al 2019 for details.
#'
#'@param K Number of components to simulate in the data (default 3)
#'@param N Number of individuals to simulate (default = 20)
#'@param mean.d Mean sampling effort per dyad (default = 80)
#'
#'
#'@return
#'Output is a list with elements:
#'#' \itemize{
#'  \item \code{d}: SRI Denominator. Number of times simulated dyad are
#'  observed in total.
#'  \item \code{n}: SRI Numerator. Number of times simulated dyad are observed
#'  together.
#'  \item \code{K}: Number of components (will match inputted K)
#'  \item \code{edge.ids}: simulated edge id names
#' }
#'
#' @examples
#' socmixmod_simulate_data()
#'
#' socmixmod_simulate_data(K = 5, N = 100, mean.d=50)
#'@export
socmixmod_simulate_data = function(
  K = 3, # number of types (K/J)
  N = 20, # number of individuals
  mean.d = 80 #sampling effort
){

  N.dyad = (N*(N-1))/2 #number of dyads
  bad.par = T
  while(bad.par){
    mu = runif(K,0,1)
    a = runif(K,0,1)
    a = a/sum(a)
    if(K>1){
    if(min(dist(mu))>=0.1 & min(a)>=0.1/K) bad.par = F
    } else {
      if(min(a)>=0.1/K) bad.par = F
    }
  }
  rho = runif(K,0,0.015) #overdispersion
  b1 = mu*(1/rho - 1) #shape parameters from means and overdispersion
  b2 = ((mu-1)*(rho-1))/rho
  k = sample(K, size = N.dyad, rep = T, prob = a) #assign classes
  p = rbeta(n=N.dyad,shape1=b1[k],shape2=b2[k]) #assign association probabilities
  d = rpois(N.dyad,mean.d) #assign denominators
  x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators
  if(requireNamespace("babynames", quietly = TRUE)){
    bn = babynames::babynames
    ids = sample(unique(bn$name), N)
  } else {
    ids = as.roman(seq(1,N, 1))
  }
  edgeids = do.call(rbind, combn(ids, 2, simplify = FALSE))
  edgeids = paste(edgeids[,1], "-", edgeids[,2], sep = "")
  edgeids = sample(edgeids)

  output = list(
    d = d,
    n = x,
    K = K,
    edge.ids = edgeids

  )
  return(output)
}


#' Extract and tidy the by edge data in a socmixmodel_model class object
#'
#' This function takes the data from a socmixmodel_model class object (produced
#' with \code{fit_binom_mixture_model}) and produces a table of data by edge. It
#' also calculates the component a given edge belongs to based on the one which
#' if has the highest probability of belonging to.
#'
#' @param obj a socmixmodel_model class object
#'
#' @return A dataframe showing each edge in the population, the probability
#' they belong to each component \code{P(ki)} and the component the edge
#' is most likely to belong to
#'
#' @examples
#' #simulate data
#' dat = socmixmod_simulate_data()
#'
#' # Run model
#' fitted.model =
#' fit_binom_mixture_model(Den = dat$d, Num = dat$n, J = dat$K, edge.ids = dat$edge.ids)
#'
#' #View output
#' fitted.model
#' summary(fitted.model)
#' by_edge(fitted.model)
#'
#'@export
#'
by_edge = function(obj){

  mu = obj$mu

  mu.df = data.frame(
    edge.id = rownames(mu),
    round(mu,6)
  )
  rownames(mu.df) = NULL
  names(mu.df)[2:ncol(mu.df)] = paste(rep.int("P(k",obj$K.out), seq.int(1,obj$K.out,1),")",sep = "")

  mu.df$likely.k = numeric(nrow(mu.df))
  for(i in 1:nrow(mu.df)){
    mu.df$likely.k[i] = which(mu.df[i, 2:ncol(mu.df)] ==
                                max(mu.df[i, 2:ncol(mu.df)]))
  }
  return(mu.df)
}



#' Get the best model method from a  socmixmod_fittedmodels object
#'
#'This function extracts the best model, as a \code{socmixmodel_model] object,
#'from the \code{socmixmod_fittedmodels] object produced by
#'\code{binom_assoc_mixt} function/
#'
#'@param obj  socmixmod_fittedmodels class object
#'@param verbose Logical. If TRUE displays summary information from the best fitting
#'model when applied. Default=TRUE.
#'
#'@return A socmixmodel_model of the best fitting model as specified by the
#'\code{criterion} choosen for model fitting.
#'
#' @examples
#'# simulate data
#' dat = socmixmod_simulate_data()
#'
#' # Run model
#' model.fitting =
#' binom_assoc_mixt(Den = dat$d, Num = dat$n, edge.ids = dat$edge.ids, run.all.J = TRUE)
#'
#' #View output
#' model.fitting
#' summary(model.fitting)
#' best.mod = get_bestmodel(model.fitting)
#' summary(best.mod)
#'
#'@seealso
#'binom_assoc_mixt
#'
#'@export
#'
get_bestmodel = function(obj, verbose = TRUE){
  bestmod.i = which.min(obj$summary[,obj$fitting.criteria])
  bestmod = obj$all.models[[bestmod.i]]
  if(verbose == TRUE){
    summary(bestmod)
  }
  return(bestmod)
}
