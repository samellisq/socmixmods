
#' Fit a \code{J} component binomial mixture model
#'
#' \code{fit_binom_mixture_model} fits a \code{J} component mixture model to social association
#' data. This function will usually be only called through \code{\link{binom_assoc_mixt}} but it
#' can be used as a stand-alone function.
#'
#' @param Den A numeric vector of the association index denominator.
#' @param Num A numeric vector of the association index numerator.
#' @param J Integer number of components to fit.
#' @param maxiter Number of iterations. Default input = 1000.
#' @param maxrep Something. Default input = 50
#' @param tolfun Somehing else. Default input = 1e-06.
#' @param minprior Something else again. Default input = 0.
#'
#' @return Output is a list with the elements (accessed via index \code{[i]} or
#' \code{$element.name}): \cr \cr
#' \code{K.in} Number of components aimed to be fitted (will match input \code{J}).
#' Integer. \cr \cr
#' \code{K.out} Number of components fitted by the model. Returned as a sanity check.
#' Should match \code{K.in} & \code{J}.
#' \code{Mean} Mean association rates of the fitted componenets Numeric vector of
#' length \code{J}. \cr \cr
#' \code{Frequency} ? \cr \cr
#' \code{S} Shannons Entropy. Numeric. \cr \cr
#' \code{rho} ? \cr \cr
#' \code{logLik}Log-likelhood of the fitted model. \cr \cr
#' \code{BIC} Model Baysien Information Criteria \cr \cr
#' \code{AIC} Model Akalie Information Criteria. \cr \cr
#' \code{ICL} Model Inform... \cr \cr
#' \code{nrep} ? \cr \cr
#' \code{mu} Probability that each association in the data (rows) belongs to each
#' fitted component (1-\code{K}: columns). If provided row names will be carried through
#' from \code{names(Y)} and/or \code{names(x)}
#'
#' @examples
#' # Generate simulated social data with a known number of components and sampling effort
#' K = 3 # number of types (K/J)
#' N = 20 # number of individuals
#' N.dyad = (N*(N-1))/2 #number of dyads
#' mean.d = 80 #sampling effort
#' bad.par = T
#' while(bad.par){
#'   mu = runif(K,0,1)
#'   a = runif(K,0,1)
#'   a = a/sum(a)
#'   if(min(dist(mu))>=0.1 & min(a)>=0.1/K) bad.par = F
#' }
#' rho = runif(K,0,0.015) #overdispersion
#' b1 = mu*(1/rho - 1) #shape parameters from means and overdispersion
#' b2 = ((mu-1)*(rho-1))/rho
#' k = sample(K, size = N.dyad, rep = T, prob = a) #assign classes
#' p = rbeta(n=N.dyad,shape1=b1[k],shape2=b2[k]) #assign association probabilities
#' d = rpois(N.dyad,mean.d) #assign denominators
#' x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators
#' names(d) = as.roman(seq(1,length(d),1)) #Add edge names
#' names(x) = as.roman(seq(1,length(x),1)) # Add edge names
#'
#' # Run model
#' fitted.model = fit_binom_mixture_model(Den = d, Num = x, J = 5)
#'
#'
#'
#' @seealso
#' \code{\link{binom_assoc_mixt}}

fit_binom_mixture_model = function(Den,Num,J,maxiter=1000,maxrep=50,tolfun=1e-6,minprior=0){

  n = Den
  Y = Num

  if(is.numeric(n) == FALSE |is.numeric(Y) == FALSE){
    stop("Den and Num must be numeric vectors")
  }
  if(anyNA(n)|anyNA(Y)){
    stop("Den and Num cannot contain NAs")
  }
  if(sum(Y > n)!=0){
    stop("Numerator (Num) of association index must be less than or equal to the denominator (Den)")
  }
  if(J%%1!=0 | J < 1){
    stop("J must be a positive integer")
  }

  fit.rho = function(Q,A,n,Y){
    K = length(Q)
    nY = length(Y)
    nm = matrix(n,ncol=K,nrow=nY)
    Ym = matrix(Y,ncol=K,nrow=nY)
    Am = matrix(A,nrow=nY,ncol=K,byrow=T)
    Qm = matrix(Q,nrow=nY,ncol=K,byrow=T)
    ll.rho = function(r){
      b1 = Qm*(1/r - 1)
      b2 = ((Qm-1)*(r-1))/r
      ll = VGAM::dbetabinom.ab(Ym,nm,b1,b2)*Am
      -sum(log(rowSums(ll)))
    }
    optimise(ll.rho,interval=c(0,1))$minimum
  } #function to fit overdispersion parameter
  mean = list() #lists to hold parameters
  freq = list()
  rho = list()
  mus = list()
  lllq = NULL
  nY = length(Y) #number of dyads
  for(mr in 1:maxrep){ # SE:this is the loop that keeps fitting the mixture model with new parameters
    ll0 = 0
    nm = matrix(n,nrow=nY,ncol=J) #a matrix for the n's
    Ym = matrix(Y,nrow=nY,ncol=J) #a matrix of Y's
    Q = 0.8*runif(J,0,1)  ###WHY 0.8
    A = rep(1/J,J)
    K = J
    for(j in 1:maxiter){ # this seems to be an inner thing which also fits new paprameters
      Q = Q[A>minprior] #drop components with no weight (typically not an issue with fuzzy clustering)
      A = A[A>minprior]
      K = length(A) #get new number of components, SE: thsi will be J if minprior = 0
      nm = matrix(n,ncol=K,nrow=nY)
      Ym = matrix(Y,ncol=K,nrow=nY)# # SE matrix of Y's. short and fat (K = J)
      Am = matrix(A,nrow=nY,ncol=K,byrow=T) #turn component paramters into matrices. SE: by row invesre the process so it fills by row rather than by columsn, makes a long thin matric rather than a short fat one here
      Qm = matrix(Q,nrow=nY,ncol=K,byrow=T) # SE: Q is the random jitter (I think, to be applied)
      ll = dbinom(Ym,nm,Qm)*Am #get likelihoods. Each row is a denometer (observed assocation)
      mu = ll/rowSums(ll) #responsibilities.
      A = colMeans(mu) #get fractions
      Q = colSums(mu*Y)/colSums(mu*n) #new parameters
      if(any(is.na(Q))|any(is.na(A))) break
      lll = sum(log(rowSums(ll))) #log-likelihood
      if(abs(lll-ll0)<tolfun & j>(maxiter/10)) break #check for convergence
      ll0 = lll
    }
    rho[[mr]] = fit.rho(Q,A,n,Y)
    mean[[mr]] = Q #save parameters
    freq[[mr]] = A
    lllq = c(lllq,lll)
    lllm = max(lllq,na.rm=T)
    if(mr>4 & sum(abs(lllm-lllq)<tolfun,na.rm = T)>4) break
  }
  mean = mean[[which.max(lllq)]]
  freq = freq[[which.max(lllq)]]
  rho = rho[[which.max(lllq)]]
  freq = freq[order(mean)]
  mean = mean[order(mean)]
  K = length(mean)
  Qm = matrix(mean,nrow=nY,ncol=K,byrow=T)
  Am = matrix(freq,nrow=nY,ncol=K,byrow=T)
  Ym = matrix(Y, nrow=nY,ncol=K)
  rownames(Ym) = names(Y)
  nm = matrix(n, nrow = nY, ncol= K)
  ll = dbinom(Ym,nm,matrix(mean,nrow=nY,ncol=K,byrow=T))*matrix(freq,nrow=nY,ncol=K,byrow=T)
  mu = ll/rowSums(ll) ### mu here is the probability is the list of probabilities that the particualr association belongs to a given segment. can be used to define the bonds. Divide by row sums to make it 'given the probability of being anywhere what is the probaiblity of it being here'
  qq = mean*freq
  qq = qq/sum(qq)
  qq = qq[qq>0] #prevents NAs in the entropy estimate
  S = -sum(qq*log(qq))
  lllm = sum(log(rowSums(ll)))
  cl = mu[mu!=0]
  vv = -sum(cl*log(cl))
  AIC = 2*(2*J-1) - 2*lllm
  BIC = (2*J-1)+log(nY)*(2*J-1) - 2*lllm
  ICL = BIC + 2*vv
  return(list(K.in = J, K.out = length(mean),Mean=mean,Frequency=freq, S = S, rho = rho, logLik=lllm, BIC = BIC, AIC = AIC, ICL = ICL, nrep=mr, mu = mu))
}

#' Fit Binomial mixture models to  social association data
#'
#' This function fits up to \code{maxJ} binomial mixture models to social data. Returning a summary
#' of the fit of each model and the output of the best fitting model. \cr \cr
#' Where J > 1 this calls function \code{fit_binom_mixture_model}. Where J = 1, fits a standard
#' (non-mixture) binomial model. If minJ=maxJ a single model is run.
#'
#' @param Den A numeric vector of the association index denominator.
#' @param Num A numeric vector of the association index numerator.
#' @param minJ Integer. Minimum component model to fit to the data. Default = \code{1}.
#' @param maxJ Integer. Maximum componenet model to fit to the data. Default = \code{9}.
#' @param nrep Something
#' @param criterion String. Model fitting criteria to use to choose best model.
#' Categorical, options are \code{"ICL"}[default], \code{"AIC"} & \code{"BIC"}.
#' @param run.all.J Logical. If \code{TRUE} fits all component models between \code{minJ}
#' and \code{maxJ}. If \code{FALSE} [default], when two component-models in a row fit less well (by
#' \code{criterion}) stops the analysis and returns output.
#'
#' @return
#' Output is a list with elements (access via index or \code{$element.name}): \cr \cr
#' \code{summary}. A summary table (data frame) summarising the fitted models. Each row is a model
#' fitted with J components with columns:
#' \itemize{
#'  \item \code{K.in}: Input number of components \code{J}
#'  \item \code{K.out}: Number of components fitted by the model (should match \code{K.in}).
#'  \item \code{S}: Shannon's Entropy.
#'  \item \code{rho}: Something
#'  \item \code{AIC}: Model AIC
#'  \item \code{BIC}: Model BIC
#'  \item \code{ICL}: Model ICL
#' }
#' \code{all.models} Output from all models (see \code{\link{fit_binom_mixture_model}})
#' in order from \code{minJ} to \code{maxJ} as list.
#' \cr \cr
#' \code{best.model} The output from the best fitting model, selected via \code{criterion}. Output
#' is list output from \code{\link{fit_binom_mixture_model}}.
#'
#' @examples
#' # Generate simulated social data with a known number of components and sampling effort
#' K = 3 # number of types (K/J)
#' N = 20 # number of individuals
#' N.dyad = (N*(N-1))/2 #number of dyads
#' mean.d = 80 #sampling effort
#' bad.par = T
#' while(bad.par){
#'   mu = runif(K,0,1)
#'   a = runif(K,0,1)
#'   a = a/sum(a)
#'   if(min(dist(mu))>=0.1 & min(a)>=0.1/K) bad.par = F
#' }
#' rho = runif(K,0,0.015) #overdispersion
#' b1 = mu*(1/rho - 1) #shape parameters from means and overdispersion
#' b2 = ((mu-1)*(rho-1))/rho
#' k = sample(K, size = N.dyad, rep = T, prob = a) #assign classes
#' p = rbeta(n=N.dyad,shape1=b1[k],shape2=b2[k]) #assign association probabilities
#' d = rpois(N.dyad,mean.d) #assign denominators
#' x = rbinom(n=N.dyad,size=d,prob=p) #assign numerators
#' names(d) = as.roman(seq(1,length(d),1)) #Add edge names
#' names(x) = as.roman(seq(1,length(x),1)) # Add edge names
#'
#' #Run the models, stopping when best fitting model is reached.
#' model.fit = binom_assoc_mixt(Den = d,Num = x,criterion="ICL")
#' model.fit$summary
#'
#' #Run all models, between minJ and maxJ
#' model.fit = binom_assoc_mixt(Den = d, Num = x, minJ = 1, maxJ = 9, run.all.J = TRUE, criterion="ICL")
#' model.fit$summary


binom_assoc_mixt = function(Den,Num,minJ = 1, maxJ=9,nrep=20,criterion="ICL", run.all.J = FALSE){

  if(is.numeric(Den) == FALSE |is.numeric(Num) == FALSE){
    stop("Den and Num must be numeric vectors")
  }
  if(anyNA(Den)|anyNA(Num)){
    stop("Den and Num cannot contain NAs")
  }
  if(sum(Num > Den)!=0){
    stop("Numerator (Num) of association index must be less than or equal to the denominator (Den)")
  }
  if(minJ%%1!=0 | minJ < 1 | maxJ%%1!=0 | maxJ < 1){
    stop("minJ and maxJ must be a positive integers")
  }
  if(minJ > maxJ){
    stop("maxJ must be larger than minJ")
  }
  if(!(criterion %in% c("AIC","BIC","ICL"))){
    stop("criterion must be AIC, BIC or ICL")
  }
  if(!(is.logical(run.all.J))){
    stop("run.all.J must be either TRUE or FALSE")
  }


  allJs = seq(minJ, maxJ,1)
  summary = matrix(ncol = 7,nrow=maxJ) #matrix for summary output
  models = list() #list to hold all models
  colnames(summary) = c("K.in", "K.out", "S", "rho", "AIC","BIC","ICL")
  worse = 0
  for(i in minJ:maxJ){
    print(paste("Fitting",i,"Component(s)"))
    if(i > 1){
      res = fit_binom_mixture_model(Den,Num,i,maxrep=nrep)
    }else{
      res <- list()
      res$K.out <- 1
      res$S <- 0
      mean <- sum(Num)/sum(Den)
      res$logLik <- sum(dbinom(Num,Den,prob=mean,log=T))
      res$AIC <- 2 - 2*res$logLik
      res$BIC <- 2+log(length(Num))*2 - 2*res$logLik
      cl <- rep(1,length(Num))
      vv = -sum(cl*log(cl))
      res$rho <- optimise(function(z)-sum(VGAM::dbetabinom(Num,Den,mean,rho=z,log=T)),upper=1,lower=0)$minimum
      res$ICL <- res$BIC + 2*vv
      res$mu <- rep.int(1, length(mean))
    }
    summary[i,"K.in"] = i
    summary[i,"K.out"] = res$K.out
    summary[i,"S"] = res$S
    summary[i, "rho"] = res$rho
    summary[i,"AIC"] = res$AIC
    summary[i,"BIC"] = res$BIC
    summary[i,"ICL"] = res$ICL
    if(i > 1){ # should be i > 1
      if(criterion == "AIC"){ #If the current model is worse than the previous, record this.
        worse = ifelse(summary[i,"AIC"]>summary[(i-1),"AIC"], worse+1, 0)
      }
      if(criterion == "BIC"){
        worse = ifelse(summary[i,"BIC"]>summary[(i-1),"BIC"], worse+1, 0)
      }
      if(criterion == "ICL"){
        worse = ifelse(summary[i,"ICL"]>summary[(i-1),"ICL"], worse+1, 0)
      }
    }
    models[[i]] = res
    if(run.all.J == FALSE & !(minJ == maxJ)){
      if(worse > 1) break #If two models in a row have gotten worse, stop fitting
    }
  }
  if(criterion == "BIC"){ #Get the best model, as chosen by the specified criteria
    best = models[[which.min(summary[,"BIC"])]]
    messege = paste("Best Model, k = ", which.min(summary[,"BIC"]), sep ="")
  }
  if(criterion == "AIC"){
    best = models[[which.min(summary[,"AIC"])]]
    messege = paste("Best Model, k = ", which.min(summary[,"AIC"]), sep ="")
  }
  if(criterion == "ICL"){
    best = models[[which.min(summary[,"ICL"])]]
    messege = paste("Best Model, k = ", which.min(summary[,"ICL"]), sep ="")
  }
  summary = summary[min(allJs):i,]

  print(messege)
  return(list(summary = data.frame(summary), all.models = models, best.model = best)) #return list containing summary table and best model
}

