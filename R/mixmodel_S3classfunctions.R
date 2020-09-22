#######################################################################
# socmixmodel class functions
#######################################################################

#' Create a soc mix model object
#'
#' Creates an S3 object of class \code{socmixmodel_model}. From the outhput of
#' a binomail association mixture model (\code{fit_binom_mixture_model})
#'
#' @return
#' The \code{socmixmodel_model} class is a wrapper for the listed output of a
#' binomial association mixture model. Individual list elements can still be accessed in
#' the usual way (\code{$} or \code{[[]]}).
#'
#' Functionality with generic methods \code{summary}
#' (\code{summary.socmixmod_model}) and plot (\code{print.socmixmod_model}).
#'
#' Also works with customised generic method \code{by_edge} to output the
#' probabilities of a each dyad beloning to each component
#'
#'
#' @examples
#' dat = socmixmod_simulate_data()
#' mod = fit_binom_mixture_model(Den = dat$d, Num = dat$n, J = dat$K, edge.ids = dat$edge.ids)
#' summary(mod)
#' print(mod)
#' by_edge(mod)
#' @export

as_socmixmod_model =
  function(K.in, K.out, Den, Num, Means, Freqs, S, rho, loglik, BIC, AIC, ICL, nrep, mu){
  output = list(K.out = K.out,
                n = Num ,
                d = Den,
                k.mean = Means,
                k.freq = Freqs,
                S = S,
                rho = rho,
                loglik = loglik,
                BIC = BIC,
                AIC = AIC,
                ICL = ICL,
                reps = nrep,
                mu.df = mu)
  class(output) <- "socmixmod_model"
  return(output)
}

#' Generic summary for soc mix model object
#' @seealso \link[socmixmods]{as_socmixmod_model}
#'
#' @export
#'
summary.socmixmod_model = function(obj){
  cat(obj$K.out, "Component Binomial Mixture Model", "\n")
  cat("K =", obj$K.out, "\n", "\n")

  cat("Shannons Entropy (Social Complexity):", "\n")
  cat("S =", round(obj$S,3), "\n", "\n")

  cat("Component Summary:", "\n")
  print(
    data.frame(component = paste("k", seq(1, obj$K.out, 1), sep = ""),
               k.mean = round(obj$k.mean, 3),
               freq = round(obj$k.freq, 3)
    )
  )

  cat("\n", "\n")


  cat("Model Fit:","\n")
  print(
    data.frame(
      measure = c("loglik:", "AIC:", "BIC:", "ICL:"),
      value = c(
        round(obj$loglik, 2),
        round(obj$AIC, 2),
        round(obj$BIC, 2),
        round(obj$ICL, 2)
        )
      )
  )
  cat("\n", "\n")

  cat("Diagnostics:", "\n")
  cat("rho = ", obj$rho, "\n")
  cat("n.reps = ", obj$reps)
  cat("\n", "\n")




}

#' Plotting method for socmixmod_model objects
#'
#' This produces a density plot of edge SRI by component to how the binomial
#' mixture model maps onto the real data. It is important to note that the
#' Binomial Mixture Models do not use SRI in their calculation and this plot is
#' intended simply as a guide.
#'
#' Plot requires ggplot2.
#'
#' @param obj a socmixmodel_model object
#'
#' @return
#' A ggplot density plot.
#'
#' @export

plot.socmixmod_model = function(obj){
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("ggplot2 needed for this function to work. Please install and try again",
         call. = FALSE)
  }
  dat = by_edge(obj)
  dat$SRI = obj$n/obj$d
  dat$likely.k = as.factor(dat$likely.k)
  gplot =
    ggplot2::ggplot(dat)+
    ggplot2::geom_density(ggplot2::aes(SRI, fill = likely.k, colour = likely.k, alpha = 0.5))+
    ggplot2::theme(legend.position = "none")
  print(gplot)
}



####################################################
# socmixmod_fittedmodels
####################################################
#' Create object of class as_socmixmod_fittedmodels
#'
#' @export
as_socmixmod_fittedmodels =
  function(summary.table, allmodels.list, criteria){
    output = list(
      summary = summary.table,
      all.models = allmodels.list,
      fitting.criteria = criteria

    )
    class(output) <- "socmixmod_fittedmodels"
    return(output)
  }

#' Summary generic method for socmixmod_fittedmodels
#'
#' @export

summary.socmixmod_fittedmodels = function(obj){

  cat("Best fitting model by", obj$fitting.criteria, ":", "\n")
  cat("k = ", which.min(obj$summary[,obj$fitting.criteria]), sep ="")
  cat("\n", "\n")

  cat("Fitted models summary: ", "\n")
  print(
    data.frame(
    round(obj$summary, 3)
    )
)
}




