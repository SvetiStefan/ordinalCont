#' Generalized logistic g function
#'
#' This function computes a parametric version of the g function following Richards (1959): 
#' \deqn{g(v) = M + \frac{1}{B} \log\left(\frac{Tv^T}{1-v^T}\right)}
#'  M is omitted as an intercept is always fitted.
#' @param v vector of standardized scores from the continuous ordinal scale, 0<v<1.
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry of the curve. 
#' @keywords Richards, generalized logistic function.
#' @details The generalized logistic functions maps from (0,1) to \eqn{(-\infty,\infty)}. 
#' B is the slope of the curve,  T is the symmetry and M is the offset. 
#' M is absorbed into the intercept of the model and so is set to zero in this function.
#' @return A vector of length equal to the length of v, with values g(v).
#' 
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#


g_glf <- function(v, par){
  return(log(par[2]*v^par[2]/(1-v^par[2]))/par[1])
}


#' Derivative of generalized logistic g function
#'
#' This function computes the derivative of the generalized logistic function as in Richards (1959): 
#' \deqn{\frac{dg(v)}{dv} = \frac{T}{B}  \frac{1}{v(1-v^{T})}}
#' @param v vector of standardized scores from the continuous ordinal scale, 0<v<1.
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry of the curve.
#' @keywords Richards, derivative, generalized logistic function.
#' @return A vector of length equal to the length of v, with values g'(v).
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#

dg_glf <- function(v, par){
  #par = c(B, T)
  return(par[2]/par[1]/v/(1-v^par[2]))
}


#' Inverse of generalized logistic g function
#'
#' This function computes the inverse of a parametric version of the g function 
#' following Richards (1959): 
#' \deqn{g^{-1}(W) = \left( \frac{e^{BW}}{T+e^{BW}}  \right)^{\frac{1}{T}}} M is omitted 
#' as it is absorbed into the 
#' model intercept.
#' @param W vector of scores on the latent scale \eqn{(-\infty,\infty)}. M, the 
#' offset of the g function, is 
#' estimated as the intercept of the model. 
#' \deqn{W=-x'\beta [without intercept] = -x'\beta - M} [as the model is 
#' estimated with an intercept].
#' @param par vector of 2 elements: B, the slope of the curve, and T, the symmetry 
#' of the curve. 
#' @return A vector of length equal to the length of W, with values \eqn{g^{-1}(W)}
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290–301.
#' @keywords Richards, generalized logistic function.

g_glf_inv <- function(W, par){
  exp.part <- exp(par[1]*W)
  return((exp.part/(par[2]+exp.part))^(1/par[2]))
}

