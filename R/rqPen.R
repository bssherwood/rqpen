#' @import quantreg
#' @import hqreg
#' @import plyr
#' @import hrqglas
#' @import data.table
#' @import lifecycle
#' @importFrom splines bs
#' @importFrom graphics lines plot par segments points legend
#' @importFrom stats coef coefficients predict quantile residuals sd xtabs fitted weighted.mean IQR
#' @importFrom Rdpack reprompt
#' @importFrom methods is
#' @useDynLib rqPen, .registration=TRUE
NULL 

#' rqPen: A package for estimating quantile regression models using penalized objective functions.
#' 
#' The package estimates a quantile regression model using LASSO, Adaptive LASSO, SCAD, MCP, elastic net, 
#' and their group counterparts, with the exception of elastic net for which there is no group penalty implementation.
#' 
#' @section rqPen functions:
#' The most important functions are rq.pen(), rq.group.pen(), rq.pen.cv() and rq.group.pen.cv(). These functions 
#' fit quantile regression models with individual or group penalties. The cv functions automate the cross-validation process for selection of tuning parameters. 
#' 
#' @docType package
#' @name rqPen
NULL