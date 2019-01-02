#' Test two multivariate binary datasets
#'
#' Peforms the test for two binary vectors
#' testing \eqn{H_0:} the underlying probability vectors are the
#' same vs. \eqn{H_1:} they are different.
#'
#'
#' Statistic:
#' T = \sum_{j=1}^d D_j^2 I( |Dj| \ge \delta(d))
#' where d is the dimension.
#' Dj = (\hat{p}_{1j} − \hat{p}_{2j} )/sqrt{ \hat{p}_j (1 − \hat{p}_j )(1/n1 + 1/n2) }
#' where \hat{p}_{cj} is the estimate of p_{cj} for the c^{th} group calculated
#' by  x_{cj} . Let \hat{p}_j be the pooled estimate for the j^{th} variable.
#' Additionally, \delta(d) = \sqrt{2 log (a_d d)} where a_d = (log d)^{-2}.
#'
#' P-value: Calculated using the permutation method.
#'
#' Note: As described in the reference below, this method does not perform
#' well on highly correlated variables. See the reference for more details.
#'
#' @param x,y Matrices (or dataframes) containing multiple
#' integer vector observations as rows. \code{x} and \code{y} must be the
#' same type and dimension. Alternatively, \code{x} can be a list of two
#' matrices (or dataframes) to be compared. In this case, \code{y} is NULL
#' by default.
#' @param numPerms Number of permutations to use to calculate the p-value.
#' Default value is 5000.
#' @return The \code{statistic} and its associated \code{p-value}.
#'
#' @seealso
#' Amanda Plunkett & Junyong Park (2017) \emph{Two-sample Tests for Sparse
#' High-Dimensional Binary Data}, Communications in Statistics - Theory and
#' Methods, 46:22, 11181-11193
#'
#' @examples
#' #Binarize the twoNewsGroups dataset:
#' binData <- list(twoNewsGroups[[1]]>0,twoNewsGroups[[2]]>0)
#' names(binData) <- names(twoNewsGroups)
#'
#' #Perform the test:
#' binary.test(binData)
#'
#' #The following are equivalent to the previous call:
#' binData %>% binary.test
#' binary.test(binData[[1]],binData[[2]])
binary.test <- function(x,y=NULL,numPerms=5000){

  if(!is.null(y) ){
    data <- list(x,y)
  }else if(class(x)=="list"){
    data <- x
  }else{
    stop("If y is NULL, x must be a list of two matrices
         or dataframes.")
  }

  #check that x and y are the same structures:
  if(class( data[[1]] ) != class( data[[2]] )){
    stop("The structures being compared must be the same class.")
  }

  #Check that the same number of categories are defined for x and y:
  #TODO: Add tibbles and sparse Matrices
  if( ((is.data.frame( data[[1]] ) || is.matrix( data[[1]] )) && any(dim( data[[1]] ) != dim( data[[2]] ))  ) ){
    stop("x and y must be matrices or dataframes and should have the same number of categories
     (i.e. same number of columns)" )
  }


}
