

#Our test to test differences in multinomial populations (based on Central Limit Theorems for Multinomial Sums):

#' Test two multinomial datasets.
#'
#' @param data A list of two matrices
#' @return The \code{statistic} and its associated \code{p-value}
#' @example
#'
multinom.test <- function(data){

  #TODO: Does this code work for both matrices and vectors?
  #TODO: Check that two datasets have the same number of columns.

  p_hat <- list() #p_hat is each element divided by the row sum.
  sum <- 0
  n <- list()

  for(g in 1:2){ #Two groups
    n[[g]] <- sum(data[[g]])
    p_hat[[g]] <- data[[g]]/n[[g]]
  }

  D <- (p_hat[[1]]-p_hat[[2]])^2 - p_hat[[1]]/n[[1]] - p_hat[[2]]/n[[2]]

  stat_numerator <- sum(D)

  D_var <- 0
  for(c in 1:2){
    D_var <- D_var + (2/n[[c]]^2)* sum(p_hat[[c]]^2 + p_hat[[c]]/n[[c]])
  }
  #integer overflow:
  #D_var <- D_var + ( 4/(n[[1]]*n[[2]]) ) * sum(p_hat[[1]] * p_hat[[2]])
  stat_var <- D_var + ( 4 * sum(p_hat[[1]] * p_hat[[2]]) /n[[1]] ) /n[[2]]

  stat_denominator <- sqrt(stat_var)

  stat <- stat_numerator/stat_denominator

  stat_pvalue <- 1-pnorm(stat)

  return(list(statistic=stat,pvalue=stat_pvalue))
}
