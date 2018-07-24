

#Our test to test differences in multinomial populations (based on Central Limit Theorems for Multinomial Sums):

#' Test two multinomial datasets
#'
#' multinom.test peforms the statistical hypothesis test for two multinomial vectors:
#' \deqn{H_o: P_1 = P_2} vs \deqn{H_1: P_1 \ne P_2}
#' where \eqn{P_1} and \eqn{P_2} are \eqn{k} dimensional probability vectors corresponding to two multinomial distributions. The test statistic is defined as
#' \deqn{TODO}
#'
#'
#'
#' @param x,y Vectors of counts or matrices containing multiple count vector observations. x and y must be the same type and dimension. If x and y are matrices (or data frames), the ith row of x will be tested against the ith row of y for all i in [1,nrow(x)].
#' @return The \code{statistic} and its associated \code{p-value}.
#' @seealso
#' Amanda Plunkett & Junyong Park (2018) \emph{Two-Sample Test for Sparse High Dimensional Multinomial Distributions}, TEST, \url{https://doi.org/10.1007/s11749-018-0600-8}
#' @examples
#' #Generate data from two different distributions:
#' data <- genMultinomialData(null_hyp=FALSE)
#'
#' #Perform test:
#' result <- multinom.test(x=data[[1]],y=data[[2]])
#' result
#'
#' #Calculate power of test:
#' mean(result$pvalue)
multinom.test <- function(x,y){

  #TODO: Finish: If given two matrices, need to test rows separately.

  #check that x and y are the same structures:
  if(class(x) != class(y)){
    stop("x and y must be the same class")
  }

  #Check that the same number of categories are defined for x and y:
  if( (is.data.frame(x) | is.matrix(x)) & ncol(x) != ncol(y) |
      (is.vector(x) | is.array(x)) & length(x) != length(y)   ){
    stop("x and y should have the same number of categories
         (i.e. same length if x and y are vectors or same number of columns if
         x and y are matrices or dataframes)")
  }

  if((is.data.frame(x) | is.matrix(x)) & nrow(x)>1 ){
    multipleSamples <- TRUE
  }

  data <- list(x,y)

  p_hat <- list() #p_hat is each element divided by the row sum.
  sum <- 0
  n <- list()

  for(g in 1:2){ #Two groups

    if(multipleSamples){
      n[[g]] <- rowSums(data[[g]])
    }else{
      n[[g]] <- sum(data[[g]])
    }

    p_hat[[g]] <- data[[g]]/n[[g]]
  }

  D <- (p_hat[[1]]-p_hat[[2]])^2 - p_hat[[1]]/n[[1]] - p_hat[[2]]/n[[2]]

  if(multipleSamples){
    stat_numerator <- rowSums(D)
  }else{
    stat_numerator <- sum(D)
  }

  D_var <- 0
  if(multipleSamples){
    #TODO
  }else{
    for(g in 1:2){
      D_var <- D_var + (2/n[[g]]^2)* sum(p_hat[[g]]^2 + p_hat[[g]]/n[[g]])
    }
    #integer overflow:
    stat_var <- D_var + ( 4 * sum(p_hat[[1]] * p_hat[[2]]) /n[[1]] ) /n[[2]]
  }

  stat_denominator <- sqrt(stat_var)

  stat <- stat_numerator/stat_denominator

  stat_pvalue <- 1-pnorm(stat)

  return(list(statistic=stat,pvalue=stat_pvalue))
}





#################### Generate multinomial data: ###################

#' Generate multinomial data
#'
#' Generate two sets of multinomially distributed vectors. Useful for hypothesis testing simulations.
#' Six different experiments with different probability vectors are available in addition to user-specified probability vector \code{p}:
#' \itemize{
#' \item Experiment1: TODO
#' \item Experiment2: TODO
#' \item Experiment3: TODO
#' \item Experiment4: TODO
#' \item Experiment5: TODO
#' \item Experiment6: TODO
#' }
#'
#' @param p An optional 2 by k matrix specifying the probabilities of the k categories for each of the two groups.
#' If defined the rest of the function parameters will not be used. Default value is NULL.
#' @param null_hyp logical; if TRUE, generate data using the same distribution.
#' @param k dimension (number of categories). Default 2000.
#' @param n Vector of length 2 specifying the parameter of the multinomial distribution used to specify the
#' total number of objects that are put into k bins in the typical multinomial experiment.
#' @param sample_size A vector of length 2 specifying the number of multinomial vectors to generate
#' @param expID Experiment number 1-6
#' @param alpha Default is 0.45. Used for experiments 5 and 6.
#' @param numzero Default is 50. Used for experiments 1-3.
#' @param beta Default is 0.25. TODO
#' @return A list containing two dataframes each having dimension sample_size by k.
#' @examples
#' #Generate data:
#' X <- genMultinomialData(null_hyp=TRUE)
#' #Look at the first 10 rows and columns of the first matrix:
#' X[[1]][1:10,1:10]
genMultinomialData <- function(p=NULL,null_hyp=TRUE,k=2000,n=c(8000,8000),sample_size=c(30,30),expID=1,alpha=0.45,numzero=50,beta=0.25){

  #Generate p if it is not given
  if(is.null(p)){
    p <- matrix(NA,2,k)
    sum <- 0

    #First group:
    if(expID==1 | expID==2 | expID==3){  #uniform
      p[1,] <- rep(1/k,times=k)
    }else if(expID==4){  #step
      d <- floor(k/2)
      a <- 8  #constant
      b <- 2  #constant
      p[1,] <- c( rep(a/k,times=d) , rep(b/k,times=k-d))
      p[1,] <- p[1,]/sum(p[1,])
    }else if(expID==5 | expID==6){  #inverse or shift
      p[1,] <- (1/(1:k)^alpha)
      p[1,] <- p[1,]/sum(p[1,])
    }

    if(null_hyp){ #null hypothesis is TRUE
      p[2,] <- p[1,]
    }else{ #null hypothesis is FALSE
      if(expID==1){ #uniform
        #Zero some out and make the rest uniform (Exp 4):
        num_nonzero = k-numzero
        p[2,1:num_nonzero] <- rep(1/num_nonzero,times=num_nonzero)
        p[2,(num_nonzero+1):k] <- 0
      }else if(expID==2){
        #Zero some out and add to a single probability (exp 3):
        p[2,] <- p[1,]
        p[2,numzero+1] <- sum(p[2,1:(numzero+1)])
        p[2,1:numzero] <- 0
      }else if(expID==3){
        #Instead of zeroing some out, use a step function.
        d <- numzero
        a <- 2  #constant
        b <- 8  #constant
        p[2,] <- c( rep(a/k,times=d) , rep(b/k,times=k-d))
        p[2,] <- p[2,]/sum(p[2,])

      }else if(expID==4){
        p[2,] <- rev(p[1,]) #Reverse order
      }else if(expID==5){ #inverse
        #Flip order two entries:
        p[2,] <- p[1,]
        w <- 1  #use 3 here for the dimension and ratio tests.
        v <- 50
        p[2,w] <- p[1,v]
        p[2,v] <- p[1,w]

      }else if(expID==6){ #shift
        p[2,] <- 1/(1:k + beta)
        p[2,] <- p[2,]/sum(p[2,])
      }

    }
  }else{ #check to verify that p is 2xk and the rows sum to 1
    if(nrow(p) != 2){
      stop("Error: p must be a 2 by k dimensional matrix.")
    }
    if(any(rowSums(p) != 1)){
      stop("Error: Row sums of p must be equal to 1.")
    }
  }

  #Generate data
  X <- list()
  for(g in 1:2){
    X[[g]] <- t(rmultinom(sample_size[g],n[g],p[g,]))
  }

  return(X)
}
