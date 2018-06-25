

#Our test to test differences in multinomial populations (based on Central Limit Theorems for Multinomial Sums):

#' Test two multinomial datasets
#'
#' TODO: Add a longer description here.
#'
#' @param data A list of two matrices
#' @return The \code{statistic} and its associated \code{p-value}
#' @example
#' #TODO: Find count dataset included with R.
multinom.test <- function(x,y){

  #TODO: Does this code work for both matrices and vectors?

  #check that x and y are the same structures:
  if(class(x) != class(y)){
    stop("x and y must be the same class")
  }

  #Check that the same number of categories are defined for x and y:
  if( (is.data.frame(x) | is.matrix(x)) & ncol(x) != ncol(y) |
      (is.vector(x) | is.array(x)) & length(x) != length(y)   ){
    stop("x and y should have the same number of categories
         (i.e. same length if x and y are vectors or same number of columsn if
         x and y are matrices or dataframes)")
  }

  data <- list(x,y)

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





#################### Generate multinomial data: ###################
#c=dimension
#n is parameter for multinomial
#p is a g by c matrix where g is the number of groups (assume 2 groups for now)

#' Generate multinomial data
#'
#' TODO: Add a longer description here.
#'
#' @param null_hyp TRUE if null hypthesis is true that the data come from the same distribution
#' @param c dimension (number of categories). Default 2000.
#' @param n Vector of length 2 specifying the parameter of the multinomial distribution used to specify the sample size
#' @param sample_size A vector of length 2 specifying the number of multinomial vectors to generate
#' @param p A vector of length c specifying the probabilities of the multinomial distribution. If defined, will be used instead of c. Default value is NULL.
#' @return A list containing two dataframes eacg having dimension sample_size by c.
#' @example
genMultinomialData <- function(null_hyp=TRUE,c=2000,n=c(8000,8000),sample_size=c(30,30),p=NULL,exp1Idx=NULL,exp2Beta=NULL,exp34numZero=NULL){

  #Generate p if it is not given
  if(is.null(p)){
    p <- matrix(NA,2,c)
    sum <- 0

    #pDist <- "uniform"
    #pDist <- "step"
    pDist <- "inverse"
    #pDist <- "shift"
    alpha <- 0.45 #Let alpha=1 for original experiment (inverse and shift experiments)
    #alpha <- 1

    #First group:
    if(pDist=="uniform"){
      p[1,] <- rep(1/c,times=c)
    }else if(pDist=="step"){
      d <- floor(c/2)
      a <- 8  #constant
      b <- 2  #constant
      p[1,] <- c( rep(a/c,times=d) , rep(b/c,times=c-d))
      p[1,] <- p[1,]/sum(p[1,])
    }else if(pDist=="inverse" | pDist=="shift"){
      p[1,] <- (1/(1:c)^alpha)
      p[1,] <- p[1,]/sum(p[1,])
    }

    if(null_hyp){ #null hypothesis is TRUE
      p[2,] <- p[1,]
    }else{ #null hypothesis is FALSE
      if(pDist=="uniform"){
        if(is.null(exp34numZero)){
          numzero = 50 #use 50 here for dimension and ratio tests.
        }else{
          numzero = exp34numZero
        }

        #Zero some out and make the rest uniform (Exp 4):
        num_nonzero = c-numzero
        p[2,1:num_nonzero] <- rep(1/num_nonzero,times=num_nonzero)
        p[2,(num_nonzero+1):c] <- 0

        #Zero some out and add to one probability (exp 3):
        #p[2,] <- p[1,]
        #p[2,numzero+1] <- sum(p[2,1:(numzero+1)])
        #p[2,1:numzero] <- 0

        #Instead of zeroing some out, use a step function.
        #d <- numzero
        #a <- 2  #constant
        #b <- 8  #constant
        #p[2,] <- c( rep(a/c,times=d) , rep(b/c,times=c-d))
        #p[2,] <- p[2,]/sum(p[2,])

      }else if(pDist=="step"){
        p[2,] <- rev(p[1,]) #Reverse order
      }else if(pDist=="inverse"){
        #p[2,] <- rev(p[1,]) #Reverse order

        #Flip order two entries:
        p[2,] <- p[1,]
        if (is.null(exp1Idx)){
          index <- 3  #use 3 here for the dimension and ratio tests.
        }else{
          index <- exp1Idx
        }
        w <- 50
        #w <- 100
        p[2,index] <- p[1,w]
        p[2,w] <- p[1,index]
        #p[2,1] <- p[1,2]
        #p[2,2] <- p[1,1]
      }else if(pDist=="shift"){
        if(is.null(exp2Beta)){
          beta <- 0.25
        }else{
          beta <- exp2Beta
        }
        p[2,] <- 1/(1:c + beta)
        p[2,] <- p[2,]/sum(p[2,])
      }

    }
  }

  #Generate data
  X <- list()
  for(g in 1:2){
    X[[g]] <- t(rmultinom(sample_size[g],n[g],p[g,]))
  }
  #TODO: What should sample size be?

  return(X)
}
