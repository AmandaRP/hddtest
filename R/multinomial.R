
#' Test two multinomial datasets
#'
#' Peforms the test for two multinomial vectors
#' testing \eqn{H_0:} the underlying multinomial probability vectors are the
#' same vs. \eqn{H_1:} they are different.
#'
#' @param x,y Integer vectors (or matrices/dataframes containing multiple
#' integer vector observations as rows). \code{x} and \code{y} must be the
#' same type and dimension. If \code{x} and \code{y} are matrices (or
#' dataframes), the \eqn{i^th} row of \code{x} will be tested against the
#' \eqn{i^th} row of \code{y} for all \eqn{i} in 1..\code{nrow(x)}.
#' Alternatively, \code{x} can be a list of two vectors, matrices, or
#' dataframes to be compared. In this case, \code{y} is NULL by default.
#' @return The \code{statistic} and its associated \code{p-value}.
#' If \code{x} and \code{y} are either matrices or dataframes, a
#' \code{statistic} and \code{p-value} will be returned for each row.
#'
#' @seealso
#' Amanda Plunkett & Junyong Park (2018) \emph{Two-Sample Test for Sparse High
#' Dimensional Multinomial Distributions}, TEST,
#' \url{https://doi.org/10.1007/s11749-018-0600-8}
#' @examples
#' #Generate two vectors from the same distribution:
#' data <- genMultinomialData(sample_size=1)
#'
#' #Perform test (the following calls are equivalent):
#' multinom.test(x=data[[1]],y=data[[2]])
#' multinom.test(data)
#' library(magrittr)
#' data %>% multinom.test()
#'
#' #Generate 1000 vectors from each of two different distributions:
#' data <- genMultinomialData(null_hyp=FALSE,sample_size=1000)
#'
#' #Perform test (compare the ith row of x to the ith row of y for all rows):
#' result <- multinom.test(x=data[[1]],y=data[[2]])
#'
#' #Return power of test at the alpha=0.05 level:
#' alpha <- 0.05
#' mean(result$pvalue<alpha)

multinom.test <- function(x,y=NULL){

  if(!is.null(y) ){
    data <- list(x,y)
  }else if(class(x)=="list"){
    data <- x
  }else{
    stop("If y is NULL, x must be a list of two vectors, matrices,
         or dataframes.")
  }

  #check that x and y are the same structures:
  if(class( data[[1]] ) != class( data[[2]] )){
    stop("The structures being compared must be the same class")
  }

  #Check that the same number of categories are defined for x and y:
  if( ((is.data.frame( data[[1]] ) || is.matrix( data[[1]] )) && any(dim( data[[1]] ) != dim( data[[2]] ))  ) ||
      ((is.vector( data[[1]] )     || is.array( data[[1]] ))  && length( data[[1]] ) != length( data[[2]] )) ){
    stop("x and y should have the same number of categories
     (i.e. same length if x and y are vectors or same number of columns if
         x and y are matrices or dataframes)")
  }


  if((is.data.frame( data[[1]] ) || is.matrix( data[[1]] )) && nrow( data[[1]] )>1  ){
    multipleSamples <- TRUE
  }else{
    multipleSamples <- FALSE
  }


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


  if(multipleSamples){
    D_var <- rep(0,nrow(data[[1]]))
    for(g in 1:2){
      D_var <- D_var + (2/n[[g]]^2)* rowSums(p_hat[[g]]^2 + p_hat[[g]]/n[[g]])
    }
    stat_var <- D_var + ( 4 * rowSums(p_hat[[1]] * p_hat[[2]]) /n[[1]] ) /n[[2]]
  }else{
    D_var <- 0
    for(g in 1:2){
      D_var <- D_var + (2/n[[g]]^2)* sum(p_hat[[g]]^2 + p_hat[[g]]/n[[g]])
    }
    stat_var <- D_var + ( 4 * sum(p_hat[[1]] * p_hat[[2]]) /n[[1]] ) /n[[2]]
  }

  stat_denominator <- sqrt(stat_var)

  stat <- stat_numerator/stat_denominator

  stat_pvalue <- 1-pnorm(stat)

  return(list(statistic=stat,pvalue=stat_pvalue))
}




#################### Neighborhood Test ############################

#' Perform the neighborhood test for multinom.test
#'
#' Peforms the test for two multinomial vectors
#' testing \eqn{H_0:} the underlying multinomial probability vectors are the
#' same vs. \eqn{H_1:} they are different.
#'
#'
#' In testing the equality of parameters from two populations
#' (as in \code{multinom.test}),
#' it frequenly happens that the null hypothesis is rejected even though the estimates
#' of effect sizes are close to each other. However, these differences are so small
#' that the parameters may not be considered different in practice. A neighborhood test
#' is useful in this situation.
#'
#' @param x,y Integer vectors (or matrices or dataframes containing multiple
#' integer vector observations as rows). \code{x} and \code{y} must be the
#' same type and dimension. If \code{x} and \code{y} are matrices (or
#' dataframes), the \eqn{i^th} row of \code{x} will be tested against the
#' \eqn{i^th} row of \code{y} for all \eqn{i} in 1..\code{nrow(x)}.
#' Alternatively, \code{x} can be a list of two vectors, matrices, or
#' dataframes to be compared. In this case, \code{y} is NULL by default.
#' @param delta A number (or vector) greater than 0. If not defined, then
#' \code{multinom.test} will be performed (which is equivalent to
#' the neighborhood test with \code{delta}=0).
#'
#' @return The \code{statistic} from \code{multinom.test} and its
#' associated \code{p_delta}, where \eqn{p_delta=pnorm(T - delta)}.
#' TODO: Why subtract from 1 in code?
#' If \code{x} and \code{y} are two dimensional (that is, they are matrices
#' or dataframes with more than one row) and/or
#' \code{delta} is a vector, then a matrix will be returned where the
#' \eqn{i,j^{th}} entry will be the \code{p.delta} associated with the
#' \eqn{i^{th}} row of \code{x} (and \code{y}) and the \eqn{j^{th}} entry of
#' the \code{delta} vector.
#'
#'
#' @seealso
#' Amanda Plunkett & Junyong Park (2018) \emph{Two-Sample Test for Sparse High
#' Dimensional Multinomial Distributions}, TEST,
#' \url{https://doi.org/10.1007/s11749-018-0600-8}
#'
#' @examples
#' require(magrittr)
#'
#' #Load the twoNewsGroups dataset
#'
#' data(twoNewsGroups)
#'
#' #Sample two sets of 200 documents from the sci.med newsGroup (to simulate
#' #    the null hypothesis being TRUE). For each of the two groups, sum the
#' #    200 term frequency vectors together. They will be the two vectors that
#' #    we test.
#'
#' num_docs <- 200
#' vecs2Test <- list(NA,2)
#' rowIDS <- 1:nrow(twoNewsGroups$sci.med)
#' group_1 <- sample(rowIDS,num_docs)
#' group_2 <- sample(rowIDS[-group_1],num_docs)
#'
#' vecs2Test[[1]] <- twoNewsGroups$sci.med[group_1,] %>%
#'                                     colSums() %>% matrix(nrow=1)
#' vecs2Test[[2]] <- twoNewsGroups$sci.med[group_2,] %>%
#'                                     colSums() %>% matrix(nrow=1)
#'
#' #Test the null that the two vectors come from the same distribution
#' #    (i.e. the same news group)
#'
#' vecs2Test %>% multinom.test()
#'
#' #The above test likely produced a significant p-value meaning that we would
#' #    reject the null. However, the difference isn't very interesting. Instead,
#' #    test that the differences are within some neighborhood:
#'
#' vecs2Test %>% multinom.neighborhood.test(delta=60)
#'
#' #How to choose the appropriate delta? The answer may come from subject
#' #    matter expertise about the problem domain. Or, run a simulation to
#' #    gain insight:
#'
#' #Simulation function:
#' simulation <- function(data,null_hyp,delta,reps=10,num_docs=c(200,200)){
#'   vecs2Test <- list( matrix(NA,reps,ncol(data[[1]])), matrix(NA,reps,ncol(data[[1]])) )
#'   for(i in 1:reps){
#'     if(null_hyp){
#'       rowIDS <- 1:nrow(data[[2]])
#'       group_1 <- sample(rowIDS,num_docs[2])
#'       group_2 <- sample(rowIDS[-group_1],num_docs[2])
#'       vecs2Test[[1]][i,] <- data[[2]][samples[group1],] %>% colSums()
#'       vecs2Test[[2]][i,] <- data[[2]][samples[-group1],] %>% colSums()
#'     }else{
#'       vecs2Test[[1]][i,] <- data[[1]][sample(1:nrow(data[[1]]),num_docs[1]),] %>%
#'                                      colSums()
#'       vecs2Test[[2]][i,] <- data[[2]][sample(1:nrow(data[[2]]),num_docs[2]),] %>%
#'                                      colSums()
#'     }
#'   }
#'   result <- vecs2Test %>% multinom.neighborhood.test(delta=delta)
#' } #end simulation function
#'
#' #Run simulation for varying values of delta:
#'
#' delta <- 1:160
#' p.delta.null <- simulation(data=twoNewsGroups,null_hyp=TRUE,delta=delta)$pvalue_delta
#' p.delta.alt  <- simulation(data=twoNewsGroups,null_hyp=FALSE,delta=delta)$pvalue_delta
#'
#' #Plot:
#' par(xpd=TRUE, mar=par()$mar+c(0,0,0,5))
#' matplot(delta, cbind(t(p.delta.null), t(p.delta.alt)), type="l",ylab="p.delta",main="P-Value Curves for Simulation",col=c( rep("red",nrow(p.delta.null)),  rep("blue",nrow(p.delta.alt)) ) )

multinom.neighborhood.test <- function(x,y=NULL,delta=NULL){

  if(is.null(delta) | any(delta <=0)){
    stop("Please specify delta>0")
  }
  result <- multinom.test(x,y)

  #subtract every value of delta from every value of the statistic:
  diff <- lapply(X=result$statistic,FUN=function(X,delta){X-delta},delta=delta)

  mat <- matrix(unlist(diff),length(diff),length(delta),byrow=TRUE)
  p_delta <- 1 - pnorm(q=mat)

  return(list(statistic=result$statistic,pvalue_delta=p_delta))

}


#################### Generate multinomial data: ###################

#' Generate multinomial data
#'
#' Generate two sets of multinomially distributed vectors using
#' \code{rmultinom}. Useful for hypothesis testing simulations. Three different
#' experiments with different probability vectors (of length \eqn{k}) are
#' available in addition to user-specified probability vector \code{p}:
#' \itemize{
#' \item Experiment 1: \eqn{p_{1i} = \frac{1/i^\alpha}{\sum_1^k 1/i^\alpha}}.
#'       When the \code{null_hyp} parameter is FALSE, the probability vector
#'       for the 2nd group is generated by switching the position of 1st and
#'       \eqn{m^th} entries.
#' \item Experiment 2: \eqn{p_{1i} = 1/k}. When the \code{null_hyp} parameter is FALSE,
#'       \eqn{p_{2i} = 0} for \eqn{i \in 1...b} and
#'       \eqn{p_{2,b+1}= \sum_{1}^{b+1} p_{1i} = (b+1)/k }.
#' \item Experiment 3: \eqn{p_{1i} = 1/k}. When the \code{null_hyp} parameter is FALSE,
#'       \eqn{p_{2i} = 0} for \eqn{i \in 1...b} and \eqn{p_{2i} = 1/(k − b)} for
#'       \eqn{i > b}.
#' }
#'
#' @param null_hyp logical; if TRUE, generate data using the same distribution.
#' Default value is TRUE.
#' @param p An optional 2 by \eqn{k} matrix specifying the probabilities of the
#' \eqn{k} categories for each of the two groups. Each row of \code{p} must sum
#' to 1. If defined, all remaining parameters in the function definition are
#' ignored. Default value is NULL.
#' @param k integer representing dimension (number of categories). Default 2000.
#' @param n Vector of length 2 specifying the parameter of each multinomial
#' distribution used to define the total number of objects that are put into
#' \eqn{k} bins in the typical multinomial experiment.
#' @param sample_size integer specifying the number of random vectors to generate
#' for each of the two groups.
#' @param expID Experiment number 1-3. Default is 1.
#' @param alpha Number between 0 and 1. Used for experiment 1. Default is 0.45.
#' @param m integer between 2 and \eqn{k}. Used in experiment 1 for the
#' alternative hypothesis. Default is 1000.
#' @param numzero integer between 1 and \eqn{k}-1. Used in experiments 2 and 3
#' for the alternative hypothesis. Default is 50.
#' @param ... Additional parameters.
#' @return A list containing two matrices each having dimension
#' \code{sample_size} by \eqn{k}.
#' @examples
#' #Generate data when the null hypothesis is FALSE:
#' X <- genMultinomialData(FALSE)
#'
#' #Dimension of the two generated datasets:
#' lapply(X,dim)
#'
#' #Proportion of entries less than 5 in the first dataset:
#' sum(X[[1]]<5)/(nrow(X[[1]])*ncol(X[[1]]))
genMultinomialData <- function(null_hyp=TRUE,p=NULL,k=2000,n=c(8000,8000),
                               sample_size=30,expID=1,alpha=0.45,m=1000,
                               numzero=50,...){

  #Generate p if it is not given
  if(is.null(p)){
    p <- matrix(NA,2,k)
    sum <- 0

    #First group:
    if(expID==1){  #inverse or shift

      if(alpha<0 || alpha>1){
        stop("Error: alpha must be between zero and one.")
      }else if(!(m %in% 1:k)){
        stop("Error: m must be an integer between 1 and k")
      }

      p[1,] <- (1/(1:k)^alpha)
      p[1,] <- p[1,]/sum(p[1,])
    }else if(expID==2 | expID==3){  #uniform

      if(!(numzero %in% 1:(k-1))){
        stop("Error: numzero must be an integer in the range 1 to k-1.")
      }

      p[1,] <- rep(1/k,times=k)
    }else{
      stop("Error: expID must be 1, 2, or 3")
    }

    if(null_hyp){ #null hypothesis is TRUE
      p[2,] <- p[1,]

    }else{ #null hypothesis is FALSE
      if(expID==1){

        #Flip order two entries:
        p[2,] <- p[1,]
        if(!exists("w")){
          w <- 1
        }
        p[2,w] <- p[1,m]
        p[2,m] <- p[1,w]

      }else if(expID==2){

        #Zero some out and add to a single probability:
        p[2,] <- p[1,]
        p[2,numzero+1] <- sum(p[2,1:(numzero+1)])
        p[2,1:numzero] <- 0

      }else if(expID==3){

        #Zero some out and make the rest uniform:
        num_nonzero = k-numzero
        p[2,1:num_nonzero] <- rep(1/num_nonzero,times=num_nonzero)
        p[2,(num_nonzero+1):k] <- 0

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
    X[[g]] <- t(rmultinom(sample_size,n[g],p[g,]))
  }

  return(X)
}


#' Document term matrix for documents sampled from two newsgroups
#'
#' A dataset containing two document term matrices for subsets
#' of two newsgroups (rec.sport.baseball and sci.med)
#' from the 20 newsgroups dataset.
#'
#' @format A list of two matrices, each having dimension 594 by 16214.
#' The (i,j) entry of each matrix is the count (term frequency) of
#' the jth word in the ith document. The first matrix in the list
#' contains 594 sampled documents from the rec.sport.baseball
#' newsgroup. The second contains 594 sampled documents from the
#' sci.med newsgroup.
#' @source \url{http://qwone.com/~jason/20Newsgroups/}
"twoNewsGroups"
