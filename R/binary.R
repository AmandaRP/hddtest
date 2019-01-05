#' Test two multivariate binary datasets
#'
#' Peforms the test for two binary vectors
#' testing \eqn{H_0:} the underlying probability vectors are the
#' same vs. \eqn{H_1:} they are different.
#'
#'
#' The statistic is \eqn{T = \sum_{j=1}^d D_j^2 I( |Dj| \ge \delta(d))}
#' where \eqn{d} is the dimension of the data. Additionally:
#' \itemize{
#'   \item \eqn{Dj = (\hat{p}_{1j} − \hat{p}_{2j} )/\sqrt{ \hat{p}_j (1 − \hat{p}_j )(1/n1 + 1/n2) } }
#'   \item \eqn{\hat{p}_{cj}} is the estimate of \eqn{p_{cj}} for the \eqn{c^{th}} group calculated by the \eqn{j^th} column mean
#'   \item \eqn{\hat{p}_j} is the pooled estimate for the \eqn{j^{th}} variable.
#'   \item \eqn{\delta(d) = \sqrt{2 log (a_d d)}} where \eqn{a_d = (log d)^{-2}}
#' }
#'
#' The p-value associated with the statistic is calculated using the
#' permutation method. The observation vectors are repeatedly shuffled
#' between groups, each time being used to re-calculate the statistic.
#' A null distribution is constructed and used to calcualate the p-value.
#'
#' @section Warning:
#' As described in the reference below, this method may not perform
#' well on highly correlated variables.
#'
#' @param x,y Matrices (or dataframes) containing multiple
#' integer vector observations as rows. \code{x} and \code{y} must be the
#' same type and dimension. Alternatively, \code{x} can be a list of two
#' matrices (or dataframes) to be compared. In this case, \code{y} is NULL
#' by default.
#' @param numPerms Number of permutations to use to calculate the p-value.
#' Default value is 5000.
#' @return A list containing the computed \code{statistic}, a list of statistics
#' (\code{null.statistics}) used to construct the null distritubution (from the
#' permutation method), and the associated \code{p-value}. The \code{p-value} is
#' the percent of \code{null.statistics} that are more extreme than the
#' \code{statistic} computed from the original dataset.
#'
#' @seealso
#' Amanda Plunkett & Junyong Park (2017), \emph{Two-sample Tests for Sparse
#' High-Dimensional Binary Data}, Communications in Statistics - Theory and
#' Methods, 46:22, 11181-11193
#'
#' @examples
#' #Binarize the twoNewsGroups dataset:
#' binData <- list(twoNewsGroups[[1]]>0,twoNewsGroups[[2]]>0)
#' names(binData) <- names(twoNewsGroups)
#'
#' #Perform the test:
#' mvbinary.test(binData)
#'
#' #The following are equivalent to the previous call:
#' mvbinary.test(binData[[1]],binData[[2]])
#' require(magrittr)
#' binData %>% mvbinary.test
mvbinary.test <- function(x,y=NULL,numPerms=5000){

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

  #Get statistic:
  stat <- get_stat(data)$statistic

  #sample size array:
  n <- as.numeric(lapply(data,FUN=nrow))
  #dimension:
  d <- ncol(data[[1]])

  #Permutation method to compute p-value:
  data <- rbind(data[[1]],data[[2]]) #combine into one dataframe
  data.perm <- list()
  perm.stats <- array(NA,dim=NumPerms)
  for(i in 1:numPerms){
    sortOrder <- sample(1:sum(n))
    data.perm[[1]] <- data[sortOrder[1:n[1]],]
    data.perm[[2]] <- data[sortOrder[(n[1]+1):(n[1]+n[2])],]
    perm.stats[i] <- get_stat(X=data.perm,n=n,d=d)$statistic
  }
  pvalue <- mean(stat <= perm.stats) #Percentage of stats from perm method that were more extreme than observed statistic

  return(list(statistic=stat, null.stats=perm.stats, pvalue=pvalue))

}


########### Function to compute statistic #############################

get_stat <- function(X,n=NULL,d=NULL){
  #X is a list of 2 dataframes
  #n is the sample size array (num samples in each dataframe of X)
  #d is the dimension (number of columns)

  #Compute p1, p2, and pooled p vector:
  pHatMatrix <- rbind(colMeans(X[[1]]), colMeans(X[[2]]), colMeans(rbind(X[[1]],X[[2]])))
  rownames(pHatMatrix) <- c("p1","p2","p")

  if(is.null(n)){
    n <- as.numeric(lapply(X,nrow)) #sample size array
  }
  if(is.null(d)){
    d <- ncol(pHatMatrix) #dimension
  }

  D.array <- apply(pHatMatrix,MARGIN=2,FUN=function(x){ if(x[3]==0 || x[3]==1){return(0)}else{(x[1]-x[2])/sqrt(x[3]*(1-x[3])*(1/n[1] + 1/n[2]))}})

  #Compute delta:
  a_d <- (log(d))^-2
  delta <- sqrt(2*log(a_d*d))

  #Comput statistic:
  stat <- (D.array^2) %*% (abs(D.array)>=delta)

  return(list(statistic=stat,pHatMatrix=pHatMatrix))
}


################### Generate binary data ############################

#' Generate multivariate binary data
#'
#' Randomly generate a list of two matrices containing multivariate binary data.
#'
#' The \eqn{(i,j)^{th}} entry of the \eqn{c^{th}} matrix is \eqn{X_{cij} = (1 - U_{ij})Y_{icj} + U_{ij}Z_{i}} where
#'
#' \itemize{
#'   \item \eqn{U_{ij} ~ Ber(r)},
#'   \item \eqn{Z_i ~ Ber(\gamma)},
#'   \item \eqn{Y_{icj} ~ Ber(p_{jc})} where
#'   \itemize{
#'     \item \eqn{p_{jc} = (1 - \beta)p_{o} + \beta h_c}
#'     \item \eqn{\beta ~ Ber(\epsilon)}
#'     \item \eqn{h_c ~ Uniform(0,\sigma_c)}
#'   }
#' }
#'
#'
#' @param n Vector of length 2 containing group size (i.e. number of samples) for
#' each group. Default value is (30,30).
#' @param d Number of variables (dimension) of the data to be generated.
#' Default value is 2000.
#' @param null_hyp Boolean indicating whether group means should be the same
#' (i.e. null hypothesis is TRUE) or different (i.e. null hypothesis is FALSE).
#' Default value is TRUE.
#' @param r Mean for distribution of of \eqn{U_{ij} ~ Ber(r)}. See details below.
#' Increase \code{r} to increase the amount of correlation among the \code{d}
#' variables. Default value is 0.3.
#' @param epsilon Used in mixture model that generates the probability vectors.
#' See details below. Sparsity can be increased by decreasing \code{epsilon}
#' and vice versa. Default value is 0.2.
#' @param sigma Used to define a uniform distribution used to generates the
#' probability vectors. See details below. Default value is (0.3,0.1).
#' @param gamma Mean for dist of \eqn{Z_i ~ Ber(gamma)}. See details below.
#' Default value is 0.3.
#' @return A list containing the following components:
#' @return X List of two n by d matrices each containing the generated datasets.
#' @return  p The probability vectors used to generate the two datasets.
#' @return null_hyp Value of the \code{null_hyp} parameter.
#' @return r Value of the \code{r} parameter.
#' @return epsilon Value of the \code{epsilon} parameter.
#'
#' @seealso
#' Amanda Plunkett & Junyong Park (2017) \emph{Two-sample tests for sparse
#' high-dimensional binary data}, Communications in Statistics - Theory and
#' Methods, 46:22, 11181-11193
#'
#' Junyong Park & J. Davis (2011) \emph{Estimating and testing conditional sums
#' of means in high dimensional multivariate binary data}, Journal of Statistical
#' Planning and Inference, 141:1021-1030
#'
#' @examples
#' binData <- genMVBinaryData(n=,d=,null_hyp=FALSE,)
#'
#' #Check the dimension of each matrix:
#' lapply(binData,dim)
#'
#' #Test whether the two datasets were generated using the same mean:
#' mvbinary.test(binData)
genMVBinaryData <- function(n=30,d=2000,null_hyp=TRUE,r=0.3,epsilon=0.2,sigma=c(0.3,0.1),gamma=0.3,p0=0.1){

  m <- length(n) #num groups

  #Generate probabilities that will be used to generate data. mxd matrix.
  if(null_hyp){
    #Using mixture dist in Dr Park paper, eqn (18). note: elementwise multiplication on purpose:
    b <- rbinom(d,size=1,prob=epsilon)
    probs <- matrix(  (rep(1,d) - b)*p0 + b * runif(d,min=0,max=sigma[1])  ,nrow=m,ncol=d,byrow=TRUE)
  }else{
    #Using mixture dist in Dr Park paper, eqn (18). elementwise multiplication on purpose:
    probs <- matrix(NA,m,d)
    for(g in 1:m){
      b <- rbinom(d,size=1,prob=epsilon)
      probs[g,] <- (rep(1,d) - b)*p0 + b * runif(d,min=0,max=sigma[g])
    }
  }

  X <- list()
  for(g in 1:m){ #loop through groups

    Y <- sapply(probs[g,],rbinom,n=n[g],size=1)  #n_g x d

    U <- matrix(rbinom(n[g]*d,size=1,prob=r),nrow=n[g],ncol=d) #n_g x d

    Z <- as.matrix(rbinom(n[g],size=1,prob=gamma)) %*% t(as.matrix(rep(1,d)))  #same for all j. uniq for each i.

    #X_ij = (1-U_{ij})Y_{ij} + U_{ij}*Z_i
    X[[g]] <- (matrix(1,n[g],d)-U)*Y + U*Z  #Note: Purposely meant to do element-wise multiplication.
  }

  return(list(X=X,p=probs,null_hyp=null_hyp,r=r,epsilon=epsilon))
}


