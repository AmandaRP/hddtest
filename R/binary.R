#' Test two multivariate binary datasets
#'
#' Peforms the test for two binary vectors
#' testing \eqn{H_0:} the underlying probability vectors are the
#' same vs. \eqn{H_1:} they are different.
#'
#'
#' Statistic:
#' \eqn{T = \sum_{j=1}^d D_j^2 I( |Dj| \ge \delta(d))}
#' where \eqn{d} is the dimension.
#' \eqn{Dj = (\hat{p}_{1j} − \hat{p}_{2j} )/sqrt{ \hat{p}_j (1 − \hat{p}_j )(1/n1 + 1/n2) } }
#' where \eqn{\hat{p}_{cj}} is the estimate of \eqn{p_{cj}} for the \eqn{c^{th}}
#' group calculated by  \eqn{x_{cj}}. Let \eqn{\hat{p}_j} be the pooled estimate
#' for the \eqn{j^{th}} variable. Additionally,
#' \eqn{\delta(d) = \sqrt{2 log (a_d d)}} where \eqn{a_d = (log d)^{-2}}.
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

  return(list(statistic=stat, perm.stats=perm.stats, pvalue=pvalue))

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

#Fuction generates nxd matrix of random correlated binary data based on methed used by Dr Park.
# n = vector containing group size for each group
# d = num variables (dimension)
# same: boolean T/F indicating whether group means are the same or not for each group.
# r: mean for dist of U_{ij} ~ Ber(r)
# gamma: mean for dist of z_i ~ Ber(gamma)
# epsilon: Used in mixture model that generates the probabilities. Based on (18) in Dr Park paper.
# sigma: used to define uniform distribution that generates the probabilities. Based on (18) in Dr Park paper. length=num_groups
genMVBinaryData <- function(n,d,same=TRUE,r=0.3,gamma=0.3,epsilon=0.2,sigma=c(0.3,0.1),p0=0.1){

  m <- length(n) #num groups

  #Generate probabilities that will be used to generate data. mxd matrix.
  if(same){
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

  return(list(X=X,p=probs,same=same,r=r,gamma=gamma,n=n,d=d))
}


