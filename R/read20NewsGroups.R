#' Read 20 Newsgroup data
#'
#' @param g1,g2 Integers. 20 newsgroups to compare. Default: g1=10 (rec.sport.baseball), g2=14 (sci.med)
#' @param rate Decimal number in the range of [0,1]. Sampling rate at which to sample records (without replacement). Default=1.
#' @param remove_stopwords Boolean. TRUE = remove stopwords. Default=FALSE.
#' @param min_threshold Integer value. Words having count < threshold will be removed. remove words that have counts less than or equal to threshold (looking across entire group of docs)
#' @param data.file File path to the file containing the 20 newsgroup data.
#' @param labels.file File path to the file containing the 20 newsgroup labels.
#' @return Two tf vectors, one for each group.
#' @examples
#' #Generate two vectors from the same distribution:
#' data <- genMultinomialData(sample_size=1)
#'

get_20newsgroup_data <- function(g1=10,g2=14,rate=1,remove_stopwords=FALSE,min_threshold=0,data.file=NA,labels.file=NA){

  #Training data: Document ID, word ID, and count.
  #Read space delimined data
  #0 counts are suppressed

  data.train <- read.delim(data.file,header=FALSE,sep=" ")
  colnames(data.train) <- c("docIdx", "wordIdx", "count")

  words <- max(data.train$wordIdx) #totalnumber of words

  #Provides news group labels for each of 11269 documents (in training set)
  #Read space delimined data
  data.train.label <- read.delim(labels.file,header=FALSE,sep=" ")

  #Add news group labels to the training data:
  data.train <- cbind(data.train.label[data.train$docIdx,1],data.train)

  colnames(data.train) <- c("newsGrpIDx", "docIdx", "wordIdx", "count")

  #Select a couple of news groups to compare:
  grps2compare <- c(g1, g2)
  data.train <- data.train[which(data.train$newsGrpIDx == grps2compare[1] | data.train$newsGrpIDx == grps2compare[2] ),]

  #Restructure data:
  X <- list()
  docsInGrp <- list()
  for(g in 1:2){
    docsInGrp[[g]] <- unique(data.train[which(data.train[,"newsGrpIDx"]==grps2compare[g]),]$docIdx)
    if(g1 == g2){ #If data coming from the same group (to test when null is TRUE), cut it in half
      if(g==1){
        half1.Indx <- sample(1:length(docsInGrp[[g]]),replace=FALSE,size=round(0.5*length(docsInGrp[[g]])))
        docsInGrp[[g]] <- docsInGrp[[g]][half1.Indx]
      }else{
        docsInGrp[[g]] <- docsInGrp[[g]][-(half1.Indx)]
      }
    }
    docsInGrp[[g]] <- sample(docsInGrp[[g]],replace=FALSE,size=round(length(docsInGrp[[g]])*rate))
    numDocsInGrp <- length(docsInGrp[[g]])
    X[[g]] <- matrix(0,nrow=numDocsInGrp,ncol=words)

    for(i in 1:numDocsInGrp){
      docID <- docsInGrp[[g]][i]
      wordsAndCount <- data.train[which(data.train$docIdx == docID),c("wordIdx","count")]
      X[[g]][i,wordsAndCount[,"wordIdx"]] <- wordsAndCount[,"count"]
    }
  }
  rm(data.train) #needed?

  #remove stopwords
  #remove_stopwords=TRUE
  if(remove_stopwords){
    if(!require(tm)){
      install.packages("tm")
      library(tm)
    }
    words <- read_words() #Read list of words (that map to columns of X)
    cols2BRemoved <- which(noquote(stopwords("SMART")) %in% words[1:nrow(words),])
    X[[1]] <- X[[1]][,-cols2BRemoved]
    X[[2]] <- X[[2]][,-cols2BRemoved]
  }

  #Keep only columns that are greater than minimum threshold
  nonzero <- which(colSums(rbind(X[[1]],X[[2]]))>min_threshold)
  X[[1]] <- X[[1]][,nonzero]
  X[[2]] <- X[[2]][,nonzero]

  print(sprintf("Groups %d and %d read. Data dimension: %d",g1,g2,dim(X[[1]])[2]))

  return(X)

}

print_newGrp_names <- function(file=NA){
  #Read name of news groups (for each of 20 groups)
  #Read space delimined data
  data.train.map <- read.delim(file,header=FALSE,sep=" ")
  data.train.map
}

print_words <- function(){
  words <- read_words()
  words
}

read_words <- function(file=NA){
  #Read name of news groups (for each of 20 groups)
  #Read space delimined data
  words <- read.delim(file,header=FALSE,sep=" ")
  return(words)
}
