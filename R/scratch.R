

myfunction <- function(data,null_hyp,delta,reps=5,num_docs=c(200,200)){

  X <- list()
  vecs2Test <- list( matrix(NA,reps,ncol(data[[1]])), matrix(NA,reps,ncol(data[[1]])) )


for(i in 1:reps){
  if(null_hyp){
    samples <- sample(1:nrow(data[[2]]),num_docs[1]+num_docs[2]) #sample row numbers for both groups
    group1 <- sample(1:length(samples),num_docs[1])
    X[[1]] <- data[[2]][samples[group1],]
    X[[2]] <- data[[2]][samples[-group1],]
  }else{
    X[[1]] <- data[[1]][sample(1:nrow(data[[1]]),num_docs[1]),]
    X[[2]] <- data[[2]][sample(1:nrow(data[[2]]),num_docs[2]),]
  }

  vecs2Test[[1]][i,] <- colSums(X[[1]])
  vecs2Test[[2]][i,] <- colSums(X[[2]])

}


result <- multinom.neighborhood.test(x=vecs2Test,delta=delta)

} #end myfunction

delta <- 1:160

p.delta.null <- myfunction(data=twoNewsGroups,null_hyp=TRUE,delta=delta)$pvalue_delta
p.delta.alt  <- myfunction(data=twoNewsGroups,null_hyp=FALSE,delta=delta)$pvalue_delta

par(xpd=TRUE, mar=par()$mar+c(0,0,0,5))
matplot(delta, cbind(t(p.delta.null), t(p.delta.alt)), type="l",ylab="p.delta",main="Title",col=c( rep("red",nrow(p.delta.null)),  rep("blue",nrow(p.delta.alt)) ) )

