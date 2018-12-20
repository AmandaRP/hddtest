library(magrittr)

simulation <- function(data,null_hyp,delta,reps=10,num_docs=c(200,200)){

  vecs2Test <- list( matrix(NA,reps,ncol(data[[1]])), matrix(NA,reps,ncol(data[[1]])) )
  for(i in 1:reps){
    if(null_hyp){
      samples <- sample(1:nrow(data[[2]]),num_docs[1]+num_docs[2])
      group1 <- sample(1:length(samples),num_docs[1])
      vecs2Test[[1]][i,] <- data[[2]][samples[group1],] %>% colSums()
      vecs2Test[[2]][i,] <- data[[2]][samples[-group1],] %>% colSums()
    }else{
      vecs2Test[[1]][i,] <- data[[1]][sample(1:nrow(data[[1]]),num_docs[1]),] %>% colSums()
      vecs2Test[[2]][i,] <- data[[2]][sample(1:nrow(data[[2]]),num_docs[2]),] %>% colSums()
    }
  }

  result <- multinom.neighborhood.test(x=vecs2Test,delta=delta)

} #end simulation function

delta <- 1:160
p.delta.null <- simulation(data=twoNewsGroups,null_hyp=TRUE,delta=delta)$pvalue_delta
p.delta.alt  <- simulation(data=twoNewsGroups,null_hyp=FALSE,delta=delta)$pvalue_delta

par(xpd=TRUE, mar=par()$mar+c(0,0,0,5))
matplot(delta, cbind(t(p.delta.null), t(p.delta.alt)), type="l",ylab="p.delta",main="P-Value Curves for Simulation",col=c( rep("red",nrow(p.delta.null)),  rep("blue",nrow(p.delta.alt)) ) )

