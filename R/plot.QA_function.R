plot.QA = function(obj, savePDF = FALSE){

  if(class(obj) == "QA"){
    n=1
  } else if (class(obj) == "list"){
    n = length(obj)
    for(i in 1:n){
      if(class(obj[[i]]) != "QA"){
        stop("obj must be either a single QA object or a list of QA objects")
      }
    }
  } else {
    stop("obj must be either a single QA object or a list of QA objects")
  }
  
  xx <- seq(0.01, 0.99, 0.01)
  
  if(n == 1){
    vali = ncol(obj$val_stations)
    niter = nrow(obj$val_stations)
    
    means.p = data.frame(xx, apply(obj$prption_byProb, 2, mean)/vali)
    means.a = data.frame(xx, apply(obj$prption_byArea, 2, mean)/vali)
    
    precision = matrix(rep(0, niter*99), ncol=niter, nrow=99)
    for (i in 1:niter){
      precision[,i] <- apply(obj$precision[[i]],1, median)
    }
    
    mean = NULL
    for(i in 1:99){
      mean = append(mean, mean(precision[i,]))
    }
    
    pre <- data.frame(xx,  1-mean)
    
    pd <- data.frame(as.numeric(obj$pd_bird_val) / obj$random_prob_density)
    
    plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2,
         xlab="Cumulative probability threshold", 
         ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
    lines(pre[,1], pre[,2], lwd=2)

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2,
         xlab="Cumulative probability threshold", 
         ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
    lines(means.p[,1], means.p[,2], lwd=2)

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, 
         xlab="Cumulative area threshold", 
         ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
    lines(means.a[,1], means.a[,2], lwd=2)
    
    boxplot(pd, ylab = "Odds ratio (known origin:random)", outline = FALSE)
    
    if(savePDF){
      pdf("QA1.pdf", width = 8, height = 3)
      layout(matrix(c(1,2,3), ncol=3))
      
      plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2,
           xlab="Cumulative probability threshold", 
           ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
      lines(pre[,1], pre[,2], lwd=2)

      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2,
           xlab="Cumulative probability threshold", 
           ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
      lines(means.p[,1], means.p[,2], lwd=2)

      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, 
           xlab="Cumulative area threshold", 
           ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
      lines(means.a[,1], means.a[,2], lwd=2)
      dev.off()
      
      pdf("QA2.pdf", width=8, height=6)      
      boxplot(pd, ylab = "Odds ratio (known origin:random)", outline = FALSE)
      dev.off()
    }
    
  } else {

    nm = rep("", n)
    vali = niter = rep(0, n)
    
    for(i in 1:n){
      if(is.null(names(obj[i]))){
        nm[i] = as.character(i)
      } else if(names(obj[i]) == "") {
        nm[i] = as.character(i)
      } else {
        nm[i] = names(obj[i])
      } 
      vali[i] = ncol(obj[[i]]$val_stations)
      niter[i] = nrow(obj[[i]]$val_stations)
    }
    
    means.p = data.frame(xx, apply(obj[[1]]$prption_byProb, 2, mean)/vali[1])
    means.a = data.frame(xx, apply(obj[[1]]$prption_byArea, 2, mean)/vali[1])
    
    precision = matrix(rep(0, niter[1] * 99), ncol=niter[1], nrow=99)
    for (i in 1:niter[1]){
      precision[,i] <- apply(obj[[1]]$precision[[i]],1, median)
    }
    
    mean = NULL
    for(i in 1:99){
      mean = append(mean, mean(precision[i,]))
    }
    
    pre = data.frame(xx,  1-mean)
    
    pd = matrix(rep(NA, n * max(niter) * max(vali)), ncol=n)
    pd[1:(niter[1]*vali[1]),1] = as.numeric(obj[[1]]$pd_bird_val) / obj[[1]]$random_prob_density
    
    for(i in 2:n){

      means.p = cbind(means.p, apply(obj[[i]]$prption_byProb, 2, mean)/vali[i])
      means.a = cbind(means.a, apply(obj[[i]]$prption_byArea, 2, mean)/vali[i])
      
      precision = matrix(rep(0, niter[i] * 99), ncol=niter[i], nrow=99)
      for (j in 1:niter[i]){
        precision[,j] <- apply(obj[[i]]$precision[[j]],1, median)
      }
      
      mean = NULL
      for(j in 1:99){
        mean = append(mean, mean(precision[j,]))
      }
      
      pre = cbind(pre,  1-mean)
      
      pd[1:(niter[1]*vali[1]),i] = as.numeric(obj[[i]]$pd_bird_val) / obj[[i]]$random_prob_density
      
    }
    
    plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2,
         xlab="Cumulative probability threshold", 
         ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
    for(i in 1:n){
      lines(pre[,1], pre[,i+1], lwd=2, col=i+1)
    }
    legend(0.01, 0.55, nm, lw=2, col=seq(2,n+1), bty="n")

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2,
         xlab="Cumulative probability threshold", 
         ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
    for(i in 1:n){
      lines(means.p[,1], means.p[,i+1], lwd=2, col=i+1)
    }
    legend(0.01, 1, nm, lw=2, col=seq(2,n+1), bty="n")
      
    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, 
         xlab="Cumulative area threshold", 
         ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
    for(i in 1:n){
      lines(means.a[,1], means.a[,i+1], lwd=2, col=i+1)
    }
    legend(0.6, 0.55, nm, lw=2, col=seq(2,n+1), bty="n")  
    
    boxplot(pd, col=seq(2,n+1), names = nm,
            ylab = "Odds ratio (known origin:random)", outline = FALSE)
    
    if(savePDF){
      pdf("QA1.pdf", width = 8, height = 3)
      layout(matrix(c(1,2,3), ncol=3))
      
      plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2,
           xlab="Cumulative probability threshold", 
           ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
      for(i in 1:n){
        lines(pre[,1], pre[,i+1], lwd=2, col=i+1)
      }
      legend(0.01, 0.55, nm, lw=2, col=seq(2,n+1), bty="n")
      
      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2,
           xlab="Cumulative probability threshold", 
           ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
      for(i in 1:n){
        lines(means.p[,1], means.p[,i+1], lwd=2, col=i+1)
      }
      legend(0.01, 1, nm, lw=2, col=seq(2,n+1), bty="n")
      
      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, 
           xlab="Cumulative area threshold", 
           ylab="Proportion of validation stations included", xlim=c(0,1), ylim=c(0,1))
      for(i in 1:n){
        lines(means.a[,1], means.a[,i+1], lwd=2, col=i+1)
      }
      legend(0.6, 0.55, nm, lw=2, col=seq(2,n+1), bty="n")        
      dev.off()
      
      pdf("QA2.pdf", width=8, height=6)      
      boxplot(pd, col=seq(2,n+1), names = nm,
              ylab = "Odds ratio (known origin:random)", outline = FALSE)
      dev.off()
    }
    
    return()
  }
}
  

