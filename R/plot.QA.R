plot.QA = function(x, ..., outDir = NULL){

  a = list(x, ...)

  if(class(a[[1]])[1] != "QA"){
    stop("x must be one or more QA objects")
  }
  if(!is.null(outDir)){
    if(class(outDir)[1] != "character"){
      stop("outDir should be a character string")
    }
    if(!dir.exists(outDir)){
      warning("outDir does not exist, creating")
      dir.create(outDir)
    }
  }
  
  n = 0
  bys = integer()
  for(i in seq_len(length(a))){
    if(class(a[[i]])[1] == "QA"){
      n = n + 1
      if(is.null(a[[i]]$by)){
        stop("plot now requires QA objects with the by element (use assignR v.1.2.1 or greater)")
      }
      bys = c(bys, a[[i]]$by)
    } 
  }
  
  if(n > 1){
    if(!all(diff(bys) == 0)){
      stop("plotting multiple QA objects requires that all have the same by increment")
    }
  }
  
  #vector of thresholds
  xx = seq(0.00, 1, a[[1]]$by/100)
  if(tail(xx, 1) != 1){
    xx = c(xx, 1)
  }
  
  if(n == 1){
    
    vali = ncol(x$val_stations)
    niter = nrow(x$val_stations)
    
    means.p = data.frame(xx, apply(x$prption_byProb, 2, mean))
    means.a = data.frame(xx, apply(x$prption_byArea, 2, mean))
    
    precision = matrix(rep(0, niter*length(xx)), 
                       ncol=niter, nrow=length(xx))
    for (i in 1:niter){
      precision[,i] = apply(x$precision[[i]], 1, median)
    }
    
    mean.pre = NULL
    for(i in 1:length(xx)){
      mean.pre = append(mean.pre, mean(precision[i,]))
    }
    
    pre = data.frame(xx, 1 - mean.pre)

    pd = data.frame(as.numeric(x$pd_val) / x$random_prob_density)
    
    plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2, lty=3,
         xlab="Probability quantile", 
         ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
    lines(pre[,1], pre[,2], lwd=2)

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
         xlab="Probability quantile", 
         ylab="Proportion of validation stations included", xlim=c(0,1), 
         ylim=c(0,1))
    lines(means.p[,1], means.p[,2], lwd=2)

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
         xlab="Area quantile", 
         ylab="Proportion of validation stations included", xlim=c(0,1), 
         ylim=c(0,1))
    lines(means.a[,1], means.a[,2], lwd=2)
    
    boxplot(pd, ylab = "Odds ratio (known origin:random)", 
            outline = FALSE)
    abline(1,0, col="dark grey", lwd=2, lty=3)
    
    if(!is.null(outDir)){
      
      png(paste0(outDir, "/QA1.png"), units = "in", width = 8, 
          height = 3, res = 600)
      p = par(no.readonly = TRUE)
      on.exit(par(p))
      par(mfrow = c(1, 3))
      
      plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2, lty=3,
           xlab="Probability quantile", 
           ylab="Proportion of area excluded", xlim=c(0,1), 
           ylim=c(0,1))
      lines(pre[,1], pre[,2], lwd=2)
      
      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
           xlab="Probability quantile", 
           ylab="Proportion of validation stations included", 
           xlim=c(0,1), ylim=c(0,1))
      lines(means.p[,1], means.p[,2], lwd=2)
      
      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
           xlab="Area quantile", 
           ylab="Proportion of validation stations included", 
           xlim=c(0,1), ylim=c(0,1))
      lines(means.a[,1], means.a[,2], lwd=2)
      
      dev.off()
      
      png(paste0(outDir, "/QA2.png"), units = "in", width = 6, 
          height = 4, res = 600)
      
      boxplot(pd, ylab = "Odds ratio (known origin:random)", 
              outline = FALSE)
      abline(1,0, col="dark grey", lwd=2, lty=3)
      
      dev.off()
    }
    
  } else{

    nm = rep("", n)
    vali = niter = rep(0, n)
    
    for(i in seq_len(n)){
      if(is.null(a[[i]]$name)){
        nm[i] = as.character(i)
      } else if(a[[i]]$name == "") {
        nm[i] = as.character(i)
      } else {
        nm[i] = a[[i]]$name
      } 
      vali[i] = ncol(a[[i]]$val_stations)
      niter[i] = nrow(a[[i]]$val_stations)
    }
    
    means.p = data.frame(xx, apply(a[[1]]$prption_byProb, 2, 
                                   mean))
    means.a = data.frame(xx, apply(a[[1]]$prption_byArea, 2, 
                                   mean))
    
    precision = matrix(ncol=niter[1], nrow=length(xx))
    for (i in 1:niter[1]){
      precision[,i] = apply(a[[1]]$precision[[i]],1, median)
    }
    
    mean.pre = NULL
    for(i in seq_along(xx)){
      mean.pre = append(mean.pre, mean(precision[i,]))
    }
    
    pre = data.frame(xx, 1 - mean.pre)
    
    pd = matrix(ncol=n, nrow = max(niter) * max(vali))
    pd[1:(niter[1] * vali[1]), 1] = as.numeric(a[[1]]$pd_val) / 
      a[[1]]$random_prob_density
    
    for(i in seq_len(n)[-1]){

      means.p = cbind(means.p, apply(a[[i]]$prption_byProb, 2, 
                                     mean))
      means.a = cbind(means.a, apply(a[[i]]$prption_byArea, 2, 
                                     mean))
      
      precision = matrix(ncol=niter[i], nrow=length(xx))
      for (j in seq(niter[i])){
        precision[,j] = apply(a[[i]]$precision[[j]], 1, median)
      }
      
      mean.pre = NULL
      for(j in seq_along(xx)){
        mean.pre = append(mean.pre, mean(precision[j,]))
      }
      
      pre = cbind(pre, 1 - mean.pre)
      
      pd[1:(niter[i] * vali[i]), i] = as.numeric(a[[i]]$pd_val) / 
        a[[i]]$random_prob_density
      
    }
    
    plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2, lty=3,
         xlab="Probability quantile", 
         ylab="Proportion of area excluded", xlim=c(0,1), ylim=c(0,1))
    for(i in seq_len(n)){
      lines(pre[,1], pre[,i+1], lwd=2, col=i+1)
    }
    legend(0.01, 0.55, nm, lwd=2, col=seq(2,n+1), bty="n")

    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
         xlab="Probability quantile", 
         ylab="Proportion of validation stations included", 
         xlim=c(0,1), ylim=c(0,1))
    for(i in seq_len(n)){
      lines(means.p[,1], means.p[,i+1], lwd=2, col=i+1)
    }
    legend(0.01, 1, nm, lwd=2, col=seq(2,n+1), bty="n")
      
    plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3, 
         xlab="Area quantile", 
         ylab="Proportion of validation stations included", 
         xlim=c(0,1), ylim=c(0,1))
    for(i in seq_len(n)){
      lines(means.a[,1], means.a[,i+1], lwd=2, col=i+1)
    }
    legend(0.6, 0.55, nm, lwd=2, col=seq(2,n+1), bty="n")  
    
    boxplot(pd, col=seq(2,n+1), names = nm,
            ylab = "Odds ratio (known origin:random)", outline = FALSE)
    abline(1, 0, col="dark grey", lwd=2, lty=3)
    
    if(!is.null(outDir)){
      png(paste0(outDir, "/QA1.png"), units = "in", width = 8, 
          height = 3, res = 600)
      p = par(no.readonly = TRUE)
      on.exit(par(p))
      par(mfrow = c(1, 3))      
      
      plot(c(0,1), c(1,0), type="l", col="dark grey", lwd=2, lty=3,
           xlab="Probability quantile", 
           ylab="Proportion of area excluded", xlim=c(0,1), 
           ylim=c(0,1))
      for(i in seq_len(n)){
        lines(pre[,1], pre[,i+1], lwd=2, col=i+1)
      }
      t = par("usr")[4] * 0.6
      b = par("usr")[3]
      yp1 = mean(c(t,b))
      yp = yp1 + (n-1) / 6 * yp1 
      legend(0, yp, nm, lwd=2, col=seq(2,n+1), bty="n")
      text(0.95, 0.95, "(a)")
      
      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3, 
           xlab="Probability quantile", 
           ylab="Proportion of validation stations included", 
           xlim=c(0,1), ylim=c(0,1))
      for(i in seq_len(n)){
        lines(means.p[,1], means.p[,i+1], lwd=2, col=i+1)
      }
      text(0.05, 0.95, "(b)")

      plot(c(0,1), c(0,1), type="l", col="dark grey", lwd=2, lty=3,
           xlab="Area quantile", 
           ylab="Proportion of validation stations included", 
           xlim=c(0,1), ylim=c(0,1))
      for(i in seq_len(n)){
        lines(means.a[,1], means.a[,i+1], lwd=2, col=i+1)
      }
      text(0.05, 0.95, "(c)")
      
      dev.off()

      png(paste0(outDir, "/QA2.png"), units = "in", width = 6, 
          height = 4, res = 600)

      boxplot(pd, col=seq(2,n+1), names = nm,
              ylab = "Odds ratio (known origin:random)", 
              outline = FALSE)
      abline(1, 0, col="dark grey", lwd=2, lty=3)
      
      dev.off()
    }
    
    return()
  }
}
  

