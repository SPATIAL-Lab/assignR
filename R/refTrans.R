refTrans = function(samples, marker = "d2H", std_scale = "VSMOW_H", 
                    niter = 1000){
  data("ham", envir = environment())
  data("oam", envir = environment())
  data("hsds", envir = environment())
  data("osds", envir = environment())
  ham = ham
  oam = oam
  hsds = hsds
  osds = osds
  
  if(class(samples) == "SOD"){
    class(samples) = "data.frame"
    #Identify values based on marker
    if(marker == "d2H"){
      #Intial standard scales
      start_scales = unique(samples$H_std_scale) 
      #Pull original sample scales to a vector for later use
      samples_scales = samples$H_std_scale
      #Remove scales from samples object
      samples = samples[,-ncol(samples)]
      #Which adjacency matrix and standard table?
      am = ham
      std_table = hsds
      #Name of the relevant SD column
      sd_col = "d2H.sd"
    } else if(marker == "d18O"){
      start_scales = unique(samples$O_std_scale) 
      samples_scales = samples$O_std_scale
      samples = samples[,-ncol(samples)]
      am = oam
      std_table = osds
      sd_col = "d18O.sd"
    }
    
    #Check that target exists in adj. matrix
    if(is.na(match(std_scale, rownames(am)))){
      warning("Standard scale not valid. Returning untransformed values.")
      return(list("data" = samples, "chains" = NULL))
    } else{
      trans_out = trans(start_scales, samples_scales, samples, std_scale, 
                        am, std_table, marker, sd_col, niter)
      return(trans_out)
    }
  } else if(class(samples) == "data.frame"){
    if(!(marker %in% c("d2H", "d18O"))){
      stop("marker must be d2H or d18O")
    }
    if(!(marker %in% names(samples))){
      stop("samples must include a data field for the selected marker")
    }
    if(!(paste0(marker, ".sd") %in% names(samples))){
      stop("samples must include a sd field for the selected marker")
    }
    if(!(paste0(marker, "_std_scale") %in% names(samples))){
      stop("samples must include a standard scale field for the selected marker")
    }
    if(marker == "d2H"){
      #Vector of all starting standard scales
      start_scales = unique(samples$d2H_std_scale) 
      #Which adjacency matrix and standard table?
      am = ham
      std_table = hsds
      #Pull starting standard scales for use
      samples_scales = samples$d2H_std_scale
      #SD column
      sd_col = "d2H.sd"
    } else {
      start_scales = unique(samples$d18O_std_scale) 
      am = oam
      std_table = osds
      samples_scales = samples$d18O_std_scale
      sd_col = "d18O.sd"
    }
    if(!is.numeric(samples[,marker])){
      stop("samples data field must be numeric")
    }
    if(!is.numeric(samples[,sd_col])){
      stop("samples sd field must be numeric")
    }
    
    #Check that target exists in adj. matrix
    if(is.na(match(std_scale, rownames(am)))){
      stop("Standard scale not valid. No transformation possible.")
    } else{
      trans_out = trans(start_scales, samples_scales, samples, std_scale, 
                        am, std_table, marker, sd_col, niter)
      return(trans_out)
    } 
    
  } else{
    stop("samples must be a data.frame")
  }

}

trans = function(start_scales, samples_scales, samples, std_scale, am, 
                 std_table, marker, sd_col, niter){
  for(i in seq_along(start_scales)){
    if(is.na(start_scales[i])){
      warning("No reference scale reported, some samples dropped from scale transformation")
      chain = NULL
    } else{
      samples_sub = samples[samples_scales == start_scales[i],]
      #samples_sub = samples_sub[!is.na(samples_sub$Sample_ID),]
      chain = std_chain(start_scales[i], std_scale, am)
    }
    if(!is.null(chain)){
      if(length(chain) > 1){
        vals = matrix(nrow = nrow(samples_sub), ncol = niter)
        for(j in seq_along(samples_sub[,1])){
          vals[j,] = rnorm(niter, samples_sub[j, marker], 
                           samples_sub[j, sd_col])
        }
        for(j in 1:(length(chain)-1)){
          ssv1 = std_vals(chain[j], std_table)
          ssv2 = std_vals(chain[j+1], std_table)
          vals = std_shift(vals, ssv1, ssv2, niter)
          if(is.null(vals)){
            chain = NULL
            break
          }
        }
        if(!is.null(chain)){
          samples_sub[, marker] = apply(vals, 1, mean)
          samples_sub[, sd_col] = apply(vals, 1, sd)
        }
      }
    }
    if(is.null(chain)){
      if(i == 1){
        samples_out = samples[0,]
        chain_out = list()
      }
    } else{
      if(i == 1){
        samples_out = samples_sub
        chain_out = list(chain)
      } else{
        samples_out = rbind(samples_out, samples_sub)
        chain_out = append(chain_out, list(chain))
      }
    }
  }
  return(list("data" = samples_out, "chains" = chain_out))
}


#Breadth first algorithm to identify shortest standard chain
std_chain = function(ss1, ss2, scs){
  #Get adj. matrix row/col index for starting standard set
  node.start = match(ss1, rownames(scs))
  #This matrix will accumulate all possible paths
  chains = matrix(node.start, nrow = 1)
  #Current path length
  l = 1
  
  #Go until the target standard set is found
  while(!(ss2 %in% rownames(scs)[chains[,l]])){
    #Store the current number of paths
    srows = nrow(chains)
    #Add a step to the path
    l = l + 1
    #A new column in the matix records that step
    chains = cbind(chains, rep(NA))
    #Each row is a unique path founds so far, cycle through them
    for(i in 1:nrow(chains)){
      #If that path has already reached a dead end, skip it
      if(!is.na(chains[i, l-1])){
        #Otherwise find all standard sets that can be reached from
        #the end of the current path
        nodes.n = unname(which(scs[chains[i, l-1],] == 1))
        #For each of those possible new nodes
        for(j in nodes.n){
          #check to make sure new node is not already in chain
          if(!(j %in% chains[i,])){
            #if not, add a new row to matrix representing the
            #newly extended path...contining old path + new node
            chains = rbind(chains, c(chains[i,1:(l-1)], j))
          }
        }
      }
    }
    #If no new rows have been added then exploration of the Adj. Mat.
    #has finished without reaching the target. Exit with error.
    if(nrow(chains) == srows){
      warning(paste("Values can not be converted from", ss1, "to",
                    ss2, "standard scale. Some samples dropped from scale transformation"))
      return(NULL)
    }
  }
  #Return the chain connecting ss1 and ss2
  chain = chains[match(ss2, rownames(scs)[chains[,l]]),]
  chain = rownames(scs)[chain]
  return(chain)
}

std_vals = function(scale, sds){
  lm = sds$Low[sds$Scale == scale]
  hm = sds$High[sds$Scale == scale]
  lse = sds$Low_se[sds$Scale == scale]
  hse = sds$High_se[sds$Scale == scale]
  ref_scale = sds$Ref_scale[sds$Scale == scale]
  ssv = list("lm" = lm, "hm" = hm, "lse" = lse, "hse" = hse, 
             "ref_scale" = ref_scale)
  class(ssv) = "ssv"
  return(ssv)
}

std_shift = function(vals, ssv1, ssv2, niter){
  if(class(ssv1) != "ssv" | class(ssv2) != "ssv"){
    stop("Standard values must be provided as class ssv")
  }
  
  if(ssv1$ref_scale == ssv2$ref_scale){
    st = "m"
  } else{
    st = "s"
  }
  
  if(st == "m"){
    if(!is.na(ssv1$lse)){
      ssv1.lv = rnorm(niter, ssv1$lm, ssv1$lse)
      ssv1.hv = rnorm(niter, ssv1$hm, ssv1$hse)
      ssv2.lv = ssv1$lm
      ssv2.hv = ssv1$hm
    } else if(!is.na(ssv2$lse)){
      ssv1.lv = ssv2$lm
      ssv1.hv = ssv2$hm
      ssv2.lv = rnorm(niter, ssv2$lm, ssv2$lse)
      ssv2.hv = rnorm(niter, ssv2$hm, ssv2$hse)      
    } else{
      warning("Measurement comparison without estimate of error, some samples dropped from scale transformation")
      return(NULL)
    }
  } else{
    ssv1.lv = ssv1$lm
    ssv1.hv = ssv1$hm
    ssv2.lv = ssv2$lm
    ssv2.hv = ssv2$hm
  }
  
  m = (ssv2.hv - ssv2.lv) / (ssv1.hv - ssv1.lv)
  b = ssv2.hv - ssv1.hv * m
  
  vals.new = vals * m + b
  
  return(vals.new)
}
