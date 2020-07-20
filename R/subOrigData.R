subOrigData = function(marker = "d2H", taxon = NULL, group = NULL, reference = NULL, 
                        age_code = NULL, mask = NULL, std_scale = "VSMOW_H",
                       niter = 1000) {
  
  #load data in funtion environ
  data("knownOrig_samples", envir = environment())
  data("knownOrig_sites", envir = environment())
  data("knownOrig_sources", envir = environment())
  data("ham", envir = environment())
  data("oam", envir = environment())
  data("hsds", envir = environment())
  data("osds", envir = environment())
  knownOrig_samples = knownOrig_samples
  knownOrig_sites = knownOrig_sites
  knownOrig_sources = knownOrig_sources
  ham = ham
  oam = oam
  hsds = hsds
  osds = osds
  
  result = knownOrig_samples

  if(length(marker) > 1){
    stop("only one marker currently allowed")
  }
  if(!marker %in% colnames(knownOrig_samples)){
    stop("marker must be column name for isotope data field")
  }
  
  if(!is.null(taxon)){
    if(!all(taxon %in% unique(knownOrig_samples$Taxon))){
      warning("One or more taxa not present in database")
    }
    result = result[result$Taxon %in% taxon,]
  }
  
  if(!is.null(group)){
    if(!all(group %in% unique(knownOrig_samples$Group))){
      warning("One or more groups not present in database")
    }
    result = result[result$Group %in% group,]
  }
  
  if(!is.null(reference)){
    if(!is.numeric(reference)){
      warning("reference is now numeric dataset ID, see knownOrig_sources.rda")
    } else if(!all(reference %in% unique(knownOrig_sources$Dataset_ID))){
      warning("One or more references not present in database")
    }
    result = result[result$Dataset_ID %in% reference,] 
  }
  
  if(!is.null(age_code)){
    if(!all(age_code %in% unique(knownOrig_samples$Age_class))){
      warning("One or more age codes not present in database")
    }
    result = result[result$Age_class %in% age_code,]
  }
  
  result = result[!is.na(result[,marker]),]
  if(nrow(result) == 0){
    stop("No samples match query")
  }

  if(!is.null(mask)) {
    if(class(mask)[1] == "SpatialPolygonsDataFrame" || 
       class(mask)[1] == "SpatialPolygons"){
      if(is.na(proj4string(mask))){
        stop("mask must have coordinate reference system")
      } else if(!identicalCRS(knownOrig_sites, mask)) {
        mask = spTransform(mask, proj4string(knownOrig_sites))
        warning("mask was reprojected")
      }
      result_sites = knownOrig_sites[mask,]
    } else {
      stop("mask should be SpatialPolygons or SpatialPolygonsDataFrame")
    }

    if(length(result_sites) > 0) {
      result = result[result$Site_ID %in% result_sites$Site_ID,]
      if(nrow(result) > 0) {
        result_sites = result_sites[result_sites$Site_ID %in% 
                                      result$Site_ID,]
        plot(mask, axes = TRUE)
        plot(result_sites, add = TRUE, col = "red")
      } else {
        stop("No samples found in mask\n")
      }
    } else {
      stop("No sites found in mask\n")
    }
    
  } else {
    result_sites = knownOrig_sites[knownOrig_sites$Site_ID %in%
                                     result$Site_ID,]
    plot(wrld_simpl, axes = TRUE)
    plot(result_sites, add = TRUE, col = "red", cex = 0.5)
    
  }
  
  message(paste(length(result[,1]),"samples are found from", 
          length(result_sites), "sites"))

  result_sources = knownOrig_sources[knownOrig_sources$Dataset_ID %in%
                                       result$Dataset_ID,]
  
  if(!is.null(std_scale)){
    if(marker == "d2H"){
      start_scales = unique(result_sources$H_std_scale) 
      am = ham
      std_table = hsds
      result_scales = merge(result, result_sources, by = "Dataset_ID", all.x = TRUE)
      result_scales = result_scales$H_std_scale
      sd_col = "d2H.sd"
    } else if(marker == "d18O"){
      start_scales = unique(result_sources$O_std_scale)
      am = oam
      std_table = osds
      result_scales = merge(result, result_sources, by = "Dataset_ID", all.x = TRUE)
      result_scales = result_scales$O_std_scale
      sd_col = "d18O.sd"
    }
    
    #Check that target exists in adj. matrix
    if(is.na(match(std_scale, rownames(am)))){
      warning("Standard scale not valid. Returning untransformed values.")
      result_data = merge(result_sites, result, by = "Site_ID", 
                          all.x = FALSE, duplicateGeoms = TRUE)
      return_obj = list("data" = result_data, "sources" = result_sources,
                        "chains" = NULL, "marker" = marker)
      class(return_obj) = "subOrigData"
    } else{
      
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
      
      std_shift = function(vals, ssv1, ssv2){
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
      
      for(i in seq_along(start_scales)){
        if(is.na(start_scales[i])){
          warning("No reference scale reported, some samples dropped from scale transformation")
          chain = NULL
        } else{
          result_sub = result[result_scales == start_scales[i],]
          result_sub = result_sub[!is.na(result_sub$Sample_ID),]
          chain = std_chain(start_scales[i], std_scale, am)
        }
        if(!is.null(chain)){
          if(length(chain) > 1){
            vals = matrix(nrow = nrow(result_sub), ncol = niter)
            for(j in seq_along(result_sub[,1])){
              vals[j,] = rnorm(niter, result_sub[j, marker], 
                               result_sub[j, sd_col])
            }
            for(j in 1:(length(chain)-1)){
              ssv1 = std_vals(chain[j], std_table)
              ssv2 = std_vals(chain[j+1], std_table)
              vals = std_shift(vals, ssv1, ssv2)
              if(is.null(vals)){
                chain = NULL
                break
              }
            }
            if(!is.null(chain)){
              result_sub[, marker] = apply(vals, 1, mean)
              result_sub[, sd_col] = apply(vals, 1, sd)
            }
          }
        }
        if(is.null(chain)){
          if(i == 1){
            result_out = result[0,]
            chain_out = list()
          }
        } else{
          if(i == 1){
            result_out = result_sub
            chain_out = list(chain)
          } else{
            result_out = rbind(result_out, result_sub)
            chain_out = append(chain_out, list(chain))
          }
        }
      }
      result_data = merge(result_sites, result_out, by = "Site_ID",
                          all.x = FALSE, duplicateGeoms = TRUE)
      row.names(result_data) = result_data$Sample_ID
      return_obj = list("data" = result_data, "sources" = result_sources,
                        "chains" = chain_out, "marker" = marker)
      class(return_obj) = "subOrigData"
      
      message(paste(length(result_data$Sample_ID),"samples from", 
                    length(unique(result_data$Site_ID)), "sites in the transformed dataset"))
    }
  } else{
    result_data = merge(result_sites, result, by = "Site_ID", 
                        all.x = FALSE, duplicateGeoms = TRUE)
    row.names(result_data) = result_data$Sample_ID
    return_obj = list("data" = result_data, "sources" = result_sources,
                      "chains" = NULL, "marker" = marker)
    class(return_obj) = "subOrigData"
  }
  
  return(return_obj)
}
