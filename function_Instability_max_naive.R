


InstabilityMax_naive <- function(G, X, Y, covariatesY=NULL, beta = 0.05, repeatsN = 100, level.precision=3, max.iter=1500, 
                                 distG, distX, distY,
                                 sampleWImportLoc=NULL, sampleWExportLoc=NULL){
  
    obs <- list()
    obs$X <- X
    obs$G <- G
    obs$y <- Y

		w <- W_svd_v5(df = obs, covariatesY=covariatesY, distG=distG, distX=distX, distY=distY, SVDdecomp = FALSE)$w
		######### setting #################
		n <- nrow(w)
		p <- ncol(w)
		outcomeVals <- unique(obs$y)
		subsampling1 <- floor(0.5 * sum(obs$y==outcomeVals[1]))
		subsampling2 <- floor(0.5 * sum(obs$y==outcomeVals[2]))
		repeats <- repeatsN
		
		
		#########################################
		cat("w:",dim(w)[1],"x",dim(w)[2],"\n")
		cat("Generating c1c2grid...")
		cat("Done\n")
		cat("computing Dhat list...\n")
		
		###############  ####################
		cat("generating sample W list...")
		totaltime11 <- Sys.time()
		if(!is.null(sampleWImportLoc)){
		  sampleW_list <- readRDS(sampleWImportLoc)
		}else {
		  sampleW_list <- lapply(X = seq(1,repeats), FUN = function(i){
		    # #Repeat to make sure that not all the same value within one outcome
		    # repeat{
		    #   sapl1 <- sample(x = which(obs$y==outcomeVals[1]),size = subsampling1, replace=FALSE)
		    #   if(sum(as.numeric(apply(obs$X[sapl1,], 2, FUN=function(x){length(unique(x))}))>3) &&
		    #      sum(as.numeric(apply(obs$G[sapl1,], 2, FUN=function(x){length(unique(x))}))>3)){
		    #     break
		    #   }
		    # }
		    # repeat{
		    #   sapl2 <- sample(x = which(obs$y==outcomeVals[2]),size = subsampling2, replace=FALSE)
		    #   if(sum(as.numeric(apply(obs$X[sapl2,], 2, FUN=function(x){length(unique(x))}))>3) &&
		    #      sum(as.numeric(apply(obs$G[sapl2,], 2, FUN=function(x){length(unique(x))}))>3)){
		    #     break
		    #   }
		    # }
		    sapl1 <- sample(x = which(obs$y==outcomeVals[1]),size = subsampling1, replace=FALSE)
		    sapl2 <- sample(x = which(obs$y==outcomeVals[2]),size = subsampling2, replace=FALSE)
		    sapl <- c(sapl1, sapl2)
		    sampleobs <- list(y = obs$y[sapl],
		                      X = obs$X[sapl,],
		                      G = obs$G[sapl,])
		    w <- W_svd_v5(df = sampleobs, covariatesY=covariatesY[sapl,], distG=distG, distX=distX, distY=distY, SVDdecomp = FALSE)$w
		  })
		  if (!is.null(sampleWExportLoc)){
		    saveRDS(sampleW_list, paste0(sampleWExportLoc, "sampleW_list.rda"))
		  } 
		}
		totaltime22 <- Sys.time()

		cat("Done ")
		cat("|time: ",
		round(difftime(totaltime22,totaltime11,"auto"),2),
		units(difftime(totaltime22,totaltime11,"auto")),"\n")
		#######################################

		cat("\n1. computing D...\n")
		
		#Fix v, search u
		c2_fixed <- mean(c(1, sqrt(p)))
		whichC1Cutoff <- sqrt(n)
		whichC1Cutoff_min <- 1
		c1_vec <- c(1, mean(c(1, sqrt(n))), sqrt(n))
		c1_vec_tmp <- c() #Initialize
		maxHit <- FALSE
		Du_vec <- c()
		repeat{
		  
		  #Intialize
		  c1_vec_toAdd <- c()
		  
		  #Loop through the c1 values
		  for(c1 in 1:length(c1_vec)){
		    
		    #Only process if it has not yet beed processed
		    if((!(c1_vec[c1] %in% c1_vec_tmp)) && (c1_vec[c1]<=whichC1Cutoff) && (c1_vec[c1]>=whichC1Cutoff_min)){
		      
		      cat(paste("Processing c1: ", c1_vec[c1], "..."))
		      
		      c1_vec_tmp <- c(c1_vec_tmp, c1_vec[c1])
		      
		      #Generate the stability value for U for the cutoff c1
		      Du_vec <- c(Du_vec, round(as.numeric(gen_disagree_UandV_separate_sum(c1=c1_vec[c1], c2=c2_fixed, W_list = sampleW_list)["DU"]), level.precision))
		      
		      #Sort the values by cl value
		      Du_vec <- Du_vec[order(c1_vec_tmp)]
		      c1_vec_tmp <- sort(c1_vec_tmp)
		      
		      #Check if the increase in Cv was found
		      currC1Idx <- which(c1_vec_tmp==c1_vec[c1])
		      if((currC1Idx>1) && (!maxHit) && (Du_vec[currC1Idx]-Du_vec[currC1Idx-1]>0)){
	          maxHit <- TRUE
	          whichC1Cutoff <- c1_vec[c1]
		      }
		      if(currC1Idx>1){
		        c1_vec_toAdd <- c(c1_vec_toAdd, mean(c(c1_vec_tmp[currC1Idx-1], c1_vec_tmp[currC1Idx])))
		      }
		      if(currC1Idx<length(c1_vec_tmp)){
		        c1_vec_toAdd <- c(c1_vec_toAdd, mean(c(c1_vec_tmp[currC1Idx+1], c1_vec_tmp[currC1Idx])))
		      }
		      #Break from the loop if the Du value is too big (to save computation time)
		      if(currC1Idx>1){
		        if(Du_vec[currC1Idx]>beta && Du_vec[currC1Idx-1]>beta){
		          whichC1Cutoff <- c1_vec_tmp[currC1Idx]
		          cat("Cutoff hit! Done\n")
		          break
		        }
		      }
		      
		      cat("Done\n")
		    } #End if this c1 was already processed
		  } #Loop to the next c1 value
		  
		  #Update the next batch of c1 values to check
		  c1_vec <- unique(sort(c(c1_vec, c1_vec_toAdd)))
		  
		  #Remove from c1_vec the ones that we do not need
		  c1_vec <- c1_vec[c1_vec>=whichC1Cutoff_min & c1_vec<=whichC1Cutoff]
		  
		  #Check if its time to break out of the loop (convergence)
		  if(maxHit){
		    if(sum(Du_vec==beta & c1_vec_tmp<=whichC1Cutoff)>0){
		      idx <- which(Du_vec==beta & c1_vec_tmp<=whichC1Cutoff)
		      #Make sure that the smallest index is larger than the one before it (to make sure on the correct side of the bump)
		      if(Du_vec[min(idx)]>Du_vec[min(idx)-1]){
		        selectC1 <- min(c1_vec_tmp[Du_vec==beta & c1_vec_tmp<=whichC1Cutoff]) 
		        break #Break out of the do while loop
		      }
		    }
		  }
		  
		  #If the maximum number of iterations has been hit, then exit loop
		  if((!is.null(max.iter) && length(c1_vec_tmp)>=max.iter) || (round(whichC1Cutoff, 14)==1.00000000000001)){
		    selectC1 <- c1_vec_tmp[which(Du_vec<=beta & c1_vec_tmp<=whichC1Cutoff)][which.max(Du_vec[Du_vec<=beta & c1_vec_tmp<=whichC1Cutoff])]
		    break
		  }
		  
		} #End loop for stability on U
		
		
		
		#Fix u, search v
		c1_fixed <- selectC1
		whichC2Cutoff <- sqrt(p)
		whichC2Cutoff_min <- 1
		c2_vec <- c(1, mean(c(1, sqrt(p))), sqrt(p))
		c2_vec_tmp <- c() #Initialize
		maxHit <- FALSE
		Dv_vec <- c()
		repeat{
		  
		  #Intialize
		  c2_vec_toAdd <- c()
		  
		  #Loop through the c2 values
		  for(c2 in 1:length(c2_vec)){
		    
		    #Only process if it has not yet beed processed
		    if((!(c2_vec[c2] %in% c2_vec_tmp)) && (c2_vec[c2]<=whichC2Cutoff) && (c2_vec[c2]>=whichC2Cutoff_min)){
		      
		      cat(paste("Processing c2: ", c2_vec[c2], "..."))
		      
		      c2_vec_tmp <- c(c2_vec_tmp, c2_vec[c2])
		      
		      #Generate the stability value for U for the cutoff c2
		      Dv_vec <- c(Dv_vec, round(as.numeric(gen_disagree_UandV_separate_sum(c1=c1_fixed, c2=c2_vec[c2], W_list = sampleW_list)["DV"]), level.precision))
		      
		      #Sort the values by cl value
		      Dv_vec <- Dv_vec[order(c2_vec_tmp)]
		      c2_vec_tmp <- sort(c2_vec_tmp)
		      
		      #Check if the increase in Cv was found
		      currC2Idx <- which(c2_vec_tmp==c2_vec[c2])
		      if((currC2Idx>1) && (!maxHit) && (Dv_vec[currC2Idx]-Dv_vec[currC2Idx-1]>0)){
		        maxHit <- TRUE
		        whichC2Cutoff <- c2_vec[c2]
		      }
		      if(currC2Idx>1){
		        c2_vec_toAdd <- c(c2_vec_toAdd, mean(c(c2_vec_tmp[currC2Idx-1], c2_vec_tmp[currC2Idx])))
		      }
		      if(currC2Idx<length(c2_vec_tmp)){
		        c2_vec_toAdd <- c(c2_vec_toAdd, mean(c(c2_vec_tmp[currC2Idx+1], c2_vec_tmp[currC2Idx])))
		      }
		      #Break from the loop if the Dv value is too big (to save computation time)
		      if(currC2Idx>1){
		        if(Dv_vec[currC2Idx]>beta && Dv_vec[currC2Idx-1]>beta){
		          whichC2Cutoff <- c2_vec_tmp[currC2Idx]
		          cat("Cutoff hit! Done\n")
		          break
		        }
		      }
		      
		      cat("Done\n")
		    } #End if this c2 was already processed
		  } #Loop to the next c2 value
		  
		  #Update the next batch of c2 values to check
		  c2_vec <- unique(sort(c(c2_vec, c2_vec_toAdd)))
		  
		  #Remove from c2_vec the ones that we do not need
		  c2_vec <- c2_vec[c2_vec>=whichC2Cutoff_min & c2_vec<=whichC2Cutoff]
		  
		  #Check if its time to break out of the loop (convergence)
		  if(maxHit){
		    if(sum(Dv_vec==beta & c2_vec_tmp<=whichC2Cutoff)>0){
		      idx <- which(Dv_vec==beta & c2_vec_tmp<=whichC2Cutoff)
		      #Make sure that the smallest index is larger than the one before it (to make sure on the correct side of the bump)
		      if(Dv_vec[min(idx)]>Dv_vec[min(idx)-1]){
		        selectC2 <- min(c2_vec_tmp[Dv_vec==beta & c2_vec_tmp<=whichC2Cutoff]) 
		        break #Break out of the do while loop
		      }
		    }
		  }
		  
		  #If the maximum number of iterations has been hit, then exit loop
		  if((!is.null(max.iter) && length(c2_vec_tmp)>=max.iter) || (round(whichC2Cutoff, 14)==1.00000000000001)){
		    selectC2 <- c2_vec_tmp[which(Dv_vec<=beta & c2_vec_tmp<=whichC2Cutoff)][which.max(Dv_vec[Dv_vec<=beta & c2_vec_tmp<=whichC2Cutoff])]
		    break
		  }
		  
		} #End loop for stability on V
		
		
		###############################################
		
		
		#Perform the final penalized matrix decomposition
		PMDout <- PMD(w,type = "standard",sumabs = NULL,
								sumabsu = selectC1,
								sumabsv = selectC2,
								trace = FALSE)
		
		stabilityC1 <- data.frame(c1=c1_vec_tmp, Du=Du_vec)
		stabilityC2 <- data.frame(c2=c2_vec_tmp, Dv=Dv_vec)
	
	return(list(U_spas = PMDout$u, V_spas = PMDout$v, d_spas = PMDout$d,
							opt_c1 = selectC1,
							opt_c2 = selectC2,
							stabilityC1 = stabilityC1,
							stabilityC2 = stabilityC2,
							w = w,
							inputVals = list(G=G, X=X, Y=Y, beta=beta, repeatsN=repeatsN, level.precision=level.precision, max.iter=max.iter),
							sampleW_list = sampleW_list))

}

