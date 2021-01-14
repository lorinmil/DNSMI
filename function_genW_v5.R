

W_svd_v5 <- function(df, covariatesY, distG, distX, distY, SVDdecomp = FALSE){
	
	ngenes <- ncol(df$G)
	ntranscripts <- ncol(df$X)
	
	#LM: UPDATE THIS!
	if(distG=="poisson"){
	  betayg_vec <- apply(X=xDF, MARGIN=2, FUN=function(x){
	    if(!is.null(covariatesY)){
	      gDF <- cbind(x=log(x+1), covariatesY)
	    }else {
	      gDF <- log(x+1)
	    }
	    OLS_cross(x=gDF, y = df$y, family_ = distY)[1]
	  
	  })
	}else {
	  betayg_vec <- apply(X=df$G, MARGIN=2, FUN=function(x){
                      	    if(!is.null(covariatesY)){
                      	      gDF <- cbind(x=x, covariatesY)
                      	    }else {
                      	      gDF <- x
                      	    }
	                          OLS_cross(x=gDF, y = df$y, family_ = distY)[1]
	                 })
	}
 
  #LM: UPDATE THIS!
	if(distG=="poisson" && distX=="poisson"){
	  betaxg_mx <- t(apply(X = df$G+1, MARGIN = 2,FUN = function(gene) {
	    apply(X = df$X+1, MARGIN=2, FUN=OLS_cross, x=log(gene), family_="poisson")
	  }))
	}else if(distG=="poisson" && distX!="poisson"){
	  betaxg_mx <- t(apply(X = df$G+1, MARGIN = 2,FUN = function(gene) {
	    apply(X = df$X, MARGIN=2, FUN=OLS_cross, x=log(gene), family_=distX)
	  }))
	}else if(distG!="poisson" && distX=="poisson"){
	  betaxg_mx <- t(apply(X = df$G, MARGIN = 2,FUN = function(gene) {
	    apply(X = df$X+1, MARGIN=2, FUN=OLS_cross, x=gene, family_="poisson")
	  }))
	}else {
	  betaxg_mx <- t(apply(X = df$G, MARGIN = 2,FUN = function(gene) {
	    apply(X = df$X, MARGIN=2, FUN=OLS_cross, x=gene, family_=distX)
	  }))
	}
  
  
	indexgrid <- as.matrix(expand.grid(
																			geneindex = seq(1,ngenes),
																			xindex = seq(1,ntranscripts)))
	indexgrid <- indexgrid[,c(2,1)]
	
	betay_xg_vec <- apply(indexgrid,1,function(cc){
									xg <- data.frame(
															var1=df$X[,cc[1],drop=FALSE],
															var2=df$G[,cc[2],drop=FALSE])
									if(!is.null(covariatesY)){
									  xg <- cbind(xg, covariatesY)
									}
									if(distX=="poisson"){
									  xg[,1] <- log(xg[,1]+1)
									}
									if(distG=="poisson"){
									  xg[,2] <- log(xg[,2]+1)
									}
									coeffs <- OLS_cross(x=xg, y=df$y, family_=distY)
									coeffs[1:2] #Exclude any covariates (if any)
									})
  
	
  betay_xg_mx <- matrix(betay_xg_vec,nrow = ngenes,ncol = ntranscripts,byrow = FALSE)
  
  
  
  betayg_mx <- matrix(rep(betayg_vec,each = ntranscripts),nrow = ngenes, ncol = ntranscripts,
											byrow = TRUE)
											
											
  w <- betayg_mx * betaxg_mx * betay_xg_mx
  
  w <- ifelse(w < 0, 0, ifelse(w > betayg_mx^2, 0,w))
  
  alphamx <- w / betayg_mx^2

  
  if(SVDdecomp == TRUE){
  	singulars <- svd(w)$d
  	u <- svd(w)$u
  	v <- svd(w)$v
  	return(list(w=w,alpha=alphamx,singulars=singulars,u=u,v=v))
  }else{
  	return(list(w=w,alpha=alphamx))
  }
}