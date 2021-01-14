
OLS <- function(x,y){
	tx <- t(x)
	#LM: THIS IS WHERE IT WILL NEED TO BE UPDATED!!
	mylogit <- glm(y ~ x, family = "binomial")
	mylogit$coefficients["x"]
	# solve(tx %*% x) %*% tx %*% y
}

OLS_qr <- function(x,y){
	solve(qr(x), y)
}




#LM: I think this is the only one necessary
OLS_cross <- function(x,y,family_="gaussian"){
  
  #TODO: MAKE THIS MORE ROBUST TO MULTIPLE DATA TYPES
  if(is.data.frame(x)){
    x <- as.matrix(x)
  }

	out <- try(glm(y ~ x, family = family_)$coefficients)
	
	if(class(out) == "try-error"){
		return(matrix(0,ncol(x),1))
	} else{
	  out <- out[2:length(out)] #Exclude the intercept
	  
	  return(out)
	 }
}