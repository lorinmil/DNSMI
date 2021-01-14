
cohenkappa <- function(Sensitivity, Specificity, TPn, TNn,Tot){
	se <- Sensitivity
	sp <- Specificity
	I <- (TPn + TNn) / Tot
	if(se != sp){
		alpha <- ((se - sp)^2) / (2 * (se + sp -1))
		Kapp <-  ((I - se) * (I - sp)) / ((I - se) * (I - sp) + alpha * (I - 1))
	}else{
		s <- se
		beta <- ((TPn - TNn)^2) / (2 * TPn * TNn)
		Kapp <- (2 * s - 1) / (1 + beta * (1 - s))
	}
	
	return(Kapp)
}