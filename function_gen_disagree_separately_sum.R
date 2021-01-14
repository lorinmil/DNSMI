


gen_disagree_UandV_separate_sum <- function(c1, c2, W_list){
	
	c1 <- c1
	c2 <- c2
	W_list <- W_list
	
	uv_array <- sapply(X = W_list, FUN = function(w){
								out <- PMD(w,type = "standard",sumabs = NULL,
														sumabsu = c1,sumabsv = c2,trace = FALSE)
								cbind(u = as.vector(out$u), v = as.vector(out$v))
					},simplify = "array")
	
	
	uv_array <- ifelse(uv_array != 0, 1,0)
	u_mx <- uv_array[,"u",][1:dim(W_list[[1]])[1],]
	v_mx <- uv_array[,"v",][1:dim(W_list[[1]])[2],]
	
	urowmeans <- rowMeans(u_mx)
	uxis <- 2 * urowmeans * (1 - urowmeans)
	Duhat <- mean(uxis)
	
	
	vrowmeans <- rowMeans(v_mx)
	vxis <- 2 * vrowmeans * (1 - vrowmeans)
	Dvhat <- mean(vxis)
	
	return(list(DU = Duhat, DV = Dvhat, 
							Duindv = uxis, Dvindv = vxis,
							Duvsum = Duhat + Dvhat, 
							Duvavg = (Duhat + Dvhat) / 2,
							Duvmax = max(Duhat,Dvhat),
							Duvmin = min(Duhat,Dvhat),
							Up1vec = urowmeans, 
							Vp1vec = vrowmeans,
							Up1mean = mean(urowmeans),
							Vp1mean = mean(vrowmeans)))

}