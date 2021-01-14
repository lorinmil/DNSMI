


kappafrPMD <- function(PMDoutvec,trueposinvec){

		SpsU <- PMDoutvec
	
		testU_pos <- rep(0,length(SpsU))
   	testU_pos[which(SpsU != 0)] <- 1
   	
   	trueU_pos <- rep(0,length(SpsU))
   	trueU_pos[trueposinvec] <- 1
   	
   	predU <- prediction(predictions = testU_pos,labels = trueU_pos)
   	prU <- performance(prediction.obj = predU, measure = "prec", x.measure = "rec")
   	ssU <- performance(prediction.obj = predU, measure = "sens", x.measure = "spec")
   	Usensi <- ssU@y.values[[1]][2]
   	Uspeci <- ssU@x.values[[1]][2]
   	Uprec <- prU@y.values[[1]][2]
   	Urecal <- prU@x.values[[1]][2]
   	UKappa <- cohenkappa(Sensitivity = Usensi,Specificity = Uspeci,
   											TPn = Usensi * length(trueposinvec),
   											TNn = Uspeci * (length(PMDoutvec) - length(trueposinvec)),
   											Tot = length(PMDoutvec))
   											
   return(list(Kappa = UKappa,Sensi = Usensi,Speci = Uspeci))
}