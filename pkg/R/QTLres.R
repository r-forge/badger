QTLres <-
function(predictedMATSNP,MATSNP){
	matchscore<-rep(NA,nrow(MATSNP))
		for(EQTL in 1:nrow(MATSNP)){
			matchscore[EQTL]<-sum((predictedMATSNP[EQTL,]
                             -as.numeric(MATSNP[EQTL,]))^2,na.rm=T)
		}
	return(matchscore)
	}
