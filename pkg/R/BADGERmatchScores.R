BADGERmatchScores <-
function(predictedSNP,BADGERSNP,useEQTL=NULL){
	if(is.null(useEQTL)){useEQTL<-1:(nrow(BADGERSNP))}
	matchscore<-matrix(NA,nrow=ncol(predictedSNP),ncol=ncol(BADGERSNP))
		for(pred in 1:ncol(predictedSNP)){
			for(obs in 1:ncol(BADGERSNP)){
				matchscore[pred,obs]<-sum((predictedSNP[useEQTL,pred]-
				as.numeric(BADGERSNP[useEQTL,obs]))^2,na.rm=T)
			}
		}
	return(matchscore)
	}
