BADGERpredSNP <-
function(MATEXP,MATSNP,BADGEREXP=NULL){

# If BADGEREXP is not specified, assume that we want to predict genotypes for each array in MATEXP
if(is.null(BADGEREXP)){BADGEREXP<-MATEXP}

# All inputs must be numeric 
if(!is.numeric(MATEXP)){stop("MATEXP must be a numeric matrix\n")}
if(!is.numeric(MATSNP)){stop("MATSNP must be a numeric matrix\n")}
if(!is.numeric(BADGEREXP)){stop("BADGEREXP must be a numeric matrix\n")}

# MATSNP entries must be 1,2,3 or NA
if(!all((unique(sort(MATSNP))) %in% c(1,2,3))){stop("Genotypes must be coded 1, 2 and 3\n")}

  preds<-matrix(NA,nrow=nrow(BADGEREXP),ncol=ncol(BADGEREXP))
# For each eQTL, work out the distributions for each genotype
	for(i in 1:(nrow(MATEXP))){
		if(sum(MATSNP[i,]==1,na.rm=T)>2){dens1<-density((MATEXP)[i,which(MATSNP[i,]==1)],na.rm=T)}
		if(sum(MATSNP[i,]==2,na.rm=T)>2){dens2<-density((MATEXP)[i,which(MATSNP[i,]==2)],na.rm=T)}
		if(sum(MATSNP[i,]==3,na.rm=T)>2){dens3<-density((MATEXP)[i,which(MATSNP[i,]==3)],na.rm=T)}
# For each array, make a prediction
		for(j in 1:(ncol(BADGEREXP))){
			useval<-BADGEREXP[i,j]
			if(is.na(useval)){useval<-mean(MATEXP[i,],na.rm=T)}
			pred1<-0
			pred2<-0
			pred3<-0
			if(sum(MATSNP[i,]==1,na.rm=T)>2){pred1<-dens1$y[which.min(abs(dens1$x-useval))]}
			if(sum(MATSNP[i,]==2,na.rm=T)>2){pred2<-dens2$y[which.min(abs(dens2$x-useval))]}
			if(sum(MATSNP[i,]==3,na.rm=T)>2){pred3<-dens3$y[which.min(abs(dens3$x-useval))]}
			predsum<-pred1+pred2+pred3
			pred1<-pred1/predsum
			pred2<-pred2/predsum
			pred3<-pred3/predsum
			preds[i,j]<-(pred1+2*pred2+3*pred3)
		}
	}
	colnames(preds)<-colnames(BADGEREXP)
	return(preds)
}
