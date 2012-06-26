\name{exampleSNP}
\alias{exampleSNP}
\docType{data}
\title{
Example genotype matrix for the BADGER package.
}
\description{
A matrix of 420 genotype callsfor 81 individuals.
}
\usage{data(exampleSNP)}
\format{
  The format is:
 num [1:420, 1:81] 1 1 1 1 1 2 2 2 2 2 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:81] "NA18500" "NA18501" "NA18503" "NA18504" ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://hapmap.ncbi.nlm.nih.gov/biomart/martview/
}
\references{
Westra HJ et al. (2011) MixupMapper: correcting sample mix-ups in genome-wide datasets increases power to detect small genetic effects.
Bioinformatics 27: 2104-2111.
Stranger BE et al. (2007) Relative impact of nucleotide and copy number variation on gene expression phenotypes. Science (New York, NY)
315: 848-853.
}
\examples{
data(exampleEXP)
data(exampleSNP)
mysampnames<-sapply(colnames(exampleEXP),function(x){unlist(strsplit(as.character(x),"_"))[2]})
mypred<-BADGERpredSNP(exampleEXP[,match(colnames(exampleSNP),mysampnames)],exampleSNP,exampleEXP)
myQTLres<-QTLres(mypred[,match(colnames(exampleSNP),mysampnames)],exampleSNP)
myMS<-BADGERmatchScores(mypred,exampleSNP,which(myQTLres<35))

mycheck<-rep(0,422) 
for(i in 1:422){mycheck[i]<-colnames(exampleSNP)[which.min(myMS[i,])]} 
mycheck[!(mysampnames \%in\% colnames(exampleSNP))]<-NA

# see that there is a problem with the expression matrices associated with NA18515

cbind(mycheck,mysampnames)[which(mycheck!=mysampnames),]
}
\keyword{datasets}
