\name{exampleEXP}
\alias{exampleEXP}
\docType{data}
\title{
Example expression matrix for the BADGER package.
}
\description{
A matrix of 420 expression values from Illumina WG6-V1 beadchip for 422 arrays from 90 individuals.
}
\usage{data(exampleEXP)}
\format{
  The format is:
 num [1:420, 1:422] 1.55 0.647 2.024 -0.169 -0.437 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:420] "Human-6_v1-1400551" "Human-6_v1-3940170" "Human-6_v1-1770358" "Human-6_v1-1500538" ...
  ..$ : chr [1:422] "GSM159664_NA18501_1_1" "GSM159664_NA18501_1_2" "GSM159664_NA18501_2_1" "GSM159664_NA18501_2_2" ...
}
\source{
Westra HJ et al. from their website "http://genenetwork.nl/wordpress/mixupmapper/"
}
\references{
Westra HJ et al. (2011) MixupMapper: correcting sample mix-ups in genome-wide datasets increases power to detect small genetic effects.
Bioinformatics 27: 2104-2111.
Stranger BE et al. (2007) Relative impact of nucleotide and copy number variation on gene expression phenotypes. Science (New York, NY)
315: 848-853.}
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

cbind(mycheck,mysampnames)[which(mycheck!=mysampnames),]}
\keyword{datasets}
