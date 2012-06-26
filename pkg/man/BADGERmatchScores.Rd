\name{BADGERmatchScores}
\alias{BADGERmatchScores}
\title{
Function to compare matrices of predicted and observed genotypes.
}
\description{
This function compares the predicted genotypes from n expression arrays, 
and compares them to m observed genotype arrays.
}
\usage{
BADGERmatchScores(predictedSNP, BADGERSNP, useEQTL = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predictedSNP}{
numeric matrix of g predicted genotypes (probably from \code{BADGERpredSNP}) 
from n expression arrays.
}
  \item{BADGERSNP}{
numeric matrix of g observed genotypes from m SNP arrays.
}
  \item{useEQTL}{
which of the g genotypes should contribute? (NULL = all)
}
}
\value{
a numeric n x m matrix containing a score for the comparison of each SNP 
array with each expression array
}
\references{
First used/described in 

Lynch et al. 2010 "Ensuring data integrity in large modern integrative genomic studies" 
ECCB, Ghent, Poster abstract. http://www.ebi.ac.uk/eccb/2010/index.php/proceedings/posters/122-other-bioinformatics-applications.html

Curtis et al. 2012 "The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups" Nature 486, 346-352

Lynch et al. "Calling sample mixups in cancer population studies" submitted

http://badger.r-forge.r-project.org/
}
\author{
Lynch AG
}
\seealso{
\code{\link{BADGER}},\code{\link{BADGERpredSNP}},\code{\link{QTLres}}
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
