\name{QTLres}
\alias{QTLres}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to assess the performance of the eQTLs that BADGER uses
}
\description{
This function compares the predicted and observed genotypes within BADGER 
with a view to informing refinement of the eQTL list
}
\usage{
QTLres(predictedMATSNP, MATSNP)
}
\arguments{
  \item{predictedMATSNP}{
A numeric matrix of predicted genotypes (probably obtained from \code{BADGERpredSNP}
}
  \item{MATSNP}{
A matching numeric matrix of observed genotypes.
}
}
\details{
Based on Euclidean distance.
}
\value{
Numeric vector, each entry giving a measure of the predictive ability of an eQTL
}
\references{
First used/described in 

Lynch et al. "Ensuring data integrity in large modern integrative genomic studies" 
ECCB 2010, Ghent, Poster abstract. http://www.ebi.ac.uk/eccb/2010/index.php/proceedings/posters/122-other-bioinformatics-applications.html

Curtis et al. "The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups" Nature 486, 346-352

Lynch et al. "Calling sample mixups in cancer population studies" submitted

http://badger.r-forge.r-project.org/%% ~put references to the literature/web site here ~
}
\author{
Lynch AG
}

\seealso{
\code{\link{BADGER}},\code{\link{BADGERmatchScores}},\code{\link{BADGERpredSNP}}
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
