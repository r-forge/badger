\name{BADGERpredSNP}
\alias{BADGERpredSNP}
\title{
Function to predict genotypes from an expression matrix
}
\description{
Given two matrices (one of expression values and one of genotype calls), 
ordered by a set of known or suspected eQTLs, this function 
defines the nature of those eQTL relationship and predicts genotypes from a 
third matrix of expression values.
}
\usage{
BADGERpredSNP(MATEXP, MATSNP, BADGEREXP=NULL)
}
\arguments{
  \item{MATEXP}{
A numeric matrix of expression values, the columns of which correspond 
to the columns of \code{MATSNP} (i.e. supposedly from the same sample). 
The rows of \code{MATEXP} are believed to correspond to the rows of 
\code{MATSNP} via an eQTL. 
}
  \item{MATSNP}{
A numeric matrix of genotype calls corresponding to \code{MATEXP} as 
described above.
}
  \item{BADGEREXP}{
A numeric matrix of expression values for which corresponding genotypes 
will be predicted. If NULL, values will be predicted for the \code{MATEXP} matrix.
}
}
\value{
A matrix of mode numeric containing predicted genotypes corresponding to \code{BADGEREXP}.
}
\references{
First used/described in 

Lynch et al. "Ensuring data integrity in large modern integrative genomic studies" 
ECCB 2010, Ghent, Poster abstract. http://www.ebi.ac.uk/eccb/2010/index.php/proceedings/posters/122-other-bioinformatics-applications.html

Curtis et al. "The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups" Nature 486, 346-352

Lynch et al. "Calling sample mixups in cancer population studies" submitted

http://badger.r-forge.r-project.org/
}
\author{
Lynch AG
}

\seealso{
\code{\link{BADGER}},\code{\link{BADGERmatchScores}},\code{\link{QTLres}}
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
