\name{BADGER-package}
\alias{BADGER-package}
\alias{BADGER}
\docType{package}
\title{
A package for identifying sample mixups in combined expression and genotyping studies.
}
\description{

BADGER takes a set of expression array data and one of genotype array data and, 
given information on how sample mappings were intended to be, determines the apparent relationships 
between genotype and expression before predicting genotypes for each expression array. From 
these predictions, mismatches can be identified.
}
\details{
\tabular{ll}{
Package: \tab BADGER\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-21\cr
License: \tab GPL-2\cr
}

\code{BADGERpredSNP} takes two matrices, one expression and one genotype, the rows of which 
correspond to eQTLs. From these, it estimates the nature of the eQTL-like relationships, and
makes predictions of genotype either for the given expression matrix or an additional expression matrix
}
\author{
Lynch, AG
}
\references{
Lynch et al. "Ensuring data integrity in large modern integrative genomic studies" 
ECCB 2010, Ghent, Poster abstract. http://www.ebi.ac.uk/eccb/2010/index.php/proceedings/posters/122-other-bioinformatics-applications.html

Curtis et al. "The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups" Nature 486, 346-352

Lynch et al. "Calling sample mixups in cancer population studies" submitted
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
