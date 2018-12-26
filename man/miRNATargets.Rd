\name{miRNATargets}
\alias{miRNATargets}
\docType{data}
\title{miRNA binding targets}
\description{These microRNA target candidates were predicted from 3 or more prediction database.}
\usage{data("miRNATargets")}
\format{
  The format is:
 chr "miRNATargets"
}
\details{Each number represents the gene. Corresponding gene symbols exist in 'genelist' data.}
\source{Seven miRNA target prediction database. DIANA-microt-CDS, miRanda-microCosm, mirSVR, miRTargets, TargetRank, PITA, TargetScan}
\references{Seven miRNA target prediction database. DIANA-microt-CDS, miRanda-microCosm, mirSVR, miRTargets, TargetRank, PITA, TargetScan}
\examples{
data(miRNATargets)
## maybe str(miRNATargets); ...
}
\keyword{datasets}
