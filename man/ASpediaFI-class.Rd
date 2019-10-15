\name{ASpediaFI-class}
\docType{class}
\alias{ASpediaFI-class}
\alias{ASpediaFI}
\alias{samples,ASpediaFI-method}
\alias{samples<-,ASpediaFI-method}
\alias{events,ASpediaFI-method}
\alias{events<-,ASpediaFI-method}
\alias{psi,ASpediaFI-method}
\alias{psi<-,ASpediaFI-method}
\alias{network,ASpediaFI-method}
\alias{gene.table,ASpediaFI-method}
\alias{as.table,ASpediaFI-method}
\alias{pathway.table,ASpediaFI-method}
\alias{samples}
\alias{samples<-}
\alias{events}
\alias{events<-}
\alias{psi}
\alias{psi<-}
\alias{network}
\alias{gene.table}
\alias{as.table}
\alias{pathway.table}

\title{ASpediaFI class}
\description{
S4 class for functional interaction analysis of AS events
}
\section{Constructor}{
\describe{
\item{}{
\code{ASpediaFI(sample.names, bam.files, conditions)}:
the constructor function requires sample names (IDs), paths to RNA-Seq BAM
files, and sample conditions.
}
}
}
\section{Accessors}{
In the code snippets below, \code{x} is an ASpediaFI object.
\describe{
\item{}{
\code{samples(x)}, \code{samples(x) <- value}:
get or set sample information. \code{value} must be a data frame.
}
\item{}{
\code{events(x)}, \code{events(x) <- value}:
get or set AS event annotations. \code{value} must be a list.}
\item{}{
\code{psi(x)}, \code{psi(x) <- value}:
get or set PSI values. \code{value} must be a SummarizedExperiment
object.}
\item{}{
\code{network(x)}, \code{gene.table(x)}, \code{as.table(x)},
\code{pathway.table(x)}:
get results of functional interaction analysis.}}
}

\examples{
bamMT <- system.file("extdata/GSM3167287.subset.bam", package = "ASpediaFI")
bamWT <- system.file("extdata/GSM3167290.subset.bam", package = "ASpediaFI")
GSE114922.ASpediaFI <- ASpediaFI(sample.names = c("GSM3167287", "GSM3167290"),
                                 bam.files = c(bamMT, bamWT),
                                 conditions = c("MT", "WT"))
}
\keyword{classes}
