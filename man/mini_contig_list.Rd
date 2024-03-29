% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mini_contig_list}
\alias{mini_contig_list}
\title{Processed subset of `contig_list`}
\format{
An R `list` of `data.frame` objects
}
\usage{
data("mini_contig_list")
}
\description{
A list of 8 data frames of T cell contigs outputted from the
`filtered_contig_annotation` files, but subsetted to  365 valid T cells
which correspond to the same barcodes found in `scRep_example`. The
data is originally derived from the following
 \href{https://pubmed.ncbi.nlm.nih.gov/33622974/}{manuscript}.
}
\seealso{
\code{\link{contig_list}}
}
\concept{Data}
\keyword{datasets}
