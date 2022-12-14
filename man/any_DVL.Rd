% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/any_DVL.R
\name{any_DVL}
\alias{any_DVL}
\title{Summarize deep sequenced data into overall infection status by well}
\usage{
any_DVL(assay)
}
\arguments{
\item{assay}{Assay data, with rows representing the distinct viral lineages (DVL) and columns representing the wells.}
}
\value{
A vector of length \code{nrow(assay)}
}
\description{
Summarize deep sequenced data into overall infection status by well
}
