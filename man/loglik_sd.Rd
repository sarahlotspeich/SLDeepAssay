% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik_sd.R
\name{loglik_sd}
\alias{loglik_sd}
\title{Log-likelihood for IUPM from single-dilution assay data}
\usage{
loglik_sd(l, M, MP, m, Y)
}
\arguments{
\item{l}{Vector of DVL-specific parameters.}

\item{M}{Total number of wells originally sequenced with the QVOA.}

\item{MP}{Number of p24-positive wells.}

\item{m}{Number of p24-positive wells that underwent the UDSA.}

\item{y}{A vector of DVL-specific infection counts.}
}
\value{
A scalar
}
\description{
Log-likelihood for IUPM from single-dilution assay data
}
