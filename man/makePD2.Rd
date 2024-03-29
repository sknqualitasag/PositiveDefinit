% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MakePD.R
\name{makePD2}
\alias{makePD2}
\title{Bending of Matrix A}
\usage{
makePD2(A)
}
\arguments{
\item{A}{input matrix}
}
\value{
Bended positive-definite matrix A
}
\description{
The input matrix A is decomposed into its eigen-values and eigen-vectors,
The negative eigen-values are projected into the range between zero and
the smallest positive eigen-value.
}
\details{
---
title: Check and Transform to Insure Positive Definite of Matrix
date:  "`r Sys.Date()`"
---
}
\examples{
G<- matrix(c(100,80,20,6,80,50,10,2,20,10,6,1,6,2,1,1), ncol = 4, byrow=TRUE)
makePD2(G)
}
