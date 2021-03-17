\name{FlowcytometricData}
\alias{FlowcytometricData}

\docType{data}
\title{Example data of hematologic marker expression.}

\description{Data set of 6 flow cytometry-based lymphoma makers from 55,843 cells from healthy subjects (class 1) and 55,843 cells from lymphoma patients (class 2).}

\usage{data("FlowcytometricData")}

\details{
Size 111686 x 6 , stored in \code{FlowcytometricData$[Var_1,Var_2,Var_3,Var_4,Var_5,Var_6]}
Classes 2, stored in \code{FlowcytometricData$Cls}
}

\examples{
data(FlowcytometricData)
str(FlowcytometricData)
}
