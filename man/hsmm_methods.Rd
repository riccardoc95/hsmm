\name{hsmm-methods}
\alias{print.hsmm_fit}
\alias{summary.hsmm_fit}
\alias{print.summary.hsmm_fit}
\alias{predict.hsmm_fit}
\alias{plot.hsmm_fit}
\alias{print.hsmm_simulation}
\alias{plot.hsmm_simulation}
\title{Methods for fitted and simulated HSMM objects}
\usage{
\method{print}{hsmm_fit}(x, ...)
\method{summary}{hsmm_fit}(object, ...)
\method{print}{summary.hsmm_fit}(x, ...)
\method{predict}{hsmm_fit}(object, type = c("state", "posterior"), ...)
\method{plot}{hsmm_fit}(x, ...)
\method{print}{hsmm_simulation}(x, ...)
\method{plot}{hsmm_simulation}(x, ...)
}
\arguments{
\item{x}{A fitted or simulated object.}
\item{object}{A fitted object.}
\item{type}{Return decoded states or posterior probabilities.}
\item{...}{Additional plotting arguments.}
}
\description{Basic print, summary, prediction and plotting methods.}
