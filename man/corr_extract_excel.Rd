% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr_extract.R
\name{corr_extract_excel}
\alias{corr_extract_excel}
\title{Creates excel output from \link[mdrcAnalysis]{corr_extract_singleton}}
\usage{
corr_extract_excel(
  .dataset,
  .by_val,
  .vars,
  .by_x,
  .alpha,
  .prob,
  .wb,
  .plot,
  .correlation_style,
  .alpha_header_style
)
}
\arguments{
\item{.dataset}{Input tibble}

\item{.vars}{Optional, character vector of variables to correlate. Defaults to all
numeric variables not specified in other args (e.g. \code{.with_vars}, \code{.by_x})}

\item{.by_x}{Optional, string of variable name to separate correlations for}

\item{.alpha}{Boolean determining if Cronbach alphas will be calculated}

\item{.prob}{Boolean determining if t-test probability will be displayed. Calculated
using an alpha of 0.05.}

\item{.wb}{Workbook to write to}

\item{.plot}{Boolean determining if correlation plot(s) included in Excel}

\item{.correlation_style}{Default styling for correlation table}

\item{.alpha_header_style}{Styling for alphas grouped header}

\item{.datset}{Output from \link[mdrcAnalysis]{corr_extract_singleton}}

\item{.byval}{Specific value of .byx from \link[mdrcAnalysis]{corr_extract}}
}
\value{
NULl
}
\description{
Creates excel output from \link[mdrcAnalysis]{corr_extract_singleton}
}
\keyword{internal}
