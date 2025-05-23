% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ex_expert_sensitivity}
\alias{ex_expert_sensitivity}
\title{Expert-based sensitivity and adaptive capacity scores for four indicators and
five pressures}
\format{
A data frame with 40 observations and 23 variables.
\describe{
\item{indicator}{Name of assessed indicator.}
\item{pressure}{Name of assessed pressure.}
\item{type}{Effect type (direct, indirect, or direct + indirect).}
\item{sens_feeding}{Sensitivity score for the feeding trait (-5 to 5).}
\item{sens_behaviour}{Sensitivity score for the behaviour trait (-5 to 5).}
\item{sens_reproduction}{Sensitivity score for the reproduction trait (-5 to 5).}
\item{sens_habitat}{Sensitivity score for the habitat trait (-5 to 5)}
\item{sens_general}{General sensitivity score (where trait-based scoring is not possible, -5 to 5).}
\item{ac_feeding}{Adaptive capacity score for the feeding trait (-1 to 1).}
\item{ac_behaviour}{Adaptive capacity score for the behaviour trait (-1 to 1).}
\item{ac_reproduction}{Adaptive capacity score for the reproduction trait (-1 to 1).}
\item{ac_habitat}{Adaptive capacity score for the habitat trait (-1 to 1).}
\item{ac_general}{General adaptive capacity score (where trait-based scoring is not possible, -1 to 1).}
\item{uncertainty_sens_feeding}{Uncertainty of sensitivity score for feeding trait (1 to 3).}
\item{uncertainty_sens_behaviour}{Uncertainty of sensitivity score for behaviour trait (1 to 3).}
\item{uncertainty_sens_reproduction}{Uncertainty of sensitivity score for reproduction trait (1 to 3).}
\item{uncertainty_sens_habitat}{Uncertainty of sensitivity score for habitat trait (1 to 3).}
\item{uncertainty_sens_general}{Uncertainty of general sensitivity score (1 to 3).}
\item{uncertainty_ac_feeding}{Uncertainty of adaptive capacity score for feeding trait (1 to 3).}
\item{uncertainty_ac_behaviour}{Uncertainty of adaptive capacity score for behaviour trait (1 to 3).}
\item{uncertainty_ac_reproduction}{Uncertainty of adaptive capacity score for reproduction trait (1 to 3).}
\item{uncertainty_ac_habitat}{Uncertainty of adaptive capacity score for habitat trait (1 to 3).}
\item{uncertainty_ac_general}{Uncertainty of general adaptive capacity score (1 to 3).}
}
}
\usage{
ex_expert_sensitivity
}
\description{
This demo dataset includes sensitivity and adaptive capacity scores for four
Baltic Sea indicators and five pressures, initialized using
\code{\link{create_template_sensitivity}}. Depending on the indicator, a
general score or trait-specific scores were assigned. This dataset serves as
input for \code{\link{calc_sensitivity}}.
}
\keyword{datasets}
