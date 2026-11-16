\name{hte_report}
\alias{hte_report}
\alias{print.hte_report}
\alias{format.hte_report}
\title{Plain-language summary of subgroup diagnostics}
\description{
Turns the output of \code{\link{compare_hte}} into short English paragraphs,
one per group, that state the group's size and share, its treated and control
counts, its difference-in-means estimate with standard error and interval, the
minimum detectable effect from a t-approximation power calculation, and
the range of estimates across estimators. Between-group contrasts are reported
in a separate table. The text reports only quantities computed from the data
and does not recommend an estimator.
}
\usage{
hte_report(x, groups = NULL, digits = 2, file = NULL)

\method{print}{hte_report}(x, width = 0.9 * getOption("width"), ...)

\method{format}{hte_report}(x, ...)
}
\arguments{
\item{x}{An object from \code{\link{compare_hte}} (or an \code{hte_report} for the methods).}
\item{groups}{Groups to describe; defaults to all, smallest first.}
\item{digits}{Digits for reported numbers.}
\item{file}{Optional path; if given, the paragraphs are also written to this file.}
\item{width}{Line width for printing.}
\item{...}{Ignored.}
}
\value{
An object of class \code{hte_report}: a list with \code{text} (one paragraph
per group), \code{contrasts} (a data frame of contrasts by estimator, with the
subgroup-means interval, range, and sign disagreement), and \code{note}.
\code{format()} returns the whole report as character lines.
}
\examples{
set.seed(1)
n <- c(White = 600, Black = 200, Asian = 40)
race <- rep(names(n), n)
treat <- unlist(lapply(n, function(m) sample(rep(0:1, length.out = m))))
y <- treat * c(White = 0.2, Black = 0.5, Asian = -0.2)[race] + rnorm(sum(n))
res <- compare_hte(data.frame(y, treat, race), "y", "treat", "race", methods = c("means", "ols"))
hte_report(res, groups = "Asian")
}
