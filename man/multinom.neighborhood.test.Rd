% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial.R
\name{multinom.neighborhood.test}
\alias{multinom.neighborhood.test}
\title{Perform the neighborhood test for multinom.test}
\usage{
multinom.neighborhood.test(x, y = NULL, delta = NULL)
}
\arguments{
\item{x, y}{Integer vectors (or matrices or dataframes containing multiple
integer vector observations as rows). \code{x} and \code{y} must be the
same type and dimension. If \code{x} and \code{y} are matrices (or
dataframes), the \eqn{i^th} row of \code{x} will be tested against the
\eqn{i^th} row of \code{y} for all \eqn{i} in 1..\code{nrow(x)}.
Alternatively, \code{x} can be a list of two vectors, matrices, or
dataframes to be compared. In this case, \code{y} is NULL by default.}

\item{delta}{A number (or vector) greater than 0.}
}
\value{
The \code{statistic} from \code{multinom.test} and its
associated \code{p_delta}, where \code{p_delta} \eqn{= 1 - pnorm(T - delta)}.
If \code{x} and \code{y} are two dimensional (that is, they are matrices
or dataframes with more than one row) and/or \code{delta} is a vector,
then a matrix will be returned where the \eqn{(i,j)^{th}} entry will be the
\code{p.delta} associated with the \eqn{i^{th}} rows of \code{x} and
\code{y} and the \eqn{j^{th}} entry of the \code{delta} vector.
}
\description{
Peforms the two sample test for two multinomial vectors
testing \eqn{H_0:} the underlying multinomial probability vectors
are within some neighborhood of one another vs. \eqn{H_1:} they are
not.
}
\details{
In testing the equality of parameters from two populations
(as in \code{multinom.test}),
it frequenly happens that the null hypothesis is rejected even though the estimates
of effect sizes are close to each other. However, these differences may be so small
that the parameters are not considered different in practice. A neighborhood test
is useful in this situation.
}
\examples{

# Load the twoNewsGroups dataset

data(twoNewsGroups)

# Sample two sets of 200 documents from the sci.med newsGroup (to simulate
# the null hypothesis being TRUE). For each of the two groups, sum the
# 200 term frequency vectors together. They will be the two vectors that
# we test.

num_docs <- 200
vecs2test <- list(NA, 2)
row_ids <- 1:nrow(twoNewsGroups$sci.med)
group_1 <- sample(row_ids, num_docs)
group_2 <- sample(row_ids[-group_1], num_docs)

vecs2test[[1]] <- twoNewsGroups$sci.med[group_1,] |>
                    colSums() |>
                    matrix(nrow=1)
vecs2test[[2]] <- twoNewsGroups$sci.med[group_2,] |>
                    colSums() |>
                    matrix(nrow=1)

# Test the null that the two vectors come from the same distribution
# (i.e. the same news group)

vecs2test |> multinom.test()

# The above test likely produced a significant p-value meaning that we would
# reject the null. However, the difference isn't very interesting. Instead,
# test that the differences are within some neighborhood:

vecs2test |> multinom.neighborhood.test(delta=60)


}
\seealso{
\code{\link{multinom.test}}, \code{vignette("multinomial-neighborhood-test-vignette")}

Amanda Plunkett & Junyong Park (2018), \emph{Two-Sample Test for Sparse High
Dimensional Multinomial Distributions}, TEST,
\url{https://doi.org/10.1007/s11749-018-0600-8}
}
