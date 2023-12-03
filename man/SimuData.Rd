\name{Generate Samples of Generalized Spiked Covariance Model}
\alias{SimuData}
\title{Sparse Heteroskedastic PCA in High Dimensions}
\usage{
SimuData(n,mu,U,lambda,sigma)
}
\arguments{
\item{n}{The sample size.}
\item{mu}{The length-\eqn{p} mean vector.}
\item{U}{A \eqn{p \times r} orthogonal matrix serving as the loading matrix.}
\item{lambda}{The spike strength, should be a length-\eqn{r} vector.}
\item{sigma}{The standard deviation of the noises, which should be a \eqn{p \times p} diagonal matrix.}
}
\value{
\item{eigvector}{The initial estimation of the loading matrix, which is a \eqn{p \times r} orthogonal matrix.}
\item{eigvalue}{The initial estimation of the spike strength, which is a length-\eqn{r} vector.}
}

\description{
This function generates samples under the generalized spiked covariance model.
}
\examples{

r <- 1
n <- 256
p <- 512
mu <- rep(0,p)
s <- floor(p*0.1)
lambda <- 2

set.seed(7)
v <- runif(p)
alpha <- 5
sigma2 <- diag(0.1*p*v^alpha/sum(v^alpha))
U <- matrix(0,nrow = p, ncol = r)
U[1:s,] <- 1/sqrt(s)

X <- SimuData(n,mu,U,lambda,sqrt(sigma2))
}
