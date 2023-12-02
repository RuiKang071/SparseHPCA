\name{Adaptive Off-diagonal Thresholding Sparse PCA}
\alias{Adap_off_diag}
\title{Sparse Heteroskedastic PCA in High Dimensions}
\usage{
Adap_off_diag <- function(X, r, alpha = 2)
}
\arguments{
\item{X}{A \eqn{n \times p} sample matrix where \eqn{n} is the sample size and \eqn{p} is the dimension. Dimension \eqn{p} is allowed to be greater than \eqn{n}.}
\item{r}{The rank of the loading matrix.}
\item{alpha}{The thresholding parameter.}
}
\value{
\item{eigvector}{The initial estimation of the loading matrix, which is a \eqn{p \times r} orthogonal matrix.}
\item{eigvalue}{The initial estimation of the spike strength, which is a length-\eqn{r} vector.}
}

\description{
This function provides an initial estimator for the main function SparseHPCA.
}
\examples{
#Generate the data

r <- 1
n <- 256
p <- 512
mu <- rep(0,p)
s <- floor(p*0.1)
lambda <- 2

set.seed(7)
v <- runif(p)
alpha <- 5#Control the degree of heteroskedasity
sigma2 <- diag(0.1*p*v^alpha/sum(v^alpha))
U <- matrix(0,nrow = p, ncol = r)
U[1:s,] <- 1/sqrt(s)

X <- SimuData(n,mu,U,lambda,sqrt(sigma2))

Sigma_hat <- cov(X)

initial <- Adap_off_diag(X,r,alpha = 2)
U0 <- initial$eigvector
lambda_hat <- initial$eigvalue
}
