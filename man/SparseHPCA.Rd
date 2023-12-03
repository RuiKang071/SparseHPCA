\name{Sparse Heteroskedastic PCA}
\alias{SparseHPCA}
\title{Sparse Heteroskedastic PCA in High Dimensions}
\usage{
SparseHPCA(Sigma_hat,n,r,U0,lambda_hat,noise_sd,T=100,gamma=1)
}
\arguments{\item{Sigma_hat}{A \eqn{p \times p} symmetric matrix of which the loading matrix is to be estimated.}
\item{n}{The sample size.}
\item{r}{The rank of the loading matrix.}
\item{U0}{A \eqn{p \times r} orthogonal matrix serving as the initial estimator.}
\item{lambda_hat}{The estimated spike strength which should be a length-\code{r} vector.}
\item{noise_sd}{The standard deviations of the noise which should be a length-\eqn{p} vector.}
\item{T}{The maximum number of iterations.}
\item{gamma}{The thresholding parameter.}}
\value{
\item{U}{The final estimation of the loading matrix, which is a \eqn{p \times r} orthogonal matrix.}
}
\description{
The main function for sparse heteroskedastic PCA.
}
\examples{
r <- 1
n <- 256
p <- 512
mu <- rep(0,p)
s <- floor(p*0.1)
lambda <- 2

set.seed(7)
alpha <- 5
v <- runif(p)
sigma2 <- diag(0.1*p*v^alpha/sum(v^alpha))
U <- matrix(0,nrow = p, ncol = r)
U[1:s,] <- 1/sqrt(s)

X <- SimuData(n,mu,U,lambda,sqrt(sigma2))

Sigma_hat <- cov(X)

initial <- Adap_off_diag(X,r,alpha = 2)
U0 <- initial$eigvector
lambda_hat <- initial$eigvalue
sin_theta_dist(U0,U)

U_hat <- SparseHPCA(Sigma_hat,n,r,U0,lambda_hat,sqrt(diag(sigma2)),
                           T=1000,gamma=0.8)
sin_theta_dist(U_hat,U)
}
