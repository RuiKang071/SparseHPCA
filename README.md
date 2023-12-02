# SparseHPCA

Sparse Heteroskedastic PCA

# Description

This package implements a novel algorithm to estimate the loading matrix under the generalized covariance model.
We assume the data $X_i \in \mathbb{R}^p$ are i.i.d. generated as follows:
$$X_i = \mu + U \Lambda^{1/2} Z_i + \epsilon_i,$$
where $U \in \mathbb{O}^{p \times r}$, $\Lambda = \textnormal{diag}(\lambda_1,\dots,\lambda_r)$, $Z_i \sim \mathcal{N}(0,I_r)$ and $\epsilon_i \sim \mathcal{N}(0, \textnormal{diag}(\sigma_1^2,\dots,\sigma_p^2))$ are independent from each other.
The loading matrix $U$ is assumed to be $s$-sparse, i.e., at most $s$ rows are non-zero.

# Installation

- Step 1: Install the devtools package

```
install.packages("devtools")
```

- Step 2: Install the package from this repository

```
library(devtools)
install_github("RuiKang071/SparseHPCA")
```

# Functions
There are three main functions in this package.
- `Generate`
- `Adaptive off-diagonal thresholding algorithm` Our algorithm requires an initial estimator that is not perpendicular to the population truth. 
We provide the following function to get the initial estimator.
- `SparseHPCA`

# Examples
- Generate the random samples.
  ```
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
  ```
- Obtain the initial estimator.
 ```
  Sigma_hat <- cov(X)

  initial <- Adap_off_diag(X,r,alpha = 2)
  U0 <- initial$eigvector
  lambda_hat <- initial$eigvalue
  sin_theta_dist(U0,U)
  ```
- Use SparseHPCA to obtain the final result.
  ```
  U_hat <- SparseHPCA(Sigma_hat,n,r,U0,lambda_hat,sqrt(diag(sigma2)),
                             T=1000,gamma=0.8)
  sin_theta_dist(U_hat,U)
  ```
