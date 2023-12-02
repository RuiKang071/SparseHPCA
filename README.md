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
- Obtain the initial estimator.
- Use SparseHPCA to obtain the final result.
