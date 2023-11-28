# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

r <- 1
n <- 256
p <- 512
mu <- rep(0,p)
s <- floor(p*0.1)
lambda <- 2

gamma_para <- 0.8
alpha <- 5

set.seed(7)
v <- runif(p)
sigma2 <- diag(0.1*p*v^alpha/sum(v^alpha))
U <- matrix(0,nrow = p, ncol = r)
U[1:s,] <- 1/sqrt(s)

Z <-  matrix(rnorm(n*p,0,1),ncol = p,nrow = n)
V <-  matrix(rnorm(n*r,0,1),ncol = r,nrow = n)
mu <- rep(mu,times = n)
mu <- matrix(mu,nrow = n,ncol = p,byrow = TRUE)
X <-  mu +  V %*%diag(sqrt(lambda),nrow = r)%*%t(U) + Z%*%sqrt(sigma2)


Sigma_hat <- cov(X)

temp <- Adap_off_diag(X,r,alpha = 2)
U0 <- temp$eigvector
lambda_hat <- temp$eigvalue
sin_theta_dist(U0,U)

U_hat <- SparseHPCA(Sigma_hat,n,r,U0,lambda_hat,sqrt(diag(sigma2)),
                           T=1000,gamma=gamma_para)
sin_theta_dist(U_hat,U)
