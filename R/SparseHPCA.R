#Hetero Sparse PCA
SparseHPCA <- function(Sigma_hat,n,r,U0,lambda_hat,noise_sd,
                             T=100,gamma=1){
  # Power method
  #Input
  #Sigma_hat      Sample covariance matrix
  #n              sample size
  #r              Target subspace dimension
  #U0             initial estimate
  #lambda_hat     r by 1 vector, estimation of eigenvalues
  #noise_sd       estimation of sigma_i, i = 1,...,p
  #gamma          threshold parameter
  #Output
  #U              principal subspace estimator

  sigma_max = max(noise_sd)
  p = nrow(U0)
  lambda_hat <- matrix(rep(lambda_hat, p),nrow = p,byrow = TRUE)
  threshold = gamma*diag(noise_sd)%*%(sigma_max + sqrt(lambda_hat))*sqrt(log(p)/n)

  N = Sigma_hat - diag(diag(Sigma_hat))

  for (t in 1:T)
  {
    U=N%*%U0
    U=U*(abs(U)>threshold)
    U=qr.Q(qr(U))
    Lambda=t(U)%*%N%*%U
    Nnew=U%*%Lambda%*%t(U)

    if (sin_theta_dist(U,U0)<1e-5)
    {
      break
    }
    else
    {
      diag(N)=diag(Nnew)
      U0=U
    }
  }
  return (U)
}
