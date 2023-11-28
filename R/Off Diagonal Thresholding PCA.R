#Off Diagonal Sparse PCA
Off_diagonal <- function (Sigma_hat,#Sample Covariance Matrix
                          n, #Sample Size
                          r, #Number of Spikes
                          alpha=2)#Threshold Parameter
{

  p <- nrow(Sigma_hat)

  Sigma_threshold <-  Sigma_hat * (abs(Sigma_hat)>alpha*sqrt(log(p)/n))

  N <- Sigma_threshold-diag(diag(Sigma_threshold))


  tryCatch(
    {

      N_eig <-  irlba(N,nu=r)

    },
    error = function(e){
      N_eig <- list(u = matrix(rnorm(p*r),ncol = r),
                    d = seq(0,r));

      N_eig$u <- qr.Q(qr(N_eig$u))
    }
  )
  #print(N_eig)
  return(list(eigvector=N_eig$u, eigvalue=N_eig$d[1:r]))

}
