#Adaptive off diagonal threshold
Adap_off_diag <- function(X, #n by p data matrix
                          r, #number of spikes
                          alpha = 2)# Thresholding level
{
  f <- function(i,j,
                X,
                Sigma_hat){
    return(mean((X[,i]*X[,j] - Sigma_hat[i,j])^2))
  }
  Sigma_hat <- cov(X)
  p <- ncol(X)
  n <- nrow(X)
  x <- seq(1,p)
  y <- seq(1,p)
  Theta_hat <- matrix(
    mapply(f,rep(x,p),rep(y,each = p),
           MoreArgs = list(X = X,Sigma_hat = Sigma_hat)),
    ncol = p,nrow = p,byrow = FALSE)
  I_signal = (abs(Sigma_hat) > alpha * sqrt(Theta_hat*log(p)/n))


  Sigma_hat <- Sigma_hat * I_signal

  N=Sigma_hat-diag(diag(Sigma_hat))
  N_eig = irlba(N,nu=r)
  return(list(eigvector=N_eig$u, eigvalue=N_eig$d[1:r]))
}
