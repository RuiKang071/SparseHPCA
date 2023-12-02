SimuData <- function(n,mu,U,lambda,sigma){
  #This function generate the required data
  #Inputs
  #mu       1 * p vector, mean vector
  #U      p * m matrix, m principal components
  #lambda   1 * m vector, spikes
  #sigma	  p * p diagonal matrix, noise sd
  #Output
  #X	      data matrix,n * p matrix, n observations


  p = nrow(U)
  m = ncol(U)#m is the rank

  Z = matrix(rnorm(n*p,0,1),ncol = p,nrow = n)
  V = matrix(rnorm(n*m,0,1),ncol = m,nrow = n)

  #-------------------Generate the data--------------------
  #lambda <- rep(sqrt(lambda),times = m)
  #lambda <- matrix(lambda,nrow = m,ncol = m,byrow = FALSE)
  lambda <- diag(sqrt(lambda),nrow = m)
  mu <- rep(mu,times = n)
  mu <- matrix(mu,nrow = n,ncol = p,byrow = TRUE)

  return(X = mu +  V %*%lambda%*%t(U) + Z%*%sigma)
}
