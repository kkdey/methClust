n.out <- 500
omega_sim <- rbind( cbind( rep(1, n.out), rep(0, n.out)),
                    cbind( rep(0, n.out), rep(1, n.out)),
                    cbind( seq(0.6, 0.4, length.out = n.out),
                           1- seq(0.6, 0.4,length.out=n.out)) )
dim(omega_sim)

K <- dim(omega_sim)[2]
barplot(t(omega_sim),
        col = 2:(K+1),
        axisnames = F, space = 0,
        border = NA,
        main=paste("No. of clusters=", K),
        las=1, ylim = c(0,1), cex.axis=1.5,cex.main=1.4)

m.out <- 200
freq <- cbind(c(rep(0.8, m.out), rep(0.2, m.out), rep(0.5, m.out), rep(0.01, m.out)),
              c(rep(0.01, m.out), rep(0.01, m.out), rep(0.5, m.out), rep(0.8, m.out)))

prob <- omega_sim %*% t(freq)

Y <- matrix(rpois(dim(prob)[1]*dim(prob)[2], 1000), dim(prob)[1], dim(prob)[2])

M <- matrix(0, dim(Y)[1], dim(Y)[2])

for(m in 1:dim(Y)[1]){
  for(n in 1:dim(Y)[2]){
    M[m,n] <- rbinom(1, Y[m,n], prob = prob[m,n])
  }
}

plot(freq[,2])

omega_init <- gtools::rdirichlet(dim(omega_sim)[1], c(10,10))
freq_init <- matrix(0.5, dim(freq)[1], dim(freq)[2])

omega_iter <- omega_init
freq_iter<- freq_init

NUM_ITER <- 100
U <- Y - M

methloglik <- function(M, U, omega, freq)
{
  prob <- omega %*% t(freq)
  return(sum(M*log(prob) + U*log(1-prob)))
}


loglik <- methloglik(M, U, omega_iter, freq_iter)

NUM_ITER <- 100

for(num in 1:NUM_ITER){

  m_lambda <- omega_iter %*% t(freq_iter)
  m_temp <- (M/m_lambda)
  m_matrix <- (m_temp %*% freq_iter)*omega_iter

  u_lambda <- omega_iter %*% t(1 - freq_iter)
  u_temp <- (U/u_lambda)
  u_matrix <- (u_temp %*% (1 - freq_iter)) * omega_iter

  omega_iter <- meth_normalizetpx(m_matrix + u_matrix + (1/(n*K)), byrow=TRUE)

  m_t_matrix <- (t(m_temp) %*% omega_iter)*freq_iter
  u_t_matrix <- (t(u_temp) %*% omega_iter)*(1-freq_iter)

  freq_iter <- m_t_matrix/(m_t_matrix + u_t_matrix)

  barplot(t(omega_iter),
          col = 2:(K+1),
          axisnames = F, space = 0,
          border = NA,
          main=paste("No. of clusters=", K),
          las=1, ylim = c(0,1), cex.axis=1.5,cex.main=1.4)

  cat("We are at iteration ", num, "\n")
  loglik_new <- methloglik(M, U, omega_iter, freq_iter)
  diff <- loglik_new - loglik
  loglik <- loglik_new

  if(num %% 5 == 0){
    cat("The posterior increase", round(diff,2), "\n")
  }

}
