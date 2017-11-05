
###########  Simulation model for methclust   ########################

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

for(num in 1:NUM_ITER){
  W <- array(0, c(dim(omega_iter)[1], dim(freq_iter)[1], dim(omega_iter)[2]))

  system.time(for(n in 1:dim(omega_iter)[1]){
      W[n, , ] <- t(replicate(dim(freq_iter)[1], omega_iter[n,])) * freq_iter
  })


  V <- aperm(replicate(dim(freq_iter)[1], omega_iter), c(1,3,2)) - W
  system.time(
  Wnorm <- aperm(apply(W, c(1,2), function(x) return(x/sum(x))), c(2,3,1)))
  system.time(Vnorm <- aperm(apply(V, c(1,2), function(x) return(x/sum(x))), c(2,3,1))
  )

  A <- replicate(dim(omega_iter)[2], M) * Wnorm
  B <- replicate(dim(omega_iter)[2], U) * Vnorm

  omega_unscaled <- apply(A, c(1,3), function(x) return(sum(x))) + apply(B, c(1,3), function(x) return(sum(x)))
  omega_iter <- t(apply(omega_unscaled, 1, function(x) return(x/sum(x))))

  system.time(freq_iter <- 1/(1 + apply(B, c(2,3), sum)/apply(A, c(2,3), sum)))

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




