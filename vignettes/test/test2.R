
#################  Test for loglik similarity   #####################

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
freq_sim <- cbind(c(rep(0.8, m.out), rep(0.2, m.out), rep(0.5, m.out), rep(0.01, m.out)),
              c(rep(0.01, m.out), rep(0.01, m.out), rep(0.5, m.out), rep(0.8, m.out)))

prob <- omega_sim %*% t(freq_sim)

Y <- matrix(rpois(dim(prob)[1]*dim(prob)[2], 1000), dim(prob)[1], dim(prob)[2])

M <- matrix(0, dim(Y)[1], dim(Y)[2])

for(m in 1:dim(Y)[1]){
  for(n in 1:dim(Y)[2]){
    M[m,n] <- rbinom(1, Y[m,n], prob = prob[m,n])
  }
}

U = Y - M

methloglik <- function(M, U, omega, freq)
{
  prob <- omega %*% t(freq)
  return(sum(M*log(prob) + U*log(1-prob)))
}


methloglik(M, U, omega_sim, freq_sim)

M[M == 0] <- 1e-10
meth_X <- slam::as.simple_triplet_matrix(M)
unmeth_X <- slam::as.simple_triplet_matrix(U)
prob_X <- slam::as.simple_triplet_matrix(prob)
wrd <- meth_X$j[order(meth_X$i)]-1

meth_tpxlpost <- function(meth_X, unmeth_X, freq, omega, admix=TRUE, grp=NULL)
{
  omega[omega==1] <- 1 - 1e-10;
  omega[omega==0] <- 1e-10;
  omega <- meth_normalizetpx(omega, byrow = TRUE)
  if(!inherits(meth_X,"simple_triplet_matrix")){ stop("meth_X needs to be a simple_triplet_matrix.") }
  if(!inherits(unmeth_X,"simple_triplet_matrix")){ stop("unmeth_X needs to be a simple_triplet_matrix.") }
  K <- ncol(freq)

  if(admix){
    prob <- meth_tpxQ(freq=freq, omega=omega, doc=X$i, wrd=X$j)
    L <- sum( meth_X$v*log(prob) +  unmeth_X$v*log(1 - prob))
  }else{ L <- sum(meth_tpxMixQ(X, omega, freq, grp)$lqlhd) }
  L <- L + sum(log(omega))/K

  return(L) }

meth_tpxlpost(meth_X, unmeth_X, freq_sim, omega_sim)








