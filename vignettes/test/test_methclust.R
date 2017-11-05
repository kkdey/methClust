
##################   Test  methClust   ######################

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

library(slam)
library(methClust)
topics <- meth_topics(M, U, K=2, tol = 10)

omega <- topics$omega

barplot(t(omega),
        col = 2:(K+1),
        axisnames = F, space = 0,
        border = NA,
        main=paste("No. of clusters=", K),
        las=1, ylim = c(0,1), cex.axis=1.5,cex.main=1.4)
