
############  test the meth_topics initialization  ######################

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


meth <- M
unmeth <- U

shape=NULL
initopics=NULL
tol=0.1
ord=TRUE
verb=1
use_squarem=TRUE
sample_init = TRUE


inifreq <- initopics
if(is.matrix(inifreq)){
  if(ncol(inifreq) != K){stop("mismatch between inifreq no. of columns and K")}
  if(max(inifreq) > 1 | min(inifreq) < 0){stop("some value of inifreq lies outside [0,1] range")}
}

if(is.null(inifreq)){ ilength <- K-1 }else{ ilength <- inifreq[1] }
if(ilength < 1){ ilength <- 1 }

nK <- length( Kseq <-  unique(ceiling(seq(2,K,length=ilength))) )

prop_methyl <- meth_X/(meth_X + unmeth_X)


freq_start <- matrix(rep(0.5, ncol(meth_X)))
omega_start <- matrix(rep(1, nrow(meth_X)))
Kpast <- ncol(freq_start)
Kdiff <- K-Kpast
n <- nrow(meth_X)
ki <- matrix(1:(n-n%%Kdiff), ncol=Kdiff)
for(i in 1:Kdiff){
  freq_start <- cbind(freq_start, col_means(prop_methyl[ki[,i],]))
}

inifreq <- freq_start




meth_X <- meth_X[index_init,]
unmeth_X <- unmeth_X[index_init,]
inifreq <- NULL
