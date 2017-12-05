
####################   test  methclust   ###################

meth <- t(read.table("../../../methClust-pages/data/meth.txt", header = TRUE))
unmeth <- t(read.table("../../../methClust-pages/data/unmeth_w_zero.txt", header = TRUE))

topics <- meth_topics(meth, unmeth, K=10, tol = 10, use_squarem = FALSE)

K = 8
shape=NULL
initopics=NULL
tol=0.1
ord=TRUE
verb=1
sample_init = TRUE
use_squarem=FALSE
NUM_INDICES_START = NULL


methx <- meth_X[index_init,]
unmethx <- unmeth_X[index_init,]

meth_X <- methx
unmeth_X <- unmethx
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

