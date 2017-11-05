
########################  test  tpxweights   ##############################

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

M[M == 0] <- 1e-10
meth_X <- slam::as.simple_triplet_matrix(M)
#meth_X$v[meth_X$v == 1e-10] = 0
unmeth_X <- slam::as.simple_triplet_matrix(U)

if(!inherits(meth_X,"simple_triplet_matrix")){ stop("meth_X needs to be a simple_triplet_matrix") }
if(!inherits(unmeth_X,"simple_triplet_matrix")){ stop("unmeth_X needs to be a simple_triplet_matrix") }

K <- ncol(freq_sim)
n <- nrow(meth_X)
p <- ncol(meth_X)

mvo <- meth_X$v[order(meth_X$i)]
uvo <- unmeth_X$v[order(unmeth_X$i)]
wrd <- meth_X$j[order(meth_X$i)]-1
doc <- c(0,cumsum(as.double(table(factor(meth_X$i, levels=c(1:nrow(meth_X)))))))

system.time(omega <- meth_tpxweights(n=n, p=p, mvo=mvo, uvo= uvo, wrd=wrd,
                                     doc=doc,
                                     start=meth_tpxOmegaStart(meth_X, unmeth_X, freq_sim),
                                     freq=freq_sim))

system.time(omega3 <- meth_tpxweights(n=n, p=p, mvo=mvo, uvo= uvo, wrd=wrd,
                                     doc=doc,
                                     start=omega_sim,
                                     freq=freq_sim))

omega2 <- meth_tpxOmegaStart(meth_X, unmeth_X, freq_sim)


prop_X <- (meth_X)/(meth_X + unmeth_X)
omega4 <- try(slam::tcrossprod_simple_triplet_matrix(prop_X, solve(t(freq_sim)%*%freq_sim)%*%t(freq_sim)), silent=TRUE )



## tracking
iter <- 0
dif <- tol+1+qn
update <- TRUE
verb <- 1
if(verb>0){
  cat("log posterior increase: " )
  digits <- max(1, -floor(log(tol, base=10))) }

Y <- NULL # only used for qn > 0
Q0 <- col_sums(X)/sum(X)
L <- meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                   freq=freq_sim, omega=omega_sim,
                   admix=TRUE, grp=grp)

wtol <- 10^{-4}
admix <- TRUE
if(admix && wtol > 0){
  Wfit <- meth_tpxweights(n=nrow(meth_X), p=ncol(meth_X), mvo=mvo, uvo = uvo,
                     wrd=wrd, doc=doc,
                     start=omega_sim, freq=freq_sim,
                     verb=0, nef=TRUE, wtol=wtol, tmax=20)
}else{
  Wfit <- omega;
}

Wfit <- meth_normalizetpx(Wfit + 1e-15, byrow=TRUE);
param_vec_in <- c(as.vector(logit(Wfit)),as.vector(logit(freq_sim)));
MAXITER_SQUAREM <- 3
res <- SQUAREM::squarem(par=as.numeric(param_vec_in),
                        fixptfn=meth_tpxsquarEM,
                        objfn= meth_tpxlpost_squarem,
                        meth = as.matrix(meth_X),
                        unmeth = as.matrix(unmeth_X),
                        K=K,
                        admix=admix,
                        grp=grp,
                        control=list(maxiter = MAXITER_SQUAREM,
                            trace = FALSE, square=TRUE, tol=1e-5));

res_omega <- inv.logit(matrix(res$par[1:(nrow(meth_X)*K)], nrow=nrow(meth_X), ncol=K));
res_freq <- inv.logit(matrix(res$par[-(1:(nrow(meth_X)*K))], nrow=ncol(meth_X), ncol=K));

move <- list("omega"=res_omega, "freq"=res_freq);

Wfit <- meth_normalizetpx(Wfit + 1e-15, byrow=TRUE);
move <- meth_tpxEM(meth=as.matrix(meth_X), unmeth = as.matrix(unmeth_X),
                   freq=freq_sim, omega=Wfit, admix=admix, grp=grp)
QNup <- move
QNup$L <-  meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                         freq=move$freq, omega=move$omega,
                         admix=admix, grp=grp)

move <- meth_tpxEM(meth= as.matrix(meth_X), unmeth = as.matrix(unmeth_X),
                   freq=freq_sim, omega=omega, admix=admix,
                   grp=grp)
QNup$L <-  meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                         freq=move$freq, omega=move$omega,
                         admix=admix, grp=grp)

dif <- (QNup$L-L)

L <- QNup$L


## check convergence
if(abs(dif) < tol){
  if(sum(abs(freq-move$freq)) < tol){ update = FALSE } }

## print
if(verb>0 && (iter-1)%%ceiling(10/verb)==0 && iter>0){
  ##if(verb>0 && iter>0){
  cat( paste( round(dif,digits), #" (", sum(abs(freq-move$freq)),")",
              ", ", sep="") ) }

## heartbeat for long jobs
if(((iter+1)%%1000)==0){
  cat(sprintf("p %d iter %d diff %g\n",
              nrow(freq), iter+1,round(dif))) }

