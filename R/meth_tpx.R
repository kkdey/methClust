
CheckCounts <- function(counts){
  if(class(counts)[1] == "TermDocumentMatrix"){ counts <- t(counts) }
  if(is.null(dimnames(counts)[[1]])){ dimnames(counts)[[1]] <- paste("doc",1:nrow(counts)) }
  if(is.null(dimnames(counts)[[2]])){ dimnames(counts)[[2]] <- paste("wrd",1:ncol(counts)) }
  empty <- row_sums(counts) == 0
  if(sum(empty) != 0){
    counts <- counts[!empty,]
    cat(paste("Removed", sum(empty), "blank documents.\n")) }
  return(as.simple_triplet_matrix(counts))
}


meth_tpxinit <- function(meth_X, unmeth_X, inifreq, K1, alpha, verb, nbundles=1,
                    use_squarem=FALSE, init.adapt){
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

    ## Solve for map omega in NEF space
    fit <- meth_tpxfit(meth_X=meth_X, unmeth_X = unmeth_X,
                       freq=inifreq, tol=100, verb=0,
                       use_squarem = FALSE, admix=TRUE, grp=NULL,
                       tmax = 10, wtol=10^{-2}, qn=100, MAXITER_SQUAREM = 3)

  return(fit$freq)
}

meth_tpxfit <- function(meth_X, unmeth_X, freq, tol, verb,
                        use_squarem, admix, grp, tmax, wtol,
                        qn, MAXITER_SQUAREM = 3){
  if(!inherits(meth_X,"simple_triplet_matrix")){ stop("meth_X needs to be a simple_triplet_matrix") }
  if(!inherits(unmeth_X,"simple_triplet_matrix")){ stop("unmeth_X needs to be a simple_triplet_matrix") }

  K <- ncol(freq)
  n <- nrow(meth_X)
  p <- ncol(meth_X)

  mvo <- meth_X$v[order(meth_X$i)]
  uvo <- unmeth_X$v[order(unmeth_X$i)]
  wrd <- meth_X$j[order(meth_X$i)]-1
  doc <- c(0,cumsum(as.double(table(factor(meth_X$i, levels=c(1:nrow(meth_X)))))))

  system.time(omega <- meth_tpxweights(n=n, p=p, mvo=mvo, uvo= uvo, wrd=wrd,
                          doc=doc,
                          start=meth_tpxOmegaStart(meth_X, unmeth_X, freq),
                          freq=freq))

  ## tracking
  iter <- 0
  dif <- tol+1+qn
  update <- TRUE
  if(verb>0){
    cat("log posterior increase: " )
    digits <- max(1, -floor(log(tol, base=10))) }

  Y <- NULL # only used for qn > 0
#  Q0 <- col_sums(X)/sum(X)
  L <- meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                     freq=freq, omega=omega,
                     alpha=alpha, admix=admix, grp=grp)

  iter <- 1;
  while( update  && iter < tmax ){
    if(admix && wtol > 0){
      Wfit <- tpxweights(n=nrow(X), p=ncol(X), mvo=mvo, uvo = uvo, wrd=wrd, doc=doc,
                         start=omega, freq=freq,  verb=0, nef=TRUE, wtol=wtol, tmax=20)
    }else{
      Wfit <- omega;
    }

    if(use_squarem){

      Wfit <- meth_normalizetpx(Wfit + 1e-15, byrow=TRUE);
      param_vec_in <- c(as.vector(logit(Wfit)),as.vector(logit(freq)));
      res <- squarem(par=as.numeric(param_vec_in),
                     fixptfn=tpxsquarEM,
                     objfn= tpxlpost_squarem,
                     meth_X=meth_X,
                     unmeth_X = unmeth_X,
                     K=K,
                     alpha=alpha,
                     admix=admix,
                     grp=grp,
                     control=list(maxiter = MAXITER_SQUAREM,
                                  trace = FALSE, square=TRUE, tol=1e-10));

      res_omega <- inv.logit(matrix(res$par[1:(nrow(X)*K)], nrow=nrow(X), ncol=K));
      res_freq <- inv.logit(matrix(res$par[-(1:(nrow(X)*K))], nrow=ncol(X), ncol=K));

      move <- list("omega"=res_omega, "freq"=res_freq);
      QNup <- list("omega"=move$omega, "freq"=move$freq, "L"=res$value.objfn, "Y"=NULL)
    }

    if(!use_squarem){
      ## joint parameter EM update
      Wfit <- meth_normalizetpx(Wfit + 1e-15, byrow=TRUE);
      move <- meth_tpxEM(meth=as.matrix(meth_X), unmeth = as.matrix(unmeth_X),
                    freq=freq, omega=Wfit, admix=admix, grp=grp)
      QNup <- move
      QNup$L <-  meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                               freq=move$freq, omega=move$omega,
                               admix=admix, grp=grp)
      # QNup <- meth_tpxQN(move=move, Y=Y, meth_X=meth_X, unmeth_X = unmeth_X,
      #               alpha=alpha, verb=verb,
      #               admix=admix, grp=grp, doqn=qn-dif)
      # Y <- QNup$Y
    }

    if(QNup$L < L){  # happens on bad Wfit, so fully reverse
      if(verb > 10){ cat("_reversing a step_") }
      move <- meth_tpxEM(meth_X=meth_X, unmeth_X = unmeth_X,
                         freq=freq, omega=omega, admix=admix, grp=grp)
      QNup$L <-  meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                               freq=move$freq, omega=move$omega,
                               admix=admix, grp=grp)
    }

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

    ## iterate
    iter <- iter+1
    freq <- move$freq
    omega <- move$omega
  }

  ## final log posterior
  L <- meth_tpxlpost(meth_X=meth_X, unmeth_X = unmeth_X,
                     freq=freq, omega=omega, alpha=alpha, admix=admix, grp=grp)

  ## summary print
  if(verb>0){
    cat("done.")
    if(verb>1) { cat(paste(" (L = ", round(L,digits), ")", sep="")) }
    cat("\n")
  }

  out <- list(freq=freq, omega=omega, K=K, L=L, iter=iter)
  invisible(out)
}


meth_tpxlpost <- function(meth_X, unmeth_X, freq, omega, admix=TRUE, grp=NULL)
{
  omega[omega==1] <- 1 - 1e-10;
  omega[omega==0] <- 1e-10;
  omega <- meth_normalizetpx(omega, byrow = TRUE)
  if(!inherits(meth_X,"simple_triplet_matrix")){ stop("meth_X needs to be a simple_triplet_matrix.") }
  if(!inherits(unmeth_X,"simple_triplet_matrix")){ stop("unmeth_X needs to be a simple_triplet_matrix.") }

  K <- ncol(freq)

  if(admix){
    prob <- meth_tpxQ(freq=freq, omega=omega, doc=meth_X$i, wrd=meth_X$j)
    L <- sum( meth_X$v*log(prob) +  unmeth_X$v*log(1 - prob))
  }else{ L <- sum(meth_tpxMixQ(X, omega, freq, grp)$lqlhd) }
  L <- L + sum(log(omega))/K

  return(L) }

meth_tpxllik <- function(meth_X, unmeth_X, freq, omega){
  wrd <- meth_X$j[order(meth_X$i)]-1
  prob <- omega %*% t(freq)
  prob_X <- slam::as.simple_triplet_matrix(prob)
  out <- .C("wllhd",
            nwrd = as.integer(length(wrd)),
            m = as.double(meth_X$v),
            u = as.double(unmeth_X$v),
            q = as.double(prob_X$v),
            package = "methClust")
  return(out)
}

meth_tpxlpost_squarem <- function(param_vec_in,  meth, unmeth,  K,
                             admix=TRUE, grp=NULL)
{
  meth_X <- CheckCounts(meth)
  unmeth_X <- CheckCounts(unmeth)
  omega_in <- inv.logit(matrix(param_vec_in[1:(nrow(meth)*K)], nrow=nrow(meth), ncol=K));
  freq_in <- inv.logit(matrix(param_vec_in[-(1:(nrow(meth)*K))], nrow=ncol(meth), ncol=K))
  return(meth_tpxlpost(meth_X, unmeth_X, freq_in, omega_in, admix, grp))
}


meth_tpxsquarEM <- function(param_vec_in, meth, unmeth, K, admix, grp){
  omega_in <- inv.logit(matrix(param_vec_in[1:(nrow(meth)*K)], nrow=nrow(meth), ncol=K));
  freq_in <- inv.logit(matrix(param_vec_in[-(1:(nrow(meth)*K))], nrow=ncol(meth), ncol=K))
  out <- meth_tpxEM(meth, unmeth, freq_in, omega_in, admix, grp);
  param_vec_out <- c(as.vector(logit(out$omega)),as.vector(logit(out$freq)))
  return(param_vec_out)
}

meth_tpxEM <- function(meth, unmeth, freq, omega, admix, grp)
{
  n <- nrow(omega)
  p <- nrow(freq)
  K <- ncol(freq)

  omega_in <- omega
  freq_in <- freq
  W <- array(0, c(dim(omega_in)[1], dim(freq_in)[1], dim(omega_in)[2]))

  for(n in 1:dim(omega_in)[1]){
    W[n, , ] <- t(replicate(dim(freq_in)[1], omega_in[n,])) * freq_in
  }

  V <- aperm(replicate(dim(freq_in)[1], omega_in), c(1,3,2)) - W
  Wnorm <- aperm(apply(W, c(1,2), function(x) return(x/sum(x))), c(2,3,1))
  Vnorm <- aperm(apply(V, c(1,2), function(x) return(x/sum(x))), c(2,3,1))

  A <- replicate(dim(omega_in)[2], meth) * Wnorm
  B <- replicate(dim(omega_in)[2], unmeth) * Vnorm

  omega_unscaled <- apply(A, c(1,3), function(x) return(sum(x))) + apply(B, c(1,3), function(x) return(sum(x)))
  omega_out <- t(apply(omega_unscaled, 1, function(x) return(x/sum(x))))

  freq_out <- 1/(1 + apply(B, c(2,3), sum)/apply(A, c(2,3), sum))

  return(list(freq=freq_out, omega=omega_out))
}


meth_tpxQ <- function(freq, omega, doc, wrd){

  omega[omega==1] <- 1 - 1e-14;
  omega[omega==0] <- 1e-14;
  omega <- meth_normalizetpx(omega, byrow = TRUE)

  if(length(wrd)!=length(doc)){stop("index mis-match in tpxQ") }
  if(ncol(omega)!=ncol(freq)){stop("freq/omega mis-match in tpxQ") }

  out <- .C("RcalcQ",
            n = as.integer(nrow(omega)),
            p = as.integer(nrow(freq)),
            K = as.integer(ncol(freq)),
            doc = as.integer(doc-1),
            wrd = as.integer(wrd-1),
            N = as.integer(length(wrd)),
            omega = as.double(omega),
            freq = as.double(freq),
            q = double(length(wrd)),
            PACKAGE="methClust" )

  return( out$q ) }


meth_tpxOmegaStart <- function(meth_X, unmeth_X, freq)
{
  if(!inherits(meth_X,"simple_triplet_matrix")){ stop("meth_X needs to be a simple_triplet_matrix.") }
  if(!inherits(unmeth_X,"simple_triplet_matrix")){ stop("unmeth_X needs to be a simple_triplet_matrix.") }
  prop_X <- (meth_X)/(meth_X + unmeth_X)
  omega <- try(slam::tcrossprod_simple_triplet_matrix(prop_X, solve(t(freq)%*%freq)%*%t(freq)), silent=TRUE )
  if(inherits(omega,"try-error")){ return( matrix( 1/ncol(freq), nrow=nrow(meth_X), ncol=ncol(freq) ) ) }
  omega[omega < 0] <- 0
  return( meth_normalizetpx(omega, byrow=TRUE) )
}

meth_tpxweights <- function(n, p, mvo, uvo, wrd, doc, start,
                            freq, verb=FALSE, nef=TRUE, wtol=10^{-5},
                            tmax=1000)
{
  K <- ncol(freq)
  start[start == 0] <- 0.1/K
  start <- start/rowSums(start)
  omega <- .C("Romega",
              n = as.integer(n),
              p = as.integer(p),
              K = as.integer(K),
              doc = as.integer(doc),
              wrd = as.integer(wrd),
              m = as.double(mvo),
              u = as.double(uvo),
              freq = as.double(freq),
              W = as.double(t(start)),
              nef = as.integer(nef),
              tol = as.double(wtol),
              tmax = as.integer(tmax),
              verb = as.integer(verb),
              PACKAGE="methClust")
  return(t(matrix(omega$W, nrow=ncol(freq), ncol=n)))
}


