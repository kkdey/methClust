

meth_topics <- function(meth,
                        unmeth,
                        K,
                        shape=NULL,
                        initopics=NULL,
                        tol=0.1,
                        ord=TRUE,
                        verb=1,
                        sample_init = TRUE,
                        NUM_INDICES_START = NULL,
                        use_squarem=FALSE){

  if(dim(meth)[1] != dim(unmeth)[1] | dim(meth)[2] != dim(unmeth)[2]){
    stop("meth and unmeth count matrices must have same dimensions")
  }

  ## sparsify the meth and unmeth count matrices

  meth[meth == 0] <- 1e-20
  unmeth[unmeth == 0] <- 1e-20

  meth_X <- CheckCounts(meth)
  unmeth_X <- CheckCounts(unmeth)
  meth_X$v[meth_X$v == 1e-20] = 0
  unmeth_X$v[unmeth_X$v == 1e-20] = 0

  zero_idx <- which(meth_X$v == 0 & unmeth_X$v == 0)
  mean_nonzero_idx <- mean(meth_X$v[-zero_idx]/unmeth_X$v[-zero_idx])
  meth_X$v[zero_idx] <- 1
  unmeth_X$v[zero_idx] <- ceiling(1/mean_nonzero_idx)

  p <- ncol(meth_X)
  if(verb>0)
    cat(sprintf("\nEstimating on a %d samples collection.\n", nrow(meth_X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }

  ##########  initialization of the methClust model   ###################

  if(sample_init==TRUE){
    if(is.null(NUM_INDICES_START)){
      pre_index_init <- 1:(max(2, min(ceiling(nrow(meth_X)*.05) + 2*K, 25)));
    }else{
      pre_index_init <- 1:NUM_INDICES_START;
    }
    samp_length <- length(pre_index_init);
    index_init <- sample(1:nrow(meth_X),samp_length, replace = FALSE);
  }else{
    if(is.null(NUM_INDICES_START)){
      index_init <- 1:(max(2, min(ceiling(nrow(meth_X)*.05) + 2*K, 25)));
    }else{
      index_init <- 1:NUM_INDICES_START;
    }
  }

  initopics <- meth_tpxinit(meth_X[index_init,], unmeth_X[index_init,],
                            initopics, (K+1), verb, use_squarem=FALSE)
  initopics2 <- initopics[,(K+1)]
  initopics1 <- initopics[,sample(1:K, K-1, replace = FALSE)]
  initopics <- cbind(initopics1, initopics2)

  tpx <- meth_tpxfit(meth_X, unmeth_X, freq = initopics,
                   tol, verb, use_squarem, admix=TRUE, grp=NULL,
                   tmax = 10000, wtol=10^{-4}, qn=100)

  K <- tpx$K
  L <- tpx$L

  ## clean up and out
  if(ord){ worder <- order(slam::col_sums(tpx$omega), decreasing=TRUE) }else{ worder <- 1:K }

  freq=matrix(tpx$freq[,worder], ncol=K, dimnames=list(phrase=dimnames(meth_X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(meth_X)){ rownames(omega) <- dimnames(meth_X)[[1]] }

  ## topic object
  out <- list(K=K, L = L,
              freq=freq, omega=omega,
              meth_X= meth_X, unmeth_X = unmeth_X)
  class(out) <- "topics"
  invisible(out) }

