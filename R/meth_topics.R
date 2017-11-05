

meth_topics <- function(meth,
                        unmeth,
                        K,
                        shape=NULL,
                        initopics=NULL,
                        tol=0.1,
                        bf=FALSE,
                        kill=2,
                        ord=TRUE,
                        verb=1,
                        use_squarem=TRUE){

  if(dim(meth)[1] != dim(unmeth)[1] | dim(meth)[2] != dim(unmeth)[2]){
    stop("meth and unmeth count matrices must have same dimensions")
  }

  ## sparsify the meth and unmeth count matrices

  meth_X <- CheckCounts(meth)
  unmeth_X <- CheckCounts(unmeth)

  p <- ncol(X)
  if(verb>0)
    cat(sprintf("\nEstimating on a %d samples collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }

  ##########  initialization of the methClust model   ###################

  index_init <- 1:(max(2, min(ceiling(nrow(X)*.05),100)));
  if(sample_init==TRUE){
    samp_length <- length(index_init);
    index_init <- sample(1:nrow(X),samp_length);
  }

  initopics <- meth_tpxinit(X[index_init,], initopics, K[1],
                       shape, verb, nbundles=1, use_squarem=FALSE, init.adapt)

  tpx <- meth_tpxfit(meth_X, unmeth_X, freq = initopics,
                   alpha=shape, tol, verb, use_squarem,
                   admix=TRUE, method_admix=1, grp=NULL,
                   tmax = 10000, wtol=10^{-4}, qn=100)

  K <- tpx$K
  L <- tpx$L

  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) }
  else{ worder <- 1:K }

  freq=matrix(tpx$freq[,worder], ncol=K, dimnames=list(phrase=dimnames(meth_X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(meth_X)){ dimnames(omega)[[1]] <- dimnames(meth_X)[[1]] }

  ## topic object
  out <- list(K=K, L = L,
              theta=theta, omega=omega,
              meth_X= meth_X, unmeth_X = unmeth_X)
  class(out) <- "topics"
  invisible(out) }

