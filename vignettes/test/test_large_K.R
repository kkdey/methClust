
####################   test  methclust   ###################

meth <- t(read.table("../methClust-pages/data/meth.txt", header = TRUE))
unmeth <- t(read.table("../methClust-pages/data/unmeth_w_zero.txt", header = TRUE))

topics <- meth_topics(meth, unmeth, K=10, tol = 0.1, use_squarem = FALSE)

K = 8
shape=NULL
initopics=NULL
tol=0.1
ord=TRUE
verb=1
sample_init = TRUE
use_squarem=FALSE


methx <- meth_X[index_init,]
unmethx <- unmeth_X[index_init,]

