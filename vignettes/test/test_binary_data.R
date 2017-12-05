
#####################  test for binary data   ###################

library(ecostructure)
data <- get(load(system.file("extdata", "HimalayanBirdsData.rda",
                             package = "ecostructure")))
taxonomic_counts <- t(exprs(data))
presence_absence <- taxonomic_counts
presence_absence[presence_absence > 0] = 1
absence_presence <- 1 - presence_absence

meth <- presence_absence
unmeth <- absence_presence

K = 2
shape=NULL
initopics=NULL
tol=0.1
ord=TRUE
verb=1
sample_init = TRUE
NUM_INDICES_START = NULL
use_squarem=FALSE

methx <- meth_X[index_init,]
unmethx <- unmeth_X[index_init,]

meth_X <- methx
unmeth_X <- unmethx
inifreq <- initopics

freq <- inifreq

grp=NULL
tmax = 10000
wtol=10^{-4}
qn=100
verb=FALSE
nef=TRUE
wtol=10^{-5}

tmax=1000
