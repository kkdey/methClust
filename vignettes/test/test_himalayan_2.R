
###############  test on himalayan data (count based)  #####################


reg_data <- readRDS("../../../methClust-pages/data/Him_site_x_cell_matrix.rds")

meth <- reg_data
reg_data_2 <- reg_data
reg_data_2[reg_data_2 > 0] = 1
unmeth <- 1 - reg_data_2

reg_data[reg_data > 0] = 1

K = 2
shape=NULL
initopics=NULL
tol=10
ord=TRUE
verb=1
sample_init = TRUE
NUM_INDICES_START = NULL
use_squarem=FALSE



freq <- initopics

grp=NULL
tmax = 10000
wtol=10^{-4}
qn=100
verb=FALSE
nef=TRUE
wtol=10^{-5}


