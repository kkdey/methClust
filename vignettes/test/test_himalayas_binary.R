
###########  test for himalayas regional asemblage   ################

reg_data <- readRDS("../../../methClust-pages/data/Him_site_x_cell_matrix.rds")
reg_data[reg_data > 0] = 1
