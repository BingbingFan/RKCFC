# # examples for KCFC, RKCFC, MKCFC and MRKCFC usage----
# # library(RKCFC)
# library(dplyr)
# dt <- pbc2
# # subsset subjects with more than 4 longitudinal measurements
# dt <- dt %>%
#   group_by(id) %>%
#   mutate(freq = n())
# dt <- subset(dt, freq >= 4)
# dt <- data.frame(dt)
# # range(dt$year)
# # options for KCFC &  RKCFC
# opt1  <-  list(methodSelectK = 2,  dataType = "Sparse", lean = TRUE)
# opt2 <-  list(nK = 4, interval = c(0,15), bw = 5, nRegGrid = 181, bwK = 5,
#               basis = fda::create.bspline.basis(c(0,15), nbasis = 18, norder = 4,
#                                                 breaks = seq(0, 15, length.out = 16)))
# # set seed
# set.seed(2023)
# # RKCFC for serBilir with 50 iterations
# RKCFC_2 <- RKCFC(data = dt, id = "id", ys = "serBilir", timevar= "year",
#                  K = 2 ,init.method =  "random",  maxIter = 50,
#                  parallel = F, nmcores = NULL, trace = 2 ,
#                  conv.prop= 0, LOO = FALSE, robust = TRUE,
#                  optns = opt2)
# # KCFC for serBilir with 50 iterations
# KCFC_2 <- RKCFC(data = dt, id = "id", ys = "serBilir", timevar= "year",
#                 K = 2 ,init.method =  "random",  maxIter = 50,
#                 parallel = F, nmcores = NULL, trace = 2 ,
#                 conv.prop= 0, LOO = FALSE, robust = FALSE,
#                 optns = opt1)
#
# table(RKCFC_2$res_clust$g)
# table(KCFC_2$res_clust$g)
# # MKCFC data
# Lt_list1 <- fdapace::MakeFPCAInputs(IDs = dt$id,tVec = dt$year,yVec = dt$serBilir)
# Lt_list2 <- fdapace::MakeFPCAInputs(IDs = dt$id,tVec = dt$year,yVec = dt$albumin)
#
# Ly <- list(Lt_list1$Ly, Lt_list2$Ly)
# Lt <- list(Lt_list1$Lt, Lt_list2$Lt)
#
# # MKCFC and MRKCFC for joint clustering of serBilir and albumin with 50 iterations
# MKCFC_2 <- MKCFC(Ly = Ly, Lt = Lt, data = NULL, id = NULL, ys= NULL, timevar= NULL,
#                   K = 2 , init.method = "random",  maxIter = 50, robust = F,
#                   parallel = FALSE, nmcores = NULL, trace = 2, conv.prop = 0, LOO = F,
#                   optns = list(methodSelectK = 2,  dataType = "Sparse", lean = TRUE) )
#
# MRKCFC_2 <- MKCFC(Ly = Ly, Lt = Lt, data = NULL, id = NULL, ys= NULL, timevar= NULL,
#                     K = 2 , init.method = "random",  maxIter = 50, robust = F,
#                     parallel = FALSE, nmcores = NULL, trace = 2, conv.prop = 0, LOO = F,
#                     optns = list(methodSelectK = 2,  dataType = "Sparse", lean = TRUE) )
#
# table(mod_KCFC$res_clust$g)
# table(mod_RKCFFC$res_clust$g)
#
#
#
