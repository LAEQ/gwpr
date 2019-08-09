# source("R/bandwith_option.R")
#
# load(file = "data/us_data/Data.Rdata")
# load(file = "data/us_data/us_dmat.rda")
#
# equation <- gsp ~ pcap + pc + emp + unemp
# sequence <- seq(1, 48, 1)
# index <- c("state", "year")
#
# kernel <- "bisquare"
# effect <- "individual"
# model <- "within"
# adaptive <- TRUE
#
#
# dmat <- compute_dmat(coordinates(USStates), p = 2, longlat = F)
# data_avg <- data_preparation(Produc, equation, "state")
# bandwidth <- bandwidth_optimisation(equation, data_avg , dmat, sequence, kernel = "bisquare", adaptive = TRUE)
#
#
# result_QX_QY <- compute_QX_QY(Produc, equation, index, model, effect)
# result_gwpr <- compute_gwpr(result_QX_QY, dmat, bandwidth[[1]], kernel, adaptive)
#
#
# result_shapefile <- new('ShapeFile')
# result_shapefile@LocalR2Mat <- compute_localR2Mat(result_QX_QY, result_gwpr)
# result_shapefile@R2 <- compute_R2(result_QX_QY, result_gwpr)
# result_SEs_TVs <- compute_std_errors_T_values(result_QX_QY, result_gwpr)
#
# result_shapefile@SEsMat <- result_SEs_TVs$SEsMat
# result_shapefile@TVsMat <- result_SEs_TVs$TVsMat
#
#
# result_final <- compute_shapefile(USStates, result_shapefile, result_gwpr, result_QX_QY)
#
#
# Plots_2 <- compute_gpwr_plot(result_final[[2]], P="0.05")
#
#
# ggarrange(plotlist = list(Plots_2[[2]][[1]], Plots_2[[2]][[2]], Plots_2[[2]][[3]], Plots_2[[2]][[4]]), ncol = 2, nrow = 2)
#
#
# result <- gwpr(  USStates, Produc, equation,
#                  index, bandwidth, us_dmat, kernel = kernel, effect = effect, model = model, adaptive = adaptive
# )
#
# nPlots <- compute_gpwr_plot(result[[2]], P="0.05")
# ggarrange(plotlist = list(Plots[[2]][[1]], Plots[[2]][[2]], Plots[[2]][[3]], Plots[[2]][[4]]), ncol = 2, nrow = 2)
