## We only need this if we want to get HEAD.
## install.packages("devtools")
library(devtools)
install_github("sbfnk/RBi",ref="c56b7b7")
install_github("sbfnk/RBi.helpers",ref="master")

rm(list = ls(all.names=TRUE))
unlink(".RData")

library('RBi')
try(detach(package:RBi, unload = TRUE), silent = TRUE)
library(RBi, quietly = TRUE)

library('RBi.helpers')

library('ggplot2', quietly = TRUE)
library('gridExtra', quietly = TRUE)

endTime <- 50

PP <- bi_model("PP.bi")
synthetic_dataset_PP <- bi_generate_dataset(endtime=endTime,
                                            model=PP,
                                            seed="42",
                                            verbose=TRUE,
                                            add_options = list(
                                                noutputs=500))

rdata_PP <- bi_read(synthetic_dataset_PP)

df <- data.frame(rdata_PP$P$nr,
                 rdata_PP$P$value,
                 rdata_PP$Z$value,
                 rdata_PP$P_obs$value)

ggplot(df, aes(rdata_PP$P$nr, y = value, color = variable)) +
    geom_line(aes(y = rdata_PP$P$value, col = "Phyto")) +
    geom_line(aes(y = rdata_PP$Z$value, col = "Zoo")) +
    geom_point(aes(y = rdata_PP$P_obs$value, col = "Phyto Obs"))

ggplot(df, aes(rdata_PP$P$value, y = value, color = variable)) +
    geom_line(aes(y = rdata_PP$Z$value, col = "Zoo"))

synthetic_dataset_PP1 <- bi_generate_dataset(endtime=endTime,
                                             model=PP,
                                             init = list(P = 100, Z=50),
                                             seed="42",
                                             verbose=TRUE,
                                             add_options = list(
                                                 noutputs=500))

rdata_PP1 <- bi_read(synthetic_dataset_PP1)

synthetic_dataset_PP2 <- bi_generate_dataset(endtime=endTime,
                                             model=PP,
                                             init = list(P = 150, Z=25),
                                             seed="42",
                                             verbose=TRUE,
                                             add_options = list(
                                                 noutputs=500))

rdata_PP2 <- bi_read(synthetic_dataset_PP2)

df1 <- data.frame(hare = rdata_PP$P$value,
                  lynx = rdata_PP$Z$value,
                  hare1 = rdata_PP1$P$value,
                  lynx1 = rdata_PP1$Z$value,
                  hare2 = rdata_PP2$P$value,
                  lynx2 = rdata_PP2$Z$value)

ggplot(df1) +
    geom_path(aes(x=df1$hare,  y=df1$lynx, col = "0")) +
    geom_path(aes(x=df1$hare1, y=df1$lynx1, col = "1")) +
    geom_path(aes(x=df1$hare2, y=df1$lynx2, col = "2")) +
    ggtitle("Test 1")
ggsave(filename="diagrams/PPviaLibBi.png",width=3,height=3)

PPInfer <- bi_model("PPInfer.bi")

bi_object_PP <- libbi(client="sample", model=PPInfer, obs = synthetic_dataset_PP)

bi_object_PP$run(add_options = list(
                     "end-time" = endTime,
                     noutputs = endTime,
                     nsamples = 2000,
                     nparticles = 128,
                     seed=42,
                     nthreads = 1),
                 verbose = TRUE,
                 stdoutput_file_name = tempfile(pattern="pmmhoutput", fileext=".txt"))

bi_file_summary(bi_object_PP$result$output_file_name)

mu <- bi_read(bi_object_PP, "mu")$value
g1 <- qplot(x = mu, y = ..density.., geom = "histogram") + xlab(expression(mu))
sigma <- bi_read(bi_object_PP, "sigma")$value
g2 <- qplot(x = sigma, y = ..density.., geom = "histogram") + xlab(expression(sigma))
grid.arrange(g1, g2)

df2 <- data.frame(hareActs = rdata_PP$P$value,
                  hareObs  = rdata_PP$P_obs$value)

ggplot(df, aes(rdata_PP$P$nr, y = value, color = variable)) +
    geom_line(aes(y = rdata_PP$P$value, col = "Phyto")) +
    geom_line(aes(y = rdata_PP$Z$value, col = "Zoo")) +
    geom_point(aes(y = rdata_PP$P_obs$value, col = "Phyto Obs"))
