install.packages("devtools")
library(devtools)
install_github("sbfnk/RBi",ref="master")
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

ggplot(df, aes(rdata_PP$P$nr, y = Population, color = variable), size = 0.1) +
    geom_line(aes(y = rdata_PP$P$value, col = "Hare"), size = 0.1) +
    geom_line(aes(y = rdata_PP$Z$value, col = "Lynx"), size = 0.1) +
    geom_point(aes(y = rdata_PP$P_obs$value, col = "Observations"), size = 0.1) +
    theme(legend.position="none") +
    ggtitle("Example Data") +
    xlab("Days") +
    theme(axis.text=element_text(size=4),
          axis.title=element_text(size=6,face="bold")) +
    theme(plot.title = element_text(size=10))
ggsave(filename="diagrams/LVdata.png",width=4,height=3)

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
    geom_path(aes(x=df1$hare,  y=df1$lynx, col = "0"), size = 0.1) +
    geom_path(aes(x=df1$hare1, y=df1$lynx1, col = "1"), size = 0.1) +
    geom_path(aes(x=df1$hare2, y=df1$lynx2, col = "2"), size = 0.1) +
    theme(legend.position="none") +
    ggtitle("Phase Space") +
    xlab("Hare") +
    ylab("Lynx") +
    theme(axis.text=element_text(size=4),
          axis.title=element_text(size=6,face="bold")) +
    theme(plot.title = element_text(size=10))
ggsave(filename="diagrams/PPviaLibBi.png",width=4,height=3)

PPInfer <- bi_model("PPInfer.bi")

bi_object_PP <- libbi(client="sample", model=PPInfer, obs = synthetic_dataset_PP)

bi_object_PP$run(add_options = list(
                     "end-time" = endTime,
                     noutputs = endTime,
                     nsamples = 4000,
                     nparticles = 128,
                     seed=42,
                     nthreads = 1),
                 ## verbose = TRUE,
                 stdoutput_file_name = tempfile(pattern="pmmhoutput", fileext=".txt"))

bi_file_summary(bi_object_PP$result$output_file_name)

mu <- bi_read(bi_object_PP, "mu")$value
g1 <- qplot(x = mu[2001:4000], y = ..density.., geom = "histogram") + xlab(expression(mu))
sigma <- bi_read(bi_object_PP, "sigma")$value
g2 <- qplot(x = sigma[2001:4000], y = ..density.., geom = "histogram") + xlab(expression(sigma))
g3 <- grid.arrange(g1, g2)
ggsave(plot=g3,filename="diagrams/LvPosterior.png",width=4,height=3)


df2 <- data.frame(hareActs = rdata_PP$P$value,
                  hareObs  = rdata_PP$P_obs$value)

ggplot(df, aes(rdata_PP$P$nr, y = value, color = variable)) +
    geom_line(aes(y = rdata_PP$P$value, col = "Phyto")) +
    geom_line(aes(y = rdata_PP$Z$value, col = "Zoo")) +
    geom_point(aes(y = rdata_PP$P_obs$value, col = "Phyto Obs"))

ln_alpha <- bi_read(bi_object_PP, "ln_alpha")$value

P <- matrix(bi_read(bi_object_PP, "P")$value,nrow=51,byrow=TRUE)
Z <- matrix(bi_read(bi_object_PP, "Z")$value,nrow=51,byrow=TRUE)

data50 <- bi_generate_dataset(endtime=endTime,
                              model=PP,
                              seed="42",
                              verbose=TRUE,
                              add_options = list(
                                  noutputs=50))

rdata50 <- bi_read(data50)

df3 <- data.frame(days = c(1:51), hares = rowMeans(P), lynxes = rowMeans(Z),
                                  actHs = rdata50$P$value, actLs = rdata50$Z$value)


ggplot(df3) +
    geom_line(aes(x = days, y = hares, col = "Est Phyto")) +
    geom_line(aes(x = days, y = lynxes, col = "Est Zoo")) +
    geom_line(aes(x = days, y = actHs, col = "Act Phyto")) +
    geom_line(aes(x = days, y = actLs, col = "Act Zoo"))
