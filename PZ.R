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

T <- 50

LV <- bi_model("LV.bi")
synthetic_dataset <- bi_generate_dataset(endtime=T,
                                         model=LV,
                                         seed="42")

rdata <- bi_read(synthetic_dataset)

df <- data.frame(rdata$y_H$nr, rdata$y_H$value, rdata$y_L$value)
ggplot(df, aes(rdata$H$nr, y = value, color = variable)) +
    geom_line(aes(y = rdata$y_H$value, col = "y1")) +
    geom_line(aes(y = rdata$y_L$value, col = "y2"))

bi_object <- libbi(client="sample", model=LV, obs = synthetic_dataset)

bi_object$run(add_options = list(
                  "end-time" = T,
                  noutputs = T,
                  nsamples = 128,
                  nparticles = 128,
                  seed=42,
                  nthreads = 1),
              verbose = TRUE,
              stdoutput_file_name = tempfile(pattern="pmmhoutput", fileext=".txt"))

bi_file_summary(bi_object$result$output_file_name)

mu <- bi_read(bi_object, "mu")$value

## The Working Model

PZ <- bi_model("PZ.bi")
synthetic_dataset_PZ <- bi_generate_dataset(endtime=T,
                                            model=PZ,
                                            seed="42",
                                            verbose=TRUE,
                                            add_options = list(
                                                noutputs=500))

rdata_PZ <- bi_read(synthetic_dataset_PZ)

df <- data.frame(rdata_PZ$P$nr,
                 rdata_PZ$P$value,
                 rdata_PZ$Z$value,
                 rdata_PZ$P_obs$value)

ggplot(df, aes(rdata_PZ$P$nr, y = value, color = variable)) +
    geom_line(aes(y = rdata_PZ$P$value, col = "Phyto")) +
    geom_line(aes(y = rdata_PZ$Z$value, col = "Zoo")) +
    geom_point(aes(y = rdata_PZ$P_obs$value, col = "Phyto Obs"))

bi_object_PZ <- libbi(client="sample", model=PZ, obs = synthetic_dataset_PZ)

bi_object_PZ$run(add_options = list(
                     "end-time" = T,
                     noutputs = T,
                     nsamples = 2^13,
                     nparticles = 2^7,
                     seed=42,
                     nthreads = 1),
                 verbose = TRUE,
                 stdoutput_file_name = tempfile(pattern="pmmhoutput", fileext=".txt"))

bi_file_summary(bi_object_PZ$result$output_file_name)

mu_PZ <- bi_read(bi_object_PZ, "mu")$value
g1_PZ <- qplot(x = mu_PZ, y = ..density.., geom = "histogram") + xlab(expression(mu_PZ))
sigma_PZ <- bi_read(bi_object_PZ, "sigma")$value
g2_PZ<- qplot(x = sigma_PZ, y = ..density.., geom = "histogram") + xlab(expression(sigma_PZ))
grid.arrange(g1_PZ, g2_PZ)

## The Working Model

PP <- bi_model("PP.bi")
synthetic_dataset_PP <- bi_generate_dataset(endtime=T,
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
