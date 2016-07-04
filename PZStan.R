install.packages("devtools")
library(devtools)
install_github("libbi/RBi",ref="master")
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

df <- data.frame(rdata_PP$P$time,
                 rdata_PP$P$value,
                 rdata_PP$Z$value,
                 rdata_PP$P_obs$value)

ggplot(df, aes(rdata_PP$P$time, color = variable), size = 0.1) +
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

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

lvStanModel <- stan_model(file = "SHO.stan",verbose=TRUE)

lvFit <- sampling(lvStanModel,
                  seed=42,
                  data=list(T = length(rdata_PP$P_obs$value),
                            y = rdata_PP$P_obs$value,
                            k1 = 2.0e2,
                            b  = 2.0e-2,
                            d  = 4.0e-1,
                            k2 = 2.0e1,
                            c  = 4.0e-3,
                            deltaT = rdata_PP$P_obs$time[2] - rdata_PP$P_obs$time[1]
                            ),
                   chains=1)

samples <- extract(lvFit)

gs1 <- qplot(x = samples$mu, y = ..density.., geom = "histogram") + xlab(expression(mu))
gs2 <- qplot(x = samples$sigma, y = ..density.., geom = "histogram") + xlab(expression(sigma))
gs3 <- grid.arrange(gs1, gs2)
ggsave(plot=gs3,filename="diagrams/LvPosteriorStan.png",width=4,height=3)
