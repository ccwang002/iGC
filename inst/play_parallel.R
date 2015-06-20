library(plyr)
library(doMC)

doMC::registerDoMC(cores=4)

system.time(ddply(iris, .(Species), function(x) {
  Sys.sleep(2)
  nrow(x)
}))

system.time(ddply(iris, .(Species), function(x) {
  Sys.sleep(2)
  nrow(x)
}, .parallel = TRUE))


system.time(ret <- as.data.table(mlply(
  sample_desc[, .(Sample, CNA_filepath)],
  function(Sample, CNA_filepath) {
    rep(paste(Sample, CNA_filepath, sep = ' | '), 10)
  },
  .expand = FALSE
)))


system.time(ret_parallel <- as.data.table(mlply(
  sample_desc[, .(Sample, CNA_filepath)],
  function(Sample, CNA_filepath) {
    rep(paste(Sample, CNA_filepath, sep = ' | '), 10)
  },
  .expand = FALSE,
  .parallel = TRUE
)))

all.equal.character(ret_parallel[, c(10), with=FALSE], ret[, c(10), with=FALSE])
