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
