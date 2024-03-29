## -----------------------------------------------------------------------------
library(momos)
# All parameters must be numerical, and all are optional
pars <- list()
# Initial values
pars$Necromasa <- 2140
pars$HLo <- 2250
pars$HSo <- 19150
pars$Co <- 55.40
# Parameters values
pars$Kvl <- 0.2070
pars$Kvs <- 0.00057
pars$fs <- 0.00002
pars$Khl <- 0.0638
pars$Khs <- 0.00077
pars$Khls <- 0.1581
pars$Kmb <- 0.001
pars$Kresp <- 0.1419
pars$Ci <- 50.04
# Time values
pars$from <- 0
pars$to <- 30
pars$by <- 1
out <- momos(params = pars)
print(out)

