#' The MOMOS function
#'
#' @param params params of the function
#'
#' @return Table with the data
#'
#' @examples
#' params <- list()
#' params$Necromasa <- 2140
#' params$HLo <- 2250
#' params$HSo <- 19150
#' params$Co <- 55.40
#' params$Kvl <- 0.2070
#' params$Kvs <- 0.00057
#' params$fs <- 0.00002
#' params$Khl <- 0.0638
#' params$Khs <- 0.00077
#' params$Khls <- 0.1581
#' params$Kmb <- 0.001
#' params$Kresp <- 0.1419
#' params$Ci <- 50.04
#' params$from <- 0
#' params$to <- 30
#' params$by <- 1
#' out <- momos(params)
#'
#' @export
momos <- function (params = NULL) {

  # Check if params are a list
  if(is.list(params)) {
    for(i in params) {
      if(!is.numeric(i)) {
        message("All elements must be numeric")
        return(NA)
      }
    }
  } else {
    if(!is.null(params)) {
      message("Parameters must be a list")
      return(NA)
    }
  }

  # Initial values
  Necromasa <- ifelse("Necromasa" %in% names(params), params$Necromasa, 2140)
  HLo <- ifelse("HLo" %in% names(params), params$HLo, 2250)
  HSo <- ifelse("HSo" %in% names(params), params$HSo, 19150)
  Co <- ifelse("Co" %in% names(params), params$Co, 55.40)

  # Parameters of the MOMOS model
  Kvl <- ifelse("Kvl" %in% names(params), params$Kvl, 0.2070)
  Kvs <- ifelse("Kvs" %in% names(params), params$Kvs, 0.00057)
  fs <- ifelse("fs" %in% names(params), params$fs, 0.00002)
  Khl <- ifelse("Khl" %in% names(params), params$Khl, 0.0638)
  Khs <- ifelse("Khs" %in% names(params), params$Khs, 0.00077)
  Khls <- ifelse("Khls" %in% names(params), params$Khls, 0.1581)
  Kmb <- ifelse("Kmb" %in% names(params), params$Kmb, 0.001)
  Kresp <- ifelse("Kresp" %in% names(params), params$Kresp, 0.1419)
  Ci <- ifelse("Ci" %in% names(params), params$Ci, 50.04)

  # Execution time
  from <- ifelse("from" %in% names(params), params$from, 0)
  to <- ifelse("to" %in% names(params), params$to, 30)
  by <- ifelse("by" %in% names(params), params$by, 1)

  # Derivs
  derivs <- function(time, y, pars) {
    with (as.list(c(pars, y)), {
      Fvl <- VL * Kvl
      Fvs <- VS * Kvs
      Fmor <- CM * Kmb
      Resp <- CM * CM * Kresp / Co
      Fhl <- HL * Khl
      Fhls <- HL * Khls
      Fhs <- HS * Khs

      dVL <- -Fvl
      dVS <- -Fvs
      dCM <- Fvl + Fvs - Fmor - Resp + Fhl + Fhs
      dHL <- Fmor - Fhl - Fhls
      dHS <- Fhls - Fhs
      dRA <- Resp

      return (list(c(dVL,  dVS, dCM, dHL, dHS, dRA)))
    })
  }

  # Set params
  pars <- c(
    Kvl <- Kvl,
    Kvs <- Kvs,
    fs <- fs,
    Khl <- Khl,
    Khs <- Khs,
    Khls <- Khls,
    Kmb <- Kmb,
    Kresp <- Kresp,
    Ci <- Ci
  )

  Co <- Co
  HLo <- HLo
  HSo <- HSo
  Necromasa <- Necromasa

  # Integral
  y <- c(
    VL = Necromasa * (1 - fs),
    VS = Necromasa * fs,
    CM = Ci,
    HL = HLo,
    HS = HSo,
    RA = 0
  )

  # Execution time
  times <- seq(from = from, to = to, by = by)

  library(deSolve)

  out <- ode(
    func = derivs,
    y = y,
    parms = pars,
    times = times,
    method = "rk4"
  )

  plot(out)

  return(as.data.frame(out))

}

momos_optimize <- function (variable = NULL) {
  # Check if variable is null
  if(is.null(variable)) {
    message("Variable cannot be null")
    return(NA)
  }

  # Check if variable is CM or RA
  if(variable != "CM" && variable != "RA") {
    message("Variable must be CM or RA")
    return(NA)
  }

  # Get simulate data
  simulate_data <- as.data.frame(list(time = out$time, CM = out$CM, RA = out$RA))

  # Get real data
  library("xlsx")
  real_data <- read.xlsx("data/momos.xlsx", sheetIndex = 1)

  data <- data.frame(Time = integer(), Real_Value = double(), Simulate_Value = double(), stringsAsFactors = FALSE)
  str(data)

  # Create data frame with all data
  library("sjmisc")
  for (var in times) {
    real <- NA
    simulate <- NA
    if (!is_empty(real_data[real_data$time == var, ])) {
      ifelse("CM" == variable, real <- real_data[real_data$time == var, ]$CM, real <- real_data[real_data$time == var, ]$RA)
    }
    if (!is_empty(simulate_data[simulate_data$time == var, ])) {
      ifelse("CM" == variable, simulate <- simulate_data[simulate_data$time == var, ]$CM, simulate <- simulate_data[simulate_data$time == var, ]$RA)
    }
    new_record <- data.frame(var, real, simulate)
    names(new_record) <- c("Time", "Real_Value", "Simulate_Value")
    data <- rbind(data, new_record)
  }

  return(data)
}

# Graphs
plot(real_data$time, real_data$RA, type = "o", main = 'MOMOS', col="dark blue", xlab="Time", ylab="RA")
lines(out$time, out$RA, col = "dark red")
legend(25, 500, legend=c("Experimental", "Simulation"), col=c("dark blue", "dark red"), lty=1:1, cex=0.8)
