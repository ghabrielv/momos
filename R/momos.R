#' The MOMOS function
#'
#' @param params params of the function
#'
#' @return Table with the data
#'
#' @examples
#' out <- momos(params)
#'
#' @export
momos <- function (params) {

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

  return(as.data.frame(out))

}
