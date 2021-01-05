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
  Necromasa <- if("Necromasa" %in% names(params)) params$Necromasa else 2140
  HLo <- if("HLo" %in% names(params)) params$HLo else 2250
  HSo <- if("HSo" %in% names(params)) params$HSo else 19150
  Co <- if("Co" %in% names(params)) params$Co else 55.40

  # Parameters of the MOMOS model
  Kvl <- if("Kvl" %in% names(params)) params$Kvl else 0.2070
  Kvs <- if("Kvs" %in% names(params)) params$Kvs else 0.00057
  fs <- if("fs" %in% names(params)) params$fs else 0.00002
  Khl <- if("Khl" %in% names(params)) params$Khl else 0.0638
  Khs <- if("Khs" %in% names(params)) params$Khs else 0.00077
  Khls <- if("Khls" %in% names(params)) params$Khls else 0.1581
  Kmb <- if("Kmb" %in% names(params)) params$Kmb else 0.001
  Kresp <- if("Kresp" %in% names(params)) params$Kresp else 0.1419
  Ci <- if("Ci" %in% names(params)) params$Ci else 50.04

  # Execution time
  from <- if("from" %in% names(params)) params$from else 0
  to <- if("to" %in% names(params)) params$to else 30
  by <- if("by" %in% names(params)) params$by else 1

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
