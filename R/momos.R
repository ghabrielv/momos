#' The MOMOS function
#'
#' @param params params of the function
#'
#' @return Table with the data
#'
#' @examples
#' out <- momos(pars)
#'
#' @export
momos <- function (pars) {

  derivs <- function(time, y, pars) {
    with (as.list(c(pars, y)), {
      Fvl<- VL*Kvl
      Fvs<- VS*Kvs
      Fmor<- CM*Kmb
      Resp<- CM*CM*Kresp/Co
      Fhl<- HL*Khl
      Fhls<- HL*Khls
      Fhs<-HS*Khs

      dVL<- -Fvl
      dVS<- -Fvs
      dCM<- Fvl+Fvs-Fmor-Resp+Fhl+Fhs
      dHL<- Fmor-Fhl-Fhls
      dHS<- Fhls-Fhs
      dRA<- Resp

      return (list(c(dVL,  dVS, dCM, dHL, dHS, dRA)))
    })
  }

  pars<- c( Kvl<- 0.2070,
            Kvs<- 0.00057,
            fs<- 0.00002,
            Khl<-0.0638,
            Khs<-0.00077,
            Khls<-0.1581,
            Kmb<-0.001,
            Kresp<-0.1419,
            Ci<-50.04)

  Co<-55.40
  HLo<-2250
  HSo<-19150
  Necromasa<-2140

  y <- c(VL=Necromasa*(1-fs),
         VS=Necromasa*fs,
         CM=Ci,
         HL=HLo,
         HS=HSo,
         RA = 0 )

  times <- seq(0, 30, by = 1)

  library(deSolve)

  out <- ode(func = derivs, y = y,
             parms = pars, times = times, method="rk4")

  return(as.data.frame(out))

}
