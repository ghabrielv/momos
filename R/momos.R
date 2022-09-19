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
#' params$from <- 1
#' params$to <- 30
#' params$at <- 1
#' out <- momos(params)
#'
#' @export
get_params = function(params = NULL) {

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
  Necromasa <<- ifelse("Necromasa" %in% names(params), params$Necromasa, 2140)
  HLo <<- ifelse("HLo" %in% names(params), params$HLo, 2250)
  HSo <<- ifelse("HSo" %in% names(params), params$HSo, 19150)
  Co <<- ifelse("Co" %in% names(params), params$Co, 55.40)

  # Parameters of the MOMOS model
  Kvl <<- ifelse("Kvl" %in% names(params), params$Kvl, 0.2070)
  Kvs <<- ifelse("Kvs" %in% names(params), params$Kvs, 0.00057)
  fs <<- ifelse("fs" %in% names(params), params$fs, 0.00002)
  Khl <<- ifelse("Khl" %in% names(params), params$Khl, 0.0638)
  Khs <<- ifelse("Khs" %in% names(params), params$Khs, 0.00077)
  Khls <<- ifelse("Khls" %in% names(params), params$Khls, 0.1581)
  Kmb <<- ifelse("Kmb" %in% names(params), params$Kmb, 0.001)
  Kresp <<- ifelse("Kresp" %in% names(params), params$Kresp, 0.1419)
  Ci <<- ifelse("Ci" %in% names(params), params$Ci, 50.04)

  # Execution time
  from <<- ifelse("from" %in% names(params), params$from, 1)
  to <<- ifelse("to" %in% names(params), params$to, 30)
  at <<- ifelse("at" %in% names(params), params$at, 1)

}

momos <- function (params = NULL) {

  # Getting parameters
  get_params(params)

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
  times <- seq(from = from, to = to, by = at)

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

# Calibrating model
calibrate_momos <- function(params = NULL) {
  # Getting parameters
  get_params(params)

  # load libraries
  library(ggplot2) #library for plotting
  library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
  library(deSolve) # library for solving differential equations
  library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
  library(xlsx)

  # Getting experimental data
  experimental_data <<- read.xlsx("data/momos.xlsx", sheetIndex = 1)
  names(experimental_data)=c("time","CM","RA")


  # plot data
  tmp=melt(experimental_data,id.vars=c("time"),variable.name="variables",value.name="values")
  ggplot(data=tmp,aes(x=time,y=values,color=variables))+geom_point(size=3)

  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  parms=c(Kresp=Kresp)
  # fitting
  fitval=nls.lm(par=parms,fn=ssq)
  summary(fitval)

  #### Graphs
  # plot of predicted vs experimental data

  # simulated predicted profile at estimated parameter values
  t=seq(from,to,at)
  parms=as.list(fitval$par)
  out=momos(params = parms)
  out=out[,c("time","CM","RA")]
  outdf=data.frame(out)
  names(outdf)=c("time","CM_pred","RA_pred")

  # Overlay predicted profile with experimental data
  tmppred=melt(outdf,id.var=c("time"),variable.name="variables",value.name="values")
  tmpexp=melt(experimental_data,id.var=c("time"),variable.name="variables",value.name="values")
  p=ggplot(data=tmppred,aes(x=time,y=values,color=variables,linetype=variables))+geom_line()
  p=p+geom_line(data=tmpexp,aes(x=time,y=values,color=variables,linetype=variables))
  p=p+geom_point(data=tmpexp,aes(x=time,y=values,color=variables))
  p=p+scale_linetype_manual(values=c(0,1,0,1))
  p=p+scale_color_manual(values=rep(c("red","blue"),each=2))+theme_bw()
  print(p)
}

ssq=function(parms){
  # Time points for which values is experimental
  # Include the points where data is available
  t=c(seq(from,to,at),experimental_data$time)
  t=sort(unique(t))

  # parameters from the parameter estimation routine
  k1=parms[1]

  # solve the equation
  out_func=momos(list(Kresp=k1))
  out_func=out_func[,c("time","CM","RA")]

  # Filter data that contains time points where data is available
  outdf=data.frame(out_func)
  outdf=outdf[outdf$time %in% experimental_data$time,]

  # Evaluate predicted vs experimental residual
  preddf=melt(outdf,id.var="time",variable.name="variables",value.name="values")
  expdf=melt(experimental_data,id.var="time",variable.name="variables",value.name="values")
  ssqres=preddf$values-expdf$values

  return(ssqres)
}


