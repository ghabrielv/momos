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
#' momos(params)
#'
#' @export
get_params = function(params = NULL) {

  # Check if params are a list
  if(is.list(params)) {
    for(i in params) {
      if(!is.numeric(i)) {
        message("Todos los elementos deben ser numericos")
        return(NA)
      }
    }
  } else {
    if(!is.null(params)) {
      message("Los parametros deben ser una lista")
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

calculate_momos <- function (params = NULL) {

  # Getting parameters
  get_params(params)

  # Checking parameters
  if(is.na(get_params(params))) {
    return(NA)
  }

  # Derivative
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

  # Create function y
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

  # Differential Equations
  out <- ode(
    func = derivs,
    y = y,
    parms = pars,
    times = times,
    method = "rk4"
  )

  out=out[,c("time","CM","RA")]

  return(as.data.frame(out))

}

# Fitting parameters
calibrate_momos <- function(params = NULL) {
  # Getting parameters
  get_params(params)

  # Checking parameters
  if(is.na(get_params(params))) {
    return(NA)
  }

  # Getting experimental data
  print("Desea cargar los datos experimentales? (s/n) (nota: si responde s debera seleccionar un archivo")
  is_file = readline();
  if(tolower(is_file) == "s") {
    experimental_data <<- read_excel(file.choose(), sheet = 1)
  } else {
    experimental_data_url <- "https://github.com/ghabrielv/momos/raw/master/data/momos.xlsx"
    GET(experimental_data_url, write_disk(file <- tempfile(fileext = ".xlsx")))
    experimental_data <<- read_excel(file, sheet = 1)
  }

  if(length(experimental_data) > 3) {
    stop("El formato de los datos experimentales son incorrectos, solo debe tener 3 columnas: time, CM, RA")
  }
  names(experimental_data)=c("time","CM_experimental","RA_experimental")

  # parameter fitting using levenberg marquart algorithm
  # initial guess for parameters
  parms=c(Kresp=Kresp)
  # fitting
  fitval=nls.lm(par=parms,fn=ssq)
  print("============ DATOS DE CALIBRACION =============")
  print(summary(fitval))

  # simulated predicted profile at estimated parameter values
  t=seq(from,to,at)
  parms=as.list(fitval$par)
  out=calculate_momos(params = parms)
  out=out[,c("time","CM","RA")]
  outdf=data.frame(out)
  names(outdf)=c("time","CM_calibrated","RA_calibrated")

  print("============ DATOS DE MOMOS CALIBRADO =============")
  print(outdf)

  return(outdf)
}

ssq=function(parms){
  # Time points for which values is experimental
  # Include the points where data is available
  t=c(seq(from,to,at),experimental_data$time)
  t=sort(unique(t))

  # parameters from the parameter estimation routine
  k1=parms[1]

  # solve the equation
  out_func=calculate_momos(list(Kresp=k1))
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

graph_momos <- function(){
  # Prepare data for graphing
  exp_data=experimental_data
  names(exp_data)=c("time","CM_experimental","RA_experimental")

  simulated_data=out_simulated
  names(simulated_data)=c("time","CM_simulated","RA_simulated")

  # Overlay calibrated profile with experimental data and simulated data
  tmp_calibrated=melt(out_calibrated,id.var=c("time"),variable.name="variables",value.name="values")
  tmp_experimental=melt(exp_data,id.var=c("time"),variable.name="variables",value.name="values")
  tmp_simulated=melt(simulated_data,id.var=c("time"),variable.name="variables",value.name="values")

  p=ggplot(data=tmp_calibrated,aes(x=time,y=values,color=variables,linetype=variables))+geom_line()
  #p=p+geom_line(data=tmp_experimental,aes(x=time,y=values,color=variables,linetype=variables)) # make lineal to experimental data
  p=p+geom_point(data=tmp_experimental,aes(x=time,y=values,color=variables))
  p=p+geom_line(data=tmp_simulated,aes(x=time,y=values,color=variables,linetype=variables))
  print(p)
}

momos <- function(params = NULL){

  # Getting parameters
  get_params(params)

  # Checking parameters
  if(is.na(get_params(params))) {
    return(NA)
  }

  # load libraries
  #require(ggplot2) #library for plotting
  #require(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
  #require(deSolve) # library for solving differential equations
  #require(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
  #require(xlsx) # library for read files with xlsx format

  # call functions
  out_simulated <<- calculate_momos(params)
  print("=========== DATOS DE MOMOS SIMULADO ===========")
  print(out_simulated)
  print("============ CALIBRANDO EL MODELO =============")
  out_calibrated <<- calibrate_momos(params)
  print("=========== DATOS DE MOMOS EXPERIMENTALES ===========")
  print(experimental_data)
  print("============ GRAFICANDO EL MODELO =============")
  graph_momos()

}
