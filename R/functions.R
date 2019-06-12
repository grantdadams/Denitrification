

#' Water density
#'
#' @description Function to calculate the water density as a function of temperature.
#'
#' @param temp Water temperature (Celcius)
#'
#' @export
#'
#' @examples
#'
#' watdens(29)
#'
watdens<-function(temp){

  t<-temp

  A <- 7.0132e-5
  B <- 7.926295e-3
  C <-  -7.575477e-5
  D<- 7.314701e-7
  E <-  -3.596363e-9
  to<- 3.9818

  dens<- (999.97358- (A*(t-to) + B*(t-to)^2 +C*(t-to)^3 + D*(t-to)^4+E*(t-to)^5) ) -4.873e-3 + 1.708e-4*t - 3.108e-6 * t^2
  dens/1000
}


#' N2 saturation
#'
#' @description Function to calculate the N2 saturation.
#'
#' @param temp Temperature (Celcius)
#' @param salinity Salinity
#' @param bp Barymetric pressure (elevation-corrected mm Hg)
#'
#' @export
#'
#' @examples
#'
#' n2satfc(15, 0, 760)
#'
n2satfc <- function(temp,salinity,bp){
  Ts=log((298.15-temp)/(273.15+temp))

  A0=6.42931
  A1=2.92704
  A2=4.32531
  A3=4.69149
  B0=-0.00744129
  B1=-0.00802566
  B2=-0.0146775

  ln.C = 6.42931 + 2.92704 * Ts + 4.32531 * Ts^2 + 4.69149 * Ts^3 + salinity*(-0.00744129 + -0.00802566 * Ts + -0.0146775 * Ts^2)
  u<-10^(8.10765-(1750.286/(235+temp)))
  C<-exp(ln.C)*((bp-u)/(760-u))
  converted<-C*watdens(temp)*(28.014/1000) # converts N2 from umolN2/kg to gN/m3
  return(converted)
}



#' Ar Saturation
#'
#' @description Function to calculate the ar saturation.
#'
#' @param temp Temperature (Celcius)
#' @param salinity Salinity
#' @param bp Barymetric pressure (elevation-corrected mm Hg)
#'
#' @export
#'
#' @examples
#'
#' arsatfc(15, 0, 760)
#'
arsatfc<-function(temp, salinity, bp){
  Ts=log((298.15-temp)/(273.15+temp))

  A0=2.79150
  A1=3.17609
  A2=4.13116
  A3=4.90379
  B0=-0.00696233
  B1=-0.00766670
  B2=-0.0116888

  ln.C = A0 + A1 * Ts + A2 * Ts^2 + A3 * Ts^3 + salinity*(B0 + B1 * Ts + B2 * Ts^2)
  u<-10^(8.10765-(1750.286/(235+temp)))
  C<-exp(ln.C)*((bp-u)/(760-u))
  converted<-C*watdens(temp)*(39.948/1000)# converts Ar from umol/kg to g/m3
  return(converted)
}


#' Data list
#'
#' @param data
#'
#' @export
#'
create_dataList <- function(data, Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, z = 0.5588){

  # Convert data
  data$arsat <- arsatfc(temp=data$temp, salinity=0, bp=data$bp)
  data$n2sat <- n2satfc(temp=data$temp, salinity=0, bp=data$bp)
  data$n2convert <- data$n2.ar * (28.014 / 1000) / (39.948 / 1000) * data$arsat

  # Subset data
  updata <- data[data$station == up,]
  downdata <- data[data$station == down,]

  # Data objects
  lag <- round(tt / 0.0636573611) # 3 for measured Ditch 2 HRT  ~4.58333/24= 0.1909720833/3 =0.063657

  # Create list of data objects
  data_list <- list()

  # Data size
  data_list$nobs <- nrow(data)
  data_list$nup  <- nrow(updata)
  data_list$ndown <- nrow(downdata)

  # Temperature
  data_list$tempup <- updata$temp[1:(data_list$nup - lag)]
  data_list$tempdown <- downdata$temp[(1+lag):data_list$ndown]

  # N2
  data_list$n2up <- updata$n2convert[1:(data_list$nup - lag)]
  data_list$n2down <- downdata$n2convert[(1+lag):data_list$ndown]

  # BP
  data_list$bpup <- updata$bp[1:(data_list$nup - lag)]
  data_list$bpdown <- downdata$bp[(1+lag):data_list$ndown]

  # N2 equilibrium concentration
  data_list$n2equilup <- n2satfc(data_list$tempup, 0, data_list$bpup )
  data_list$n2equildown <- n2satfc(data_list$tempdown, 0, data_list$bpdown)

  # Light
  data_list$light<-downdata$light

  # Parameters
  data_list$Kmean = Kmean
  data_list$Ksd = Ksd

  return(data_list)

}



#' Estimate
#'
#' @param data data.frame object with columnds for station, temperature, n2.ar, light, bp, z, and dtime, arsat, and n2sat
#' @param model Model number for fitting algorithm.
#'
#' @details Model determines which model to estimate. 3 being the two-station model without N consumption (DN base), and 4 being the two station model with N consumption (DN N consume).
#'
#' @export
#'
#' @examples
#'
#' dataList <- create_dataList(InitialData, , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, z)
fitmod <- function(dataList, model = 3, nChains = 2, niter = 3000, burnin = 1000, verbose = FALSE, warmup){

# Data set up

  dataList$model = model

  # StanFit <- stan('src/nn2_model.stan', data = dataList, iter = niter, chains = nChains, verbose = FALSE, warmup = burnin)


}

