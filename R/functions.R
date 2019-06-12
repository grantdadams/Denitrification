

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
#' @description Creates a list of data objects to be used by stan
#'
#' @param data data.frame object with columnds for station, temperature, n2.ar, light, bp, z, and dtime, arsat, and n2sat
#' @param Kmean Mean of normal prior distribution for K600
#' @param Ksd Sd of normal prior distribution for K600
#' @param up Name indicating the up-river station name
#' @param down Name indicating the down-river station name
#' @param tt Time between stations
#' @param z Depth
#'
#' @export
#'
create_dataList <- function(data, Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, z = 0.5588){



  # Convert data
  data$arsat <- arsatfc(temp=data$temp, salinity=0, bp=data$bp)
  data$n2sat <- n2satfc(temp=data$temp, salinity=0, bp=data$bp)
  data$n2convert <- data$n2.ar * (28.014 / 1000) / (39.948 / 1000) * data$arsat
  data$dtime<-chron::chron(dates=as.character(data$date), times=as.character(data$time))

  # Subset data
  updata <- data[data$station == up,]
  downdata <- data[data$station == down,]

  # Create list of data objects
  data_list <- list()

  # Data size
  data_list$nobs <- nrow(data)
  nup  <- nrow(updata)
  ndown <- nrow(downdata)

  # Data objects
  data_list$lag <- round(tt / 0.0636573611) # 3 for measured Ditch 2 HRT  ~4.58333/24= 0.1909720833/3 =0.063657

  # Data-subsets
  up_sub <- 1:(nup - data_list$lag)
  down_sub <- (1 + data_list$lag):ndown

  # Temperature
  data_list$tempup <- updata$temp[up_sub]
  data_list$tempdown <- downdata$temp[down_sub]

  # N2
  data_list$n2up <- updata$n2convert[up_sub]
  data_list$n2down <- downdata$n2convert[down_sub]

  # BP
  data_list$bpup <- updata$bp[up_sub]
  data_list$bpdown <- downdata$bp[down_sub]

  # N2 equilibrium concentration
  data_list$n2equilup <- n2satfc(data_list$tempup, 0, data_list$bpup )
  data_list$n2equildown <- n2satfc(data_list$tempdown, 0, data_list$bpdown)

  # Time object
  data_list$timedown <- downdata$dtime[down_sub]
  data_list$timeup<-updata$dtime[up_sub]

  # Light
  data_list$light<-downdata$light

  # Parameters
  data_list$Kmean = Kmean
  data_list$Ksd = Ksd

  # Depth and time
  data_list$z = z
  data_list$tt = tt

  data_list$nup  <- length(up_sub)
  data_list$ndown <- length(down_sub)

  return(data_list)
}



#' Estimate model
#'
#' @description fits a stan model and returns a stan model object
#'
#' @param dataList dataList created by \code{create_dataList}
#' @param model Model number for fitting algorithm.
#' @param nChains Number of chains to run
#' @param niter Number of iterations to run for each chain (including burnin)
#' @param burnin A positive integer specifying the number of warmup (aka burnin) iterations per chain.
#' @param verbose TRUE or FALSE: flag indicating whether to print intermediate output from Stan on the console, which might be helpful for model debugging.
#'
#' @details Model determines which model to estimate. 3 being the two-station model without N consumption (DN base), and 4 being the two station model with N consumption (DN N consume).
#'
#' @export
#'
#' @examples
#'
#' # Run model Eq. 5
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, z = 0.5588)
#' mod <- fitmod(dataList, model = 3)
#'
#' # Run model Eq. 6
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.1909720833, z = 0.5588)
#' mod <- fitmod(dataList, model = 3)
#' plotmod(mod, dataList = dataList, model = 3, file = NULL)
#'
#' Run model Eq. 7
#' data(InitialData)
#' dataList <- create_dataList(InitialData , Kmean = 4.03, Ksd = 4.0, up = "up1", down = "down1", tt = 0.19097290833, z = 0.5588)
#' mod2 <- fitmod(dataList, model = 4, verbose = FALSE)
#' plotmod(mod2, dataList = dataList, model = 4, file = NULL)
#'
fitmod <- function(dataList, model = 3, nChains = 2, niter = 3000, burnin = 1000, verbose = FALSE){

  # Model set up
  dataList$mod = model
  if(model == 3){
    dataList$nparam = 2
    params_names <- c("DN", "K600", "sigma2", "logPost")
  }
  if(model == 4){
    dataList$nparam = 3
    params_names <- c("DN", "K600", "Nfix", "sigma2", "logPost")
  }

  # Get model
  stan_directory <- system.file("executables",package="Denitrification")
  old_wd <- getwd()
  setwd(stan_directory)

  # Estimate
  StanFit <- rstan::stan("nn2_model.stan", data = dataList, iter = niter, chains = nChains, verbose = verbose, warmup = burnin)
  setwd(old_wd)

  return(StanFit)
}


#' Plot and diagnost
#'
#' @param StanFit Standmodel output by \code{fitmod}
#' @param dataList dataList created by \code{create_dataList}
#' @param model Model number for fitting algorithm.
#' @param file filname to save the parameter estimates. Will not save if NULL.
#'
#' @details Model determines which model to estimate. 3 being the two-station model without N consumption (DN base), and 4 being the two station model with N consumption (DN N consume).
#' @export
#'
#' @examples
plotmod <- function(StanFit, dataList, model = 3, file = NULL){

  # Results of parameters to return
  results <- summary(StanFit)$summary
  rows_sub <- c(grep( "n2hat", rownames(results)), grep( "n2pred", rownames(results)))
  results <- results[-rows_sub,]

  if(model == 3){
    params_names <- c("DN", "K600", "sigma2", "logPost")
  }
  if(model == 4){
    params_names <- c("DN", "K600", "Nfix", "sigma2", "logPost")
  }

  rownames(results) <- params_names
  print(results)

  if(!is.null(file)){
    write.csv(results, file = paste0(file, ".csv"))
  }

  # Get list of returned objects
  returned_objects <- extract(StanFit)

  n2hat <- returned_objects$n2hat
  n2pred <- returned_objects$n2pred


  # Plot it
  for(i in 1:length(1+!is.null(file))){
    if(!is.null(file)){
    png(filename = paste0(file, ".png"), height = 4, width = 6, units = "in")
    }
    plot(dataList$timedown,colMeans(n2hat), type="l",xlab="Time", ylab="Nitrogen  (mg/L)", ylim=c(11,15),cex.lab=1.5, cex.axis=1.5, lwd=3, col="black" )
    points(dataList$timedown,dataList$n2down)

    if(!is.null(file)){
     dev.off()
    }
  }


  return(results)
}

