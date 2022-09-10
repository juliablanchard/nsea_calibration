##########------------------------------------------------------
########## Calculate error function for time series calibration 
##########------------------------------------------------------

getErrorTime <- function(vary,params,effort,dat,sumsquares=T,variable="reproduction_level") {
  
  if (variable=="erepro"){
  params <- setBevertonHolt(params, erepro = vary[1:dim(species_params(params))[1]])
  }
  if (variable=="R_max"){
    params <- setBevertonHolt(params, R_max = 10^vary[1:dim(species_params(params))[1]])
  }
  if (variable=="reproduction_level"){
    params <- setBevertonHolt(params, reproduction_level = vary[1:dim(species_params(params))[1]])
  }
  
  gear_params(params)$catchability<-vary[(1+dim(species_params(params))[1]):(2*dim(species_params(params))[1])]
  
  params<-projectToSteady(params,effort=effort[1,])
  
  #run time-varying effort model though time from the new initial steady state

  simt <- project(params, effort = effort)

  #get yield through time from model:
  
  yield_species<-getYield(simt)
  
  yield_frame <- melt(yield_species)
  
  # leave out spin up 
  
  y<-yield_frame[yield_frame$time >= 1947,]
  
  # disregard zeroes - these were NAs only filled in to run the model   
  
  obs<-dat$Yield/1e9 
  pred<-y$value/1e9
  
  # sum of squared errors, could use  log-scale of predictions and data (could change this or use other error or likelihood options)
  
  sumsq <- sum((log(pred) - log(obs))^2,na.rm=T)
  
  # total relative error
  tre= sum(1-pred/obs,na.rm=T)
  
  # can use a strong penalty on the error to ensure we reach a minimum of 10% of the data (biomass or catch) for each species
  # if(any(pred < 0.1*dat)) discrep <- discrep + 1e10
  
  if (sumsquares==T) error<-sumsq
  if (sumsquares==F) error<-tre
  
  return(error)
  
}

#vary<-c(as.numeric(getReproductionLevel(params2)),gear_params(params)$catchability)
# ## test it
#err<-getErrorTime(vary, params = params2, effort=effort,yields_obs)
# 
# 
# err

##########------------------------------------------------------
########## Plot the outputs of the time series calibration 
##########------------------------------------------------------


### Helpful functions



# function for calibration of reproduction parameters
fastOptim <- function(params,effort,dat,variable,lower_val,upper_val)
{
  # create set of params for the optimisation process
  params_optim <- params
  
  if (variable=="erepro"){
  vary <-  species_params(params_optim)$erepro # variable to explore
  }
  if (variable=="R_max"){
    vary <-  log10(species_params(params_optim)$R_max) # variable to explore
  }
  if (variable=="reproduction_level"){
    vary <-  getReproductionLevel(params_optim) # variable to explore
  }
  
  # set up workers
  noCores <- parallel::detectCores() - 1 # keep some spare core
  cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, varlist = "cl",envir=environment())
  clusterEvalQ(cl, {
    library(mizerExperimental)
    library(optimParallel)
  })

  optim_result <- optimParallel::optimParallel(par=vary,getErrorTime,params=params_optim, effort =effort,dat = dat, method   ="L-BFGS-B", lower=c(rep(lower_val,dim(params_optim@species_params)[1])), upper= c(rep(upper_val,dim(params_optim@species_params)[1])),
                                           parallel=list(loginfo=TRUE, forward=TRUE))
  stopCluster(cl)
  
  
  if (variable=="erepro"){
    params_optim <- setBevertonHolt(params_optim, erepro = optim_result$par)
     }
  if (variable=="R_max"){
    params_optim <- setBevertonHolt(params_optim, R_max = 10^optim_result$par)
     }
  if (variable=="reproduction_level"){
    params_optim <- setBevertonHolt(params_optim, reproduction_level = optim_result$par)
     }
  
  params_optim<- projectToSteady(params_optim, effort = effort[1,])
  sim_optim<- project(params_optim , effort = effort)
  return(sim_optim)
}




# function running tuneParams function in a row for a quick start to a calibration
fastCalib <- function(params, match = F)
{
  params <- calibrateBiomass(params) # changes kappa and rmax
  if(match) params <- matchBiomasses(params) # set rmax to inf and adjust erepro
  params <- steady(params, tol = 0.001)
  sim <- project(params, t_max = 1000)
  return(sim)
}

# removing effort = 1 so always using intial effort and removing /1e6 so everything is in grams
getError2 <- function (vary, params, dat, data_type = "catch", tol = 0.1, 
                       timetorun = 10) 
{
  params@species_params$R_max[] <- 10^vary[1:length(params@species_params$R_max)]
  params <- setParams(params)
  params <- projectToSteady(params, distance_func = distanceSSLogN, 
                            tol = tol, t_max = 200, return_sim = F)
  sim <- project(params, t_max = timetorun, progress_bar = F)
  if (data_type == "SSB") {
    output <- getSSB(sim)[timetorun, ]
  }
  if (data_type == "catch") {
    output <- getYield(sim)[timetorun, ]
  }
  pred <- log(output)
  dat <- log(dat)
  discrep <- pred - dat
  discrep <- (sum(discrep^2))
  return(discrep)
}



