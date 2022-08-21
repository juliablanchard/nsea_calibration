##########------------------------------------------------------
########## Calculate error function for time series calibration 
##########------------------------------------------------------

getErrorTime <- function(vary,params,effort,dat,env=state,tol = 0.1) {
  
  species_params(params)$R_max<-10^(vary[1:dim(species_params)[1]])
  species_params(params)$erepro<-vary[(dim(species_params)[1]+1):(dim(species_params)[1]*2)]
  
  # params <- setParams(params)
  # run to steady state and update params
  # env$params<- projectToSteady(env$params, distance_func = distanceSSLogN,
  #                 tol = tol, t_max = 200,return_sim = F)
  
  #params_steady<- projectToSteady(params, distance_func = distanceSSLogN,
  #                               tol = tol, t_max = 200,return_sim = F)
  
  #run time-varying effort model tthough time with new erepro
  
  simt <- project(params, effort = effort)
  
  # get biomass through time
  biomass <- sweep(simt@n, 3, simt@params@w * simt@params@dw, "*")
  
  #get yield through time from model:
  
  f_gear<-getFMortGear(params,effort)
  yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                              c(1, 2, 3), sum)
  yield_species_gear
  
  yield_species <-apply(yield_species_gear, c(1, 3), sum)
  
  yield_frame <- melt(yield_species)
  
  # leave out spin up    
  y<-yield_frame[yield_frame$time >= 1947,]
  
  # disregard zeroes - these were NAs only filled in to run the model   
  
  obs<-dat$Yield/1e9 
  pred<-y$value/1e9
  
  # sum of squared errors, could use  log-scale of predictions and data (could change this or use other error or likelihood options)
  
  error <- sum((log(pred) - log(obs))^2,na.rm=T)
  
  # can use a strong penalty on the error to ensure we reach a minimum of 10% of the data (biomass or catch) for each species
  # if(any(pred < 0.1*dat)) discrep <- discrep + 1e10
  
  return(error)
  
}

#vary<-c(log10(simt@params@species_params$R_max),simt@params@species_params$erepro,log10(5e11),4)

# ## test it
# err<-getErrorTime(vary = vary, params = params, dat = obsy)
# 
# 
# err

##########------------------------------------------------------
########## Plot the outputs of the time series calibration 
##########------------------------------------------------------


plotFittedTime<-function(sim=simt,obsy=obsy,allSpecies=T,plotSpecies=NULL,startyr=1947,endyr=2019){
  
  # output modelled yields and reshape for plotting - dont know why built-in getYield function doesn't woprk
  
   y <- getYield(simt)
   y <- reshape2::melt(y)
  
  
  # plot these
  
  if (allSpecies ==T) { 
    p<-ggplot(y) + geom_line(data=y, aes(x = time, y = (value), 
                                         colour = sp)) +
      geom_point(data=obsy,aes(x = time, y = (value), 
                               colour = sp),size=0.1) +
      facet_wrap(~sp,scales="free_y") +
      scale_y_continuous(name = "yield [g/year]")  +
      scale_colour_manual(values = sim@params@linecolour) +
      xlim(startyr, endyr)
  }
  
  # look only at  one species at a time and examine on linear 
  if (allSpecies ==F){
    p<-ggplot(y) + geom_line(data=filter(y,sp==plotSpecies), aes(x = time, y = value,colour = sp)) +
      geom_point(data=filter(obsy,sp=="Cod"),aes(x = time, y = value, 
                                                 colour = sp),size=0.6) +
      #facet_wrap(~sp) +
      scale_y_continuous(name = "Yield [g/year]")  +
      scale_colour_manual(values = sim@params@linecolour) +
      xlim(startyr, endyr)
  }
  
  return(p)
}

### Helpful functions



# function for calibration of Rmax
fastOptim <- function(params)
{
  # create set of params for the optimisation process
  params_optim <- params
  vary <-  c(log10(params_optim@species_params$R_max)) # variable to explore
  # set up workers
  noCores <- parallel::detectCores() - 1 # keep some spare core
  cl <- parallel::makeCluster(noCores, setup_timeout = 0.5)
  setDefaultCluster(cl = cl)
  clusterExport(cl, varlist = "cl",envir=environment())
  clusterEvalQ(cl, {
    library(mizerExperimental)
    library(optimParallel)
  })
  optim_result <- optimParallel::optimParallel(par=vary,getErrorTime,params=params_optim, dat = yields_obs, method   ="L-BFGS-B", lower=c(rep(3,dim(params_optim@species_params)[1])), upper= c(rep(15,dim(params_optim@species_params)[1])),
                                               parallel=list(loginfo=TRUE, forward=TRUE))
  stopCluster(cl)
  species_params(params_optim)$R_max <- 10^optim_result$par 
  sim_optim <- project(params_optim, t_max = 2000)
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



