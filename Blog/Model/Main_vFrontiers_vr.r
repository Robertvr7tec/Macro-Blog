

# Prep --------------------------------------------------------------------
library(dplyr)

#cloud
 root<-"C:/Users/Roberto/Desktop/Maestria/Trimestre 6/Macro/Blog/Experimental Design/"
 Number.Cores <- parallel::detectCores() - 2


  dir.exp.inputs<-paste(root,"RDM Inputs\\",sep="")
  Limits.File<-"Limits.csv"
  Policies.File<-"Policies.csv"
  Climate.File<-"Climate.csv"
  sample.size<-20
  Policy.Switch<-TRUE
  Climate.Switch<-TRUE
  source(paste(dir.exp.inputs,"create_experiment_function-1.r",sep=""))
  Exp.design<-exp.design.table(dir.exp.inputs,Limits.File,sample.size,Policies.File,Policy.Switch,Climate.File,Climate.Switch)
  write.csv(Exp.design, paste(dir.exp.inputs, "Exp.design.csv", sep=""), row.names=FALSE)

#Source Model
dir.model<-paste(root,"TechChange Model/",sep="")
model.version<-"InternationalGreenTechChangeModel_10_22_2015.r"
source(paste(dir.model,model.version,sep=""))

  
  
  #Source Experimental Design
  
  dir.exp<-paste(root,"RDM Inputs\\",sep="")
  experiment.version<-"Exp.design.csv"
  Exp.design<-read.csv(paste(dir.exp,experiment.version,sep=""))
  Exp.design<- Exp.design %>% filter(policy.name =="Nordhaus+R&DGreenClimateFund" | policy.name =="FWA")
  detach("package:dplyr", unload = TRUE)
#Define directory to print output files
  dir.harness<-paste(root,"RDM Harness//",sep="")
#Clean output folder
  csv_files <- list.files(dir.harness, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) > 0) {
    removed <- file.remove(csv_files)
    print(removed)
  } else {
    message("No CSV files found to delete.")
  }
  

  
  
  
# Run Model in Parallel ---------------------------------------------------


#Set up parallel environment
#Run Model in Parallel
  library(parallel)
  library(deSolve)
  library(optimx)
  library(tictoc)
  
  # Define policy config function BEFORE parLapply (in global env so workers can use it)
  get_policy_config <- function(policy.name,
                                ceN.0, ceN.m, ceN.M,
                                ceS.0, ceS.m, ceS.M,
                                tN.0, tN.m, tN.M,
                                sN.0, sN.m, sN.M,
                                tS.0, tS.m, tS.M,
                                sS.0, sS.m, sS.M) {
    switch(policy.name,
           "Nordhaus" = list(
             start = c(ceN.0, ceS.0),
             lower = c(ceN.m, ceS.m),
             upper = c(ceN.M, ceS.M),
             parscale = c(10, 10)
           ),
           "Nordhauds+TechnologyPolicy" = list(
             start = c(ceN.0, ceS.0, tN.0, sN.0),
             lower = c(ceN.m, ceS.m, tN.m, sN.m),
             upper = c(ceN.M, ceS.M, tN.M, sN.M),
             parscale = c(10, 10, 3, 3)
           ),
           "Nordhaus+TraditionalGreenClimateFund" = list(
             start = c(ceN.0, tN.0, sN.0, tS.0),
             lower = c(ceN.m, tN.m, sN.m, tS.m),
             upper = c(ceN.M, tN.M, sN.M, tS.M),
             parscale = c(10, 3, 3, 3)
           ),
           "Nordhaus+R&DGreenClimateFund" = list(
             start = c(ceN.0, tN.0, sN.0, tS.0, sS.0),
             lower = c(ceN.m, tN.m, sN.m, tS.m, sS.m),
             upper = c(ceN.M, tN.M, sN.M, tS.M, sS.M),
             parscale = c(10, 3, 3, 3, 3)
           ),
           "Nordhaus+TechnologyPolicy.Both" = list(
             start = c(ceN.0, ceS.0, tN.0, sN.0, tS.0, sS.0),
             lower = c(ceN.m, ceS.m, tN.m, sN.m, tS.m, sS.m),
             upper = c(ceN.M, ceS.M, tN.M, sN.M, tS.M, sS.M),
             parscale = c(10, 10, 3, 3, 3, 3)
           ),
           "Nordhaus+TraditionalGreenClimateFund+R&DS" = list(
             start = c(ceN.0, tN.0, sN.0, tS.0, sS.0),
             lower = c(ceN.m, tN.m, sN.m, tS.m, sS.m),
             upper = c(ceN.M, tN.M, sN.M, tS.M, sS.M),
             parscale = c(10, 3, 3, 3, 3)
           ),
           "Nordhaus+CoR&DGreenClimateFund" = list(
             start = c(ceN.0, tN.0, sN.0, sS.0),
             lower = c(ceN.m, tN.m, sN.m, sS.m),
             upper = c(ceN.M, tN.M, sN.M, sS.M),
             parscale = c(10, 3, 3, 3)
           ),
           "Nordhaus+CoR&DGreenClimateFund+TecS" = list(
             start = c(ceN.0, tN.0, sN.0, tS.0, sS.0),
             lower = c(ceN.m, tN.m, sN.m, tS.m, sS.m),
             upper = c(ceN.M, tN.M, sN.M, tS.M, sS.M),
             parscale = c(10, 3, 3, 3, 3)
           ),
           NULL
    )
  }
  
  # Split dataframe into a list of rows (each row = one experiment)
  exp_list <- split(Exp.design, seq_len(nrow(Exp.design)))
  
  # Create a cluster with specified number of cores
  cl <- makeCluster(Number.Cores, type = "PSOCK")  # use FORK if on macOS/Linux
  
  # Automatically stop cluster on exit/error
  on.exit(stopCluster(cl))
  
  # Export all necessary functions and variables to workers
  clusterExport(cl, varlist = c("TechChangeMod", "dir.harness", "dede", "lagderiv", "lagvalue", "optimx", "get_policy_config"), envir = environment())
  
  tic("Parallel model run")
  results <- parLapply(cl, exp_list, function(x) {
    # Extract parameters
    params <- c(
      S.0 = x$S.0, TimeStep = x$TimeStep, EndTime = x$EndTime,
      alfa = x$alfa, epsilon = x$epsilon,
      Gamma.re = x$Gamma.re, k.re = x$k.re,
      Gamma.ce = x$Gamma.ce, k.ce = x$k.ce,
      Eta.re = x$Eta.re, Eta.ce = x$Eta.ce,
      Nu.re = x$Nu.re, Nu.ce = x$Nu.ce,
      qsi = x$qsi, Delta.S = x$Delta.S,
      Delta.Temp.Disaster = x$Delta.Temp.Disaster,
      Beta.Delta.Temp = x$Beta.Delta.Temp,
      CO2.base = x$CO2.base, CO2.Disaster = x$CO2.Disaster,
      labor.growth_N = x$labor.growth_N,
      labor.growth_S = x$labor.growth_S,
      lambda.S = x$lambda.S,
      sigma.utility = x$sigma.utility,
      rho = x$rho,
      Yre.0_N = x$Yre.0_N, Yce.0_N = x$Yce.0_N,
      Yre.0_S = x$Yre.0_S, Yce.0_S = x$Yce.0_S,
      size.factor = x$size.factor, Run.ID = x$Run.ID,
      policy.name = as.character(x$policy.name),
      dir.harness = dir.harness
    )
    
    # Tax/subsidy bounds
    ceN.0 <- if (x$Climate.Model %in% c("GFDL-ESM2G", "GFDL-ESM2M")) 0.30 else if (x$epsilon < 8) 0.25 else if (x$epsilon < 9) 0.20 else if (x$epsilon < 10) 0.15 else 0.10
    ceS.0 <- ceN.0
    ceN.m <- ceS.m <- 0.05
    ceN.M <- ceS.M <- 0.5
    tN.0 <- 0.10; tN.m <- 0.01; tN.M <- 0.15
    tS.0 <- 0.05; tS.m <- 0.01; tS.M <- 0.15
    sN.0 <- 2.0;  sN.m <- 0.5;  sN.M <- 3.0
    sS.0 <- 0.5;  sS.m <- 0.01; sS.M <- 3.0
    
    policy <- x$policy.name
    
    result <- tryCatch({
      if (policy == "FWA") {
        TechChangeMod(c(0.0, 0.0, 0.0, 0.0), params)
      } else {
        cfg <- get_policy_config(policy, ceN.0, ceN.m, ceN.M, ceS.0, ceS.m, ceS.M,
                                 tN.0, tN.m, tN.M, sN.0, sN.m, sN.M,
                                 tS.0, tS.m, tS.M, sS.0, sS.m, sS.M)
        
        if (!is.null(cfg)) {
          optimx(par = cfg$start, fn = TechChangeMod,
                 lower = cfg$lower, upper = cfg$upper,
                 method = "L-BFGS-B",
                 control = list(fnscale = -1,
                                ndeps = rep(0.01, length(cfg$start)),
                                parscale = cfg$parscale,
                                maxit = 2),
                 params = params)
        } else {
          paste("Unknown policy:", policy)
        }
      }
    }, error = function(e) list(error = e$message))
    
    return(list(run_id = x$Run.ID, result = result))
  })
  toc()
   #Stop cluster
   stopCluster(cl)



   
   
### We dont have the file to continue code from here.
   
   
   
   
## =====================================================================================================
## This section reads the output of simulations and reshapes into a format ready for scenario discovery
## =====================================================================================================
  Number.Cores<-18
 #Define directory parameters
  dir.inputs<-paste(root,"RDM Inputs\\",sep="")
  dir.harness<-paste(root,"RDM Harness\\",sep="")
  dir.output<-paste(root,"RDM Outputs\\",sep="")
 #load needed libraries
  library(reshape2)
  library(data.table)
  library(snow)
#create vector with file names
  filenames <- list.files(dir.harness, pattern="*.csv", full.names=FALSE)
#source function to process harnessed output data
  source(paste(dir.inputs,"harness_processing.r",sep=""))
#run post-processing in parallel
  nCore<-Number.Cores
  cl <- makeSOCKcluster(names = rep('localhost',nCore))
  global.elements<-list("dir.harness","process.prim.data","data.table")
  clusterExport(cl,global.elements,envir=environment())
  prim.data <- parLapply(cl,filenames, function(x){data.table(process.prim.data(x,dir.harness))} )
  stopCluster(cl)
  prim.data<-rbindlist(prim.data)
#merge data with experimental design
  experiment.version<-"Exp.design.csv"
  prim.data<-merge.exp.design(dir.inputs,experiment.version,prim.data)
#create future without action consumption
  prim.data.fwa<-subset(prim.data,prim.data$policy.name=="FWA")
  prim.data.fwa<-prim.data.fwa[,c("Future.ID","Y.Total_N","Y.Total_S","Consumption.Total_N","Consumption.Total_S","Climate.Coef","CO2.GrowthRate","Delta.Temp.GrowthRate","Consumption.Total_N.300","Consumption.Total_S.300"),with=FALSE]
  setnames( prim.data.fwa,c("Y.Total_N","Y.Total_S","Consumption.Total_N","Consumption.Total_S","Climate.Coef","CO2.GrowthRate","Delta.Temp.GrowthRate","Consumption.Total_N.300","Consumption.Total_S.300"),c("Y.Total_N.fwa","Y.Total_S.fwa","Consumption.Total_N.fwa","Consumption.Total_S.fwa","Climate.Coef.fwa","CO2.GrowthRate.fwa","Delta.Temp.GrowthRate.fwa","Consumption.Total_N.300.fwa","Consumption.Total_S.300.fwa"))
  prim.data<-merge(prim.data,prim.data.fwa,by="Future.ID")
#standarize relative values
  prim.data$Z.Relative.Gamma<-prim.data$Gamma.re/prim.data$Gamma.ce
  prim.data$Z.Relative.Gamma<-scale(prim.data$Z.Relative.Gamma, center=TRUE, scale=TRUE)
  prim.data$Z.Relative.Eta<-prim.data$Eta.re/prim.data$Eta.ce
  prim.data$Z.Relative.Eta<-scale(prim.data$Z.Relative.Eta, center=TRUE, scale=TRUE)
  prim.data$Z.Relative.Nu<-prim.data$Nu.re/prim.data$Nu.ce
  prim.data$Z.Relative.Nu<-scale(prim.data$Z.Relative.Nu, center=TRUE, scale=TRUE)
  prim.data$Z.epsilon<-scale(prim.data$epsilon, center=TRUE, scale=TRUE)

  #write.csv(prim.data, paste(dir.output, "prim.data_7_06_2015.csv", sep=""), row.names=FALSE)
  write.csv(prim.data, paste(dir.output, "prim.data_extras_seminar.csv", sep=""), row.names=FALSE)

## =====================================================================================================
## This section reads the output of simulations and reshapes it into time series split by region,
## =====================================================================================================
Number.Cores<-14
#Define directory parameters
 dir.inputs<-paste(root,"RDM Inputs\\",sep="")
 dir.harness<-paste(root,"RDM Harness\\",sep="")
 dir.output<-paste(root,"RDM Outputs\\",sep="")

#crate vector with file names
 experiment.version<-"Exp.design.csv"
 #experiment.version<-"Exp.design_with_control_runs.csv"
 #experiment.version<-"Exp.design_defense_seminar_extras.csv"
 filenames <- list.files(dir.harness, pattern="*.csv", full.names=FALSE)
#source function to process harnessed output data
 source(paste(dir.inputs,"harness_processing.r",sep=""))
#run post-processing in parallel
  library(data.table)
  library(snow)
  nCore<-Number.Cores
  cl <- makeSOCKcluster(names = rep('localhost',nCore))
  global.elements<-list("dir.inputs","experiment.version","dir.harness","process.harness.data")
  clusterExport(cl,global.elements,envir=environment())
  modelruns <- parLapply(cl,filenames, function(x){process.harness.data(x,dir.inputs,experiment.version,dir.harness)} )
  stopCluster(cl)
  modelruns<-rbindlist(modelruns)
#print time series for model
  write.csv(modelruns, paste(dir.output, "model.runs_7_09_2015.csv", sep=""), row.names=FALSE)
