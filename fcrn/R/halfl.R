
#' Builds the jacobian of the ODE model.
#'
#' @param p List of parameters with which to solve.
#'
#' @return Jacobian of the ODE system.
jacFunc <- function(p) {
  with(as.list(c(p)), {
    matrix(c(-Q,  Q/Vp,          0,
             Q,  -(Q + Qu)/Vp,  Qu/Vin,
             0,  Qu*sortF*releaseF/Vp,  Qu/Vin*((1 - releaseF)*sortF - 1)
    ), 3, 3)
  })
}

#' Solves for the model half-life given a list of paraemeters.
#'
#' @param th List of parameters with which to solve.
#'
#' @return The half-life given this parameter set.
#' @export
halfl_fcrn <- function(th) {
  C0 <- 20.0
  
  jac <- jacFunc(th)
  
  outt <- stats::uniroot(function(t) expm::expAtv(jac, c(C0, 0.0, 0.0), t)$eAtv[1] - C0/2,
                         tol = 1E-1,
                         lower = 0,
                         upper = 10000)
  
  return(outt$root)
}

models <- c("diff", "scarlette", "marlene")

#' Title
#'
#' @param name Name of the in vivo model to use.
#' @export
#'
#' @return Nothing.
#' @export
runsample <- function(name) {
  if (name == "diff") {         # WT, DHS, LS, YTE
    dataIn <- list(halflData = c(49.3, 335.9, 106.9, 204.3), 
                   halflStd = c(2.7, 14.9, 4.3, 5.2))
  } else if (name == "scarlette") {
    dataIn <- list(halflData = c(101.1, 323.0, 284.1, 294.7),
                   halflStd = c(11.4, 24.1, 9.7, 17.8))
  } else if (name == "marlene") {
    dataIn <- list(halflData = c(148.9, 417.9, 412.3, 359.9),
                   halflStd = c(14.3, 24.0, 22.6, 30.5))
  }
  
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  fit <- rstan::stan(system.file("extdata", "model.stan", package = "fcrn"),
                     cores = parallel::detectCores(),
                     data = dataIn,
                     iter = 400,
                     save_dso = T,
                     chains = 4,
                     control = list(adapt_delta = 0.99,
                                    max_treedepth = 20));
  
  print(paste("Writing to: ", samples_filename, sep = ""))
  
  save(fit, file = samples_filename)
}

#' Run the sampling for all the in vivo models.
#'
#' @return Nothing.
#' @export
sampleAll <- function() {
  runsample("diff")
  
  runsample("scarlette")
  
  runsample("marlene")
}

#' Title
#'
#' @param name The name of the model to load.
#'
#' @return None
#' @export
shinyModel <- function(name) {
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  samples_path <- system.file("extdata", samples_filename, package = "fcrn")
  
  load(samples_path)
  
  shinystan::launch_shinystan(fit)
}


#' Title
#'
#' @param name Name of the in vivo model fit to load.
#' @param as.df Do you want the return to be a data.frame?
#'
#' @return Loaded samples from a fit
#' @export
loadsample <- function(name, as.df = F) {
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  samples_path <- system.file("extdata", samples_filename, package = "fcrn")
  
  load(samples_path)
  
  if (as.df) {
    return(as.data.frame(rstan::extract(fit)))
  } else {
    return(fit)
  }
}




#' Title
#'
#' @return Return all the sampling across models
#' @export
loadAll <- function() {
  loadData <- function(name) {
    data <- as.data.frame(rstan::summary(loadsample(name))$summary)
    data$param <- row.names(data)
    data$model <- factor(name)
    return(data)
  }
  
  all <- rbind(loadData("scarlette"), loadData("diff"), loadData("marlene"))
  
  assertthat::assert_that(max(all$Rhat) < 1.1)
  assertthat::assert_that(min(all$n_eff) > 100)
  
  return(all)
}

