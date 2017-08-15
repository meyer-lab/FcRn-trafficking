
#' Title
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

#' Title
#'
#' @param th List of parameters with which to solve.
#'
#' @return The half-life given this parameter set.
#' @export
halfl_fcrn <- function(th) {
  C0 <- 20.0
  
  jac <- jacFunc(th)
  
  outt <- stats::uniroot(function(t) expm::expAtv(jac, c(C0, 0.0, 0.0), t)$eAtv[1] - C0/2,
                         lower = 0,
                         upper = 10000)
  
  return(outt$root)
}

#' Title
#' 
#' TODO: Spread data is stdev — check this is right
#' TODO: Need values for FcRn KO
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

shinyModel <- function(name) {
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  load(samples_filename)
  
  shinystan::launch_shinystan(fit)
}


#' Title
#'
#' @param name 
#'
#' @return
#' @export
#'
#' @examples
loadsample <- function(name) {
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  samples_path <- system.file("extdata", samples_filename, package = "fcrn")
  
  load(samples_path)
  
  return(fit)
}