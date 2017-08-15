
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
  if (name == "diff") {
    dataIn <- list(halflData = c(49.3, 335.9, 106.9, 204.3, 24), 
                   halflStd = c(2.7, 14.9, 4.3, 5.2, 1.0))
  } else if (name == "scarlette") {
    dataIn <- list(halflData = c(101.1, 323.0, 284.1, 294.7, 24),
                   halflStd = c(11.4, 24.1, 9.7, 17.8, 1.0))
  } else if (name == "marlene") {
    dataIn <- list(halflData = c(0, 0, 0, 0, 24),
                   halflStd = c(0, 0, 0, 0, 1.0))
  }
  
  samples_filename <- paste(name, "samples.rds", sep = "_")
  
  fit <- rstan::stan(system.file("extdata", "model.stan", package = "fcrn"),
                     cores = parallel::detectCores(),
                     data = dataIn,
                     iter = 500,
                     save_dso = T,
                     chains = 4,
                     verbose = T,
                     control = list(adapt_delta = 0.99));
  
  save(fit, file = system.file("extdata", samples_filename, package = "fcrn"))
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