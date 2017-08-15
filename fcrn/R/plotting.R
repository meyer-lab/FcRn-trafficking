globalVariables("%>%")

#' Title
#'
#' @param fit The rstan fit object with all the samples.
#'
#' @return The ggplot2 object.
#' @export
plot_halfls <- function(name) {
  best_fit <- getML(loadsample(name))
  
  output <- expand.grid(sortF = seq(0.5, 0.98, length.out = 50),
                        releaseF = seq(0.25, 1, length.out = 50))
  
  output$Ve <- best_fit$Ve
  
  # Run the predictions
  output$halfl <- apply(output, 1, fcrn::halfl_fcrn)
  
  e <- ggplot2::ggplot(output, ggplot2::aes_(x = 'sortF', y = 'releaseF', z = 'halfl')) + 
    ggplot2::geom_contour(ggplot2::aes_(color = '..level..'), breaks = c(48, 72, 96, 200, 300, 400, 500, 700)) +
    ggplot2::annotate("text", label = "WT", x = best_fit$actual_sortF_wt, y = 0.99, color = "black") +
    ggplot2::annotate("text", label = "DHS", x = best_fit$actual_sortF_dhs, y = 0.99, color = "black") +
    ggplot2::annotate("text", label = "LS", x = best_fit$actual_sortF_ls,
                      y = best_fit$actual_release_ls, color = "black") +
    ggplot2::annotate("text", label = "YTE", x = best_fit$sortF_yte,
                      y = best_fit$releaseF_yte, color = "black") +
    ggplot2::xlab("Endosomal Sorting Fraction") +
    ggplot2::ylab("Fraction Released at Surface") +
    ggplot2::scale_x_continuous(trans = "atanh",
                                limits = c(0.5, 0.98),
                                breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98)) +
    ggplot2::scale_colour_gradient(low = "gray", high = "gray")
  
  e <- directlabels::direct.label(e, "top.pieces")
  
  return(e)
}


#' Title
#'
#' @param fit Stan fit chain.
#'
#' @return Maximum likelihood point from sampling.
getML <- function(fit) {
  best_fit <- as.data.frame(rstan::extract(fit)) %>%
    dplyr::filter_('lp__' == max('lp__'))
}