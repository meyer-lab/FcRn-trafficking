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
  
  output$Vin <- best_fit$Vin
  output$Q <- best_fit$Q
  output$Qu <- best_fit$Qu
  output$Vp <- best_fit$Vp
  
  # Run the predictions
  output$halfl <- apply(output, 1, fcrn::halfl_fcrn)
  
  e <- ggplot2::ggplot(output, ggplot2::aes_(x = ~sortF, y = ~releaseF, z = ~halfl)) + 
    ggplot2::geom_contour(ggplot2::aes_(color = ~..level..), breaks = c(48, 72, 96, 200, 300, 400, 500, 700)) +
    ggplot2::annotate("text", label = "WT", x = best_fit$actual_sortF_wt, y = 0.99, color = "black", size = 2) +
    ggplot2::annotate("text", label = "DHS", x = best_fit$actual_sortF_dhs, y = 0.99, color = "black", size = 2) +
    ggplot2::annotate("text", label = "LS", x = best_fit$actual_sortF_ls,
                      y = best_fit$actual_release_ls, color = "black", size = 2) +
    ggplot2::annotate("text", label = "YTE", x = best_fit$sortF_yte,
                      y = best_fit$releaseF_yte, color = "black", size = 2) +
    ggplot2::xlab("Endosomal Sorting Fraction") +
    ggplot2::ylab("Fraction Released at Surface") +
    ggplot2::scale_x_continuous(trans = "atanh",
                                limits = c(0.5, 0.98),
                                breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.97, 0.98)) +
    ggplot2::scale_colour_gradient(low = "gray", high = "gray") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  e <- directlabels::direct.label(e, list("top.pieces", cex=0.5))
  
  return(e)
}

#' Title
#'
#' @return ggplot object of plot.
#' @export
plot_otherPs <- function() {
  globs <- dplyr::mutate_(loadAll(), 'param = as.factor(param)') %>%
    dplyr::filter(param == "Q" | param == "Qu" | param == "Vin" | param == "Vp")
  
  gg <- ggplot2::ggplot(globs, ggplot2::aes_(x = ~`param`, y = ~`50%`, color = ~`model`)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(.9)) +
    ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~`25%`, ymax = ~`75%`),
                           position = ggplot2::position_dodge(.9),
                           width = 0.5) +
    ggplot2::theme(legend.position = c(0.5, 0.9),
                   legend.title=ggplot2::element_blank()) +
    ggplot2::ylab("Parameter Value") +
    ggplot2::xlab("")
  
  return(gg)
}

#' Title
#'
#' @param fit Stan fit chain.
#'
#' @return Maximum likelihood point from sampling.
getML <- function(fit) {
  best_fit <- as.data.frame(rstan::extract(fit)) %>%
    dplyr::filter_('lp__ == max(lp__)')
  
  return(best_fit)
}


#' Produce a plot showing the posteriors for each sorting function across the in vivo models.
#'
#' @return Plot object
#' @export
getSortingPosterior <- function() {
  # Load all the sampling data
  all <- loadAll()
  
  all <- dplyr::select_(all, '-mean', '-se_mean', '-sd', '-n_eff', '-Rhat') %>%
    dplyr::mutate(param = as.factor(param)) %>%
    dplyr::filter(param != "lp__") %>%
    dplyr::filter(param != "sortF_wt") %>%
    dplyr::filter(param != "sortF_dhs") %>%
    dplyr::filter(param != "sortF_ls") %>%
    dplyr::filter(param != "releaseF_ls")
  
  sorts <- all %>%
    dplyr::filter(!(param == "Q" | param == "Qu" | param == "Vin" | param == "Vp")) %>%
    dplyr::mutate(param = gsub("actual_", "", param)) %>%
    dplyr::mutate(param = gsub("release_", "releaseF_", param)) %>%
    tidyr::separate(param, into = c("params", "IgG"), sep = "_") %>%
    reshape2::melt(id.vars = c("params", "IgG", "model"), variable.name = "quantile") %>%
    reshape2::dcast(IgG + model ~ params + quantile, value.var = "value")
  
  sorts[is.na(sorts)] <- 1.0

  print(sorts)
  
  g <- ggplot2::ggplot(sorts, ggplot2::aes_(x = ~`sortF_50%`,
                                            y = ~`releaseF_50%`,
                                            color = ~`IgG`)) +
    ggplot2::geom_point() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = c(0, 1),
                   legend.title=ggplot2::element_blank()) +
    ggplot2::facet_wrap(~model) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::geom_errorbar(ggplot2::aes_(ymin = ~`releaseF_2.5%`, ymax = ~`releaseF_97.5%`)) +
    ggplot2::geom_errorbarh(ggplot2::aes_(xmin = ~`sortF_2.5%`, xmax = ~`sortF_97.5%`)) +
    ggplot2::xlab("Endosomal Sorting Fraction") +
    ggplot2::ylab("Surface Release Fraction")
  
  return(g)
}



#' Title
#'
#' @param save Boolean for whether this should save a PDF of the figure.
#'
#' @return Object with all the plots assembled.
#' @export
full_plot <- function(save = T) {
  ggplot2::theme_set(cowplot::theme_cowplot(font_size = 8))
  
  gg <- cowplot::ggdraw() +
    cowplot::draw_plot(plot_otherPs(), 0.37, 0.5, 0.3, 0.45) +
    cowplot::draw_plot(getSortingPosterior(), 0.7, 0.5, 0.3, 0.45) +
    cowplot::draw_plot(plot_halfls("diff") + ggplot2::ggtitle("hemizygous Tg276"), 0.04, 0.0, 0.3, 0.45) +
    cowplot::draw_plot(plot_halfls("marlene") + ggplot2::ggtitle("Marlene"), 0.37, 0.0, 0.3, 0.45) +
    cowplot::draw_plot(plot_halfls("scarlette") + ggplot2::ggtitle("Scarlette"), 0.7, 0.0, 0.3, 0.45) +
    cowplot::draw_plot_label(label = c("A", "B", "C", "D", "E", "F"),
                             size = 15,
                             x = c(0, 0.33, 0.66, 0, 0.33, 0.66),
                             y = c(1, 1, 1, 0.5, 0.5, 0.5))
  
  if (save) {
    cowplot::save_plot(filename = "output.pdf",
                       plot = gg,
                       ncol = 1,
                       base_width = 7,
                       base_height = 4)
  }
  
  return(gg)
}