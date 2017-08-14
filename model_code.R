jacFunc <- function(t, y, p) {
  with(as.list(c(p)), {
    matrix(c(-Q,  Q/Vp,          0,
             Q,  -(Q + Qu)/Vp,  Qu/Vin,
             0,  Qu*sortF*releaseF/Vp,  Qu/Vin*((1 - releaseF)*sortF - 1)
    ), 3, 3)
  })
}

halfl_fcrn <- function(th) {
  C0 <- 20.0
  
  jac <- jacFunc(0.0, 0.0, th)
  
  outt <- uniroot(function(t) expAtv(jac, c(C0, 0.0, 0.0), t)$eAtv[1] - C0/2,
                  lower = 0,
                  upper = 10000)
  
  return(outt$root)
}

plot_halfls <- function(fit) {
  best_fit <- rstan::extract(fit) %>%
    as.data.frame %>%
    filter(lp__ == max(lp__))
  
  sorts <- expand.grid(id = 1,
                       sortF = seq(0.5, 0.98, length.out = 50),
                       releaseF = seq(0.25, 1, length.out = 50))
  
  output <- best_fit %>%
    mutate(id = 1) %>%
    full_join(sorts, by = c("id"))
  
  # Run the predictions
  output$halfl <- apply(output, 1, halfl_fcrn)
  
  e <- ggplot(output, aes(x = sortF, y = releaseF, z = halfl)) + 
    geom_contour(aes(color = ..level..), breaks = c(48, 72, 96, 200, 300, 400, 500, 700)) +
    theme_bw() +
    annotate("text", label = "WT", x = best_fit$actual_sortF_wt, y = 0.99, color = "black") +
    annotate("text", label = "DHS", x = best_fit$actual_sortF_dhs, y = 0.99, color = "black") +
    annotate("text", label = "LS", x = best_fit$actual_sortF_ls,
             y = best_fit$actual_release_ls, color = "black") +
    annotate("text", label = "YTE", x = best_fit$sortF_yte,
             y = best_fit$releaseF_yte, color = "black") +
    xlab("Endosomal Sorting Fraction") +
    ylab("Fraction Released at Surface") +
    scale_x_continuous(trans = "atanh",
                       limits = c(0.5, 0.98),
                       breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98)) +
    scale_colour_gradient(low = "gray", high = "gray")
  
  e <- directlabels::direct.label(e, "top.pieces")
  
  return(e)
}
