
# function to plot KM from survival fit
plot_KM <- function(fit, gene, colors, labels = NULL, lty = NULL) {
  plotdat <- data.frame(time = fit$time,
                        surv = fit$surv,
                        group = rep(names(fit$strata), fit$strata))
  if (is.null(labels)) {
    labels <- paste0(gene, " ", sub("^.+\\=(.+)$", "\\1", names(fit$strata)),
                     " (n=", fit$n, ")")
    plotdat$group <- factor(plotdat$group, levels = names(fit$strata),
                            labels = labels)
  } else {
    plotdat$group <- factor(plotdat$group, levels = names(fit$strata),
                            labels = labels)
  }
  plotdat <- plotdat[rep(seq_len(nrow(plotdat)), each = 2), ]
  plotdat$surv <- c(1, plotdat$surv[1:(length(plotdat$surv)-1)])
  switchidx <- which(diff(plotdat$time)<0)+1
  plotdat[c(1,switchidx),]$surv <- 1
  tmp <- plotdat[c(1,switchidx),]
  tmp$time <- 0
  plotdat <- rbind(plotdat, tmp)
  maxx <- max(plotdat$time)
  xpos <- maxx - maxx/6
  pval <- data.frame(x=xpos, y=0.85, label=fit$pval_txt)
  if (is.null(lty)) lty <- rep(1, length(levels(plotdat$group)))
  ggplot(plotdat, aes(x = time, y = surv)) +
    geom_line(aes(col = group, lty = group), lwd = 1.1, alpha = .6) +
    theme_classic(base_size = 14, base_family = "Arial") +
    geom_text(data=pval, aes(x=x, y=y, label=label),
              size = 5, family="Arial") +
    scale_color_manual(values = colors, name=NULL) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0), 
                       name = "OS probability") +
    scale_x_continuous(expand = c(0,0), name = "Time [years]") +
    scale_linetype_manual(values=lty, name=NULL) +
    theme(strip.background = element_blank(),
          panel.spacing = unit(1, "line"),
          plot.title = element_text(size=14),
          legend.position="bottom")
}