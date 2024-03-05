# ES Function ---------------------------------------------------------------
#' Calculate running enrichment score
#'
#' \code{get_es} returns the running enrichment scores for an ordered set of
#' input/test features compared against a refence set of features. E.g. can be
#' used to compare a list of genes ordered by log2 fold change between two
#' conditions of interest against a set of genes involved in a specific pathway.
#'
#' @param input_features named & ordered list of input feature scores
#' @param ref_features vector of feature IDs to compare input features against
#' @return returns a dataframe of running enrichment score accross all input features
#' @export
get_es <- function(input_features, ref_features) {

  # get lengths
  nref <- length(ref_features)
  ninp <- length(input_features)

  # calc running ES
  hits <- misses <- numeric(ninp)
  hit_idx <- names(input_features) %in% ref_features
  hits[hit_idx] <- abs(input_features[hit_idx])^1
  hits <- cumsum(hits / sum(hits))
  misses[!hit_idx] <-  1 / (ninp - nref)
  misses <- cumsum(misses)
  running_es <- hits - misses

  # return df
  df <- data.frame(row.names = NULL,
                   idx = seq_along(running_es),
                   feature = names(input_features),
                   feature_score = input_features,
                   running_score = running_es,
                   in_set = names(input_features) %in% ref_features)
  return(df)
}

# Plotting Function ---------------------------------------------------------------
#' Plot fgsea results
#' \code{plot_fgsea} provides formatting plotting method for results from fgsea package
#' @param fgsea_res an fgsea results object
#' @param feature_list ordered, named list of features used as input to fgsea
#' @param refsets reference feature sets used as input to fgsea
#' @param set_name name of set to plot - must be present in fgsea results
#' @param ylab character string - used for y-axis title
#' @return returns a ggplot2 object.
#' @export
plot_fgsea <- function(fgsea_res, feature_list, refsets, setname,
                       barwidth = NA, linecol = "#bb333399", ylab = "Fold change [log2]",
                       ylim = NA, scale_factor = NA) {
  if (!setname %in% fgsea_res$pathway) {
    stop("Given set name is not present in the FSEA object")
  }
  set <- refsets[[setname]]
  plotdat <- get_es(input_features = feature_list,
                    ref_features = set)
  ymin <- ifelse(is.na(ylim[1]), max(plotdat$feature_score), ylim[[1]])
  ymax <- ifelse(is.na(ylim[2]), min(plotdat$feature_score), ylim[[2]])
  plotdat <- rbind(data.frame(idx = 0, feature = NA, feature_score = ymax,
                              running_score = 0, in_set = FALSE),
                   plotdat[plotdat$feature %in% set, ],
                   data.frame(idx = max(plotdat$idx)+1, feature = NA, feature_score = ymin,
                              running_score = 0, in_set = FALSE))
  scale_factor <- ifelse(!is.na(scale_factor), scale_factor,
                         (max(abs(plotdat$feature_score))) / (max(abs(plotdat$running_score))*2))
  set_res <- fgsea_res[fgsea_res$pathway == setname, , drop = F]
  y_range <- ifelse(set_res$NES > 0, -diff(range(plotdat$feature_score)),
                    diff(range(plotdat$feature_score)))
  x_range <- ifelse(set_res$NES > 0, -diff(range(plotdat$idx)),
                    diff(range(plotdat$idx)))
  y_lim <- ifelse(set_res$NES > 0, max(plotdat$feature_score),
                  min(plotdat$feature_score))
  x_lim <- ifelse(set_res$NES > 0, max(plotdat$idx),
                  min(plotdat$idx))
  # nes_text <- paste0("italic('NES =')~", round(set_res$NES, 2))
  # sig_text <- paste0("italic('q =')~", round(set_res$padj, 3))
  nes_text <- paste0("NES = ", round(set_res$NES, 2))
  sig_text <- paste0("q = ", round(set_res$padj, 3))
  stats <- data.frame(x = rep((x_lim + (x_range/100)), 2),
                      y = c((y_lim + (y_range/5)),
                            (y_lim + (y_range/19))),
                      text = c(nes_text,
                               sig_text))
  if (set_res$NES > 0) stats$y <- rev(stats$y)
  hjust <- ifelse(set_res$NES > 0, 1, 0)
  if (is.na(barwidth)) barwidth <- 5000/sum(plotdat$feature_score != 0)
  p <- ggplot(plotdat, aes(x = idx)) +
    geom_hline(yintercept = 0, alpha = .5) +
    geom_bar(data = dplyr::filter(plotdat, in_set == F),
             aes(y = feature_score), 
             width  = barwidth,
             stat = "identity", alpha = .1,
             color = "white", fill = "white") +
    geom_bar(data = dplyr::filter(plotdat, in_set == T),
             aes(y = feature_score),
             stat = "identity", alpha = .5,
             width = barwidth) +
    geom_line(aes(y = running_score * scale_factor),
              col = linecol, size = 1.5) +
    geom_text(data = stats, aes(x = x, y = y, label = text),
              size = 4, hjust = hjust, parse = F) +
    scale_y_continuous(sec.axis = sec_axis(~. / scale_factor,
                                           name = "Enrichment score"),
                       name = ylab) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black"),
          axis.ticks.y = element_line(color = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y.left = element_text(color = "black"),
          axis.title.y.left = element_text(size = 11, color = "black"),
          axis.text.y.right = element_text(color = linecol),
          axis.title.y.right = element_text(size = 11, color = linecol))
  withCallingHandlers({
    return(p)
  }, warning=function(w) {
    if (startsWith(conditionMessage(w), "position_stack"))
      invokeRestart("muffleWarning")
  })
}

