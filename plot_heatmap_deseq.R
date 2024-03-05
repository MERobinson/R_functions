library(tidyverse)
library(DESeq2)

# plotting function
plot_deseq2_hm <- function(res, dds, lfc=NULL, n=100, row.names=T,
                           colors=c("steelblue4","white","firebrick4")) {
  resdf <- as.data.frame(res) %>% 
    mutate(gene = rownames(res))
  if (is.null(n) | is.na(n)) stop("n required - plotting all genes is not recommnended")
  if (!is.null(lfc)) {
    resdf <- resdf %>% dplyr::filter(abs(log2FoldChange) > lfc)
  }
  resdf <- resdf %>%
    mutate(dir=ifelse(log2FoldChange < 0, "neg", "pos")) %>%
    dplyr::filter(!is.na(dir) & !is.na(padj)) %>%
    group_by(dir) %>%
    arrange(-abs(log2FoldChange)) %>%
    slice_head(n=n) %>%
    ungroup() %>%
    arrange(log2FoldChange)
  rlogdds <- rlog(dds)
  plotdat <- assay(rlogdds)[match(resdf$gene, rownames(rlogdds)),] %>%
    as.data.frame() %>%
    cbind(resdf, .) %>%
    mutate(gene=factor(gene, levels=unique(gene))) %>%
    gather("sample", "expr", -c(1:8)) %>%
    group_by(gene) %>%
    mutate(expr = scale(as.numeric(expr)))
  message(paste0("Plotting ", length(unique(plotdat$gene)), " genes"))
  p <- ggplot(plotdat, aes(x=sample, y=gene, fill=expr)) +
    geom_tile() +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradient2(low = colors[[1]], mid = colors[[2]], high = colors[[3]],
                         name = "Expression level\n[z-score]") +
    theme(panel.background = element_blank(),
          panel.border = element_rect(color="black", fill=NA),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  if (row.names == F) {
    p <- p + theme(axis.text.y = element_blank())
  }
  return(p)
}

# load example data
dds <- makeExampleDESeqDataSet(m=4)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","B","A"))

# plot
plot_deseq2_hm(res=res, dds=dds, n=250, row.names=F)
plot_deseq2_hm(res=res, dds=dds, n=250, lfc=4, row.names=T)
