library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

read_numeric=function(dirs, filename) {
  files=file.path(dirs,filename)
  use=file.exists(files)
  values=lapply(which(use), function(i) {
    message(dirs[i])
    scan(files[i],what="") %>% as.numeric()
  })
  nmin=sapply(values,length) %>% min()
  values %<>% lapply(., "[", 1:nmin) %>% do.call("cbind",.)
  colnames(values)=dirs[use]
  values
}

get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

quantile_traceplot <- function(datasets, label, filename, block_size=100) {
  data_mat = read_numeric(datasets, filename)
  df = data.frame(data_mat)
  colnames(df) = paste0("seed", 1:ncol(df))
  df["block"] = rep(1:(nrow(df) / block_size), each=block_size)
  df = df%>% tidyr::pivot_longer(cols=1:(ncol(df) - 1))
  
  quantiled = df %>%
    group_by(block, name) %>%
    summarise(value = quantile(value, c(0.1, 0.25, 0.5, 0.75, 0.9)),
              quantile = c(0.1, 0.25, 0.5, 0.75, 0.9))
  
  g = ggplot(quantiled, aes(x=block, y=value, colour=name)) +
    geom_vline(xintercept=max(df["block"])/ 2, colour="red") + 
    scale_color_manual(values=cbbPalette) +
    geom_line(data=subset(quantiled,quantile==0.5)) +
    geom_line(data=subset(quantiled,quantile==0.1), alpha = 0.2, linetype='dashed') +
    geom_line(data=subset(quantiled,quantile==0.9), alpha = 0.2, linetype='dashed') +
    geom_line(data=subset(quantiled,quantile==0.25), alpha = 0.2) +
    geom_line(data=subset(quantiled,quantile==0.75), alpha = 0.2) +
    labs(y=label)
  return(g)
}

quantile_traceplots_dataset <- function(dataset, datasets, block_size=100) {
  g_alpha = quantile_traceplot(datasets, "alpha", "alpha.tsv")
  g_K = quantile_traceplot(datasets, "K", "K.txt") + theme(legend.position = "none")
  g_latent = quantile_traceplot(datasets, "pLatents", "pLatentsGivenClusters.tsv")
  g = grid.arrange(g_alpha + theme(legend.position = "none"),
                   g_K + theme(legend.position = "none"),
                   g_latent + theme(legend.position = "none"),
                   get_only_legend(g_latent),
                   top = textGrob(paste0("Quantile traceplot - ", dataset),gp=gpar(fontsize=20,font=3)),
                   ncol=4,
                   widths=c(2,2,2,1))
  ggsave(paste0("./plots/trace_quantiles_", dataset, ".png"), g, width=12, height=7, units="in")
  return(g)
}

burren_traceplots <- function(directory) {
  dir.create(paste0(directory, "/plots/"), recursive = TRUE)
  datasets = paste0(directory, "/trimmed_results/with_subtypes_001/seed", 1:5, "/")
  quantile_traceplots_dataset("with_subtypes_001", datasets)
}

burren_traceplots("./")
