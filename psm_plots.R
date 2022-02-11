source("utils.R")
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(mcclust)
library(mclust)
library(pheatmap)

psm_plots <- function(dataset, datasets, name, focus_dataset=NULL) {
    obsData = read.table(paste0("data/", dataset, "/beta.tsv"),
                         header=1, row.names=1)
    obsVars = read.table(paste0("data/", dataset, "/var.tsv"),
                         header=1, row.names=1)

    result = calc_psms(dataset, datasets)
    bigpsm = result$bigpsm
    psms = result$psms

    print(isSymmetric(bigpsm))
    print(max(bigpsm))
    print(min(bigpsm))
    print(diag(bigpsm))
    rownames(bigpsm) = rownames(obsData)
    colnames(bigpsm) = rownames(obsData)

    calls=maxpear(bigpsm, method="comp") ## calls

    # MClust solution
    BIC <- mclustBIC(obsData)
    mclust_solution <- Mclust(obsData, x=BIC)

    annotations = get_ann_colors(calls$cl, mclust_solution$classification, obsData)
    annotations$ann$Var.PC1 = obsVars[, 1]

    # Same hclust calculation that maxpear does prior to choosing optimal cut point
    hclust.comp <- hclust(as.dist(1 - bigpsm), method = "complete")
    psm_heatmap = pheatmap(bigpsm,
                           show_rownames = TRUE,
                           show_colnames = FALSE,
                           cluster_rows = hclust.comp,
                           cluster_cols = hclust.comp,
                           annotation_names_row = FALSE,
                           treeheight_col=0,
                           fontsize_row=6,
                           annotation_row = annotations$ann,
                           annotation_col = annotations$ann,
                           color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                            name = "Blues")))(100),
                           annotation_colors = annotations$colors,
                           width=20,
                           height=14,
                           filename=paste0("plots/psm_heatmap_", name, ".png"))

#    if (! is.null(focus_dataset) ) {
#        print("Plotting focus dataset")
#        cluster = calls$cl[focus_dataset]
#        print("Dataset is in cluster")
#        print(cluster)
#        rows_in_cluster = rownames(obsData)[calls$cl == cluster]
#        print(rows_in_cluster)
#        psm_heatmap = pheatmap(bigpsm[rows_in_cluster, rows_in_cluster],
#                               show_rownames = FALSE,
#                               show_colnames = TRUE,
#                               annotation_row = annotations$ann,
#                               annotation_col = annotations$ann,
#                               fontsize_col=3,
#                               annotation_legend=FALSE,
#                               color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
#                                                                                name = "Blues")))(100),
#                               annotation_colors = annotations$colors,
#                               filename=paste0("plots/psm_heatmap_focus_", name, ".png"))
#    }

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }

    customColours = generate_balanced_colours(obsData)
    obs_heatmap = pheatmap(obsData,
                           clustering_method="complete",
                           cluster_rows = hclust.comp,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           cluster_col = FALSE,
                           color = customColours$colours,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0("plots/obs_heatmap_", name, ".png"))

    if (!grepl("novar", name)) {
        customColours = generate_balanced_colours(obsVars)
        obs_vars_heatmap = pheatmap(obsVars,
                                    clustering_method="complete",
                                    cluster_rows = hclust.comp,
                                    cluster_col = FALSE,
                                    annotation_colors = annotations$colors,
                                    annotation_row = annotations$ann,
                                    color = customColours$colours,
                                    fontsize_col = 8,
                                    fontsize_row = 6,
                                    width=9,
                                    height = 14,
                                    breaks = customColours$breaks,
                                    filename=paste0("plots/obs_vars_heatmap_", name, ".png"))
    }

    heatmaps = lapply(psms, function(x) pheatmap(x,
                                                 legend=FALSE,
                                                 color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                                                  name = "Blues")))(100),
                                                 cluster_rows = hclust.comp,
                                                 cluster_cols = hclust.comp,
                                                 treeheight_col=0,
                                                 border_color = NA,
                                                 treeheight_row=0))
    heatmap_grid = grid.arrange(grobs=lapply(heatmaps, function(x) x[[4]]))
    ggsave(paste0("plots/psm_heatmap_grid_", name, ".png"), heatmap_grid)

    individual_calls = do.call(cbind, lapply(psms, function(x) adjust_labels_B_to_match_A(calls$cl, maxpear(x, method="comp")$cl)))
    all_calls = data.frame(cbind(calls$cl, individual_calls))
    colnames(all_calls) = c("Overall", paste("Seed", 1:length(datasets)))

    print(all_calls)

    map_to_call_counts <- function(calls) {
        call_labels = paste0(LETTERS[1:max(calls)], " (", table(calls), ")")
        return(plyr::mapvalues(calls, from=1:max(calls), to=call_labels))
    }
    cell_labels = apply(all_calls, MARGIN=2, map_to_call_counts)

    calls_heatmap = pheatmap(all_calls,
                             display_numbers = cell_labels,
                             fontsize_number = 4,
                             number_color = "#CCCCCC",
                             show_rownames = TRUE,
                             show_colnames = TRUE,
                             color = palette[1:max(all_calls)],
                             breaks = seq(0.5, max(all_calls) + 0.5, by=1),
                             cluster_col = FALSE,
                             cluster_row = hclust.comp,
                             fontsize_row = 6,
                             width=8,
                             height=14,
                             filename=paste0("plots/calls_heatmap_", name, ".png"))

    return(list(name=name, bigpsm=bigpsm, calls=calls, hclust.comp=hclust.comp, mclust_calls=mclust_solution$classification))
}

burren_heatmaps <- function(directory) {
  dir.create(paste0(directory, "/plots/"), recursive = TRUE)

  psm_results = list()
  for (dataset in c("with_subtypes_001", "with_subtypes_noGA_001", "with_subtypes_noGA_001_novar", "with_subtypes_001_novar")) {
       datasets = paste0("./trimmed_results/", dataset, "/seed", 1001:1010) 
       result = psm_plots(dataset, datasets, name=dataset)
       psm_results[[dataset]] = result
  }
  save(psm_results, file="psm_results.Rda")
}

burren_heatmaps("./")
