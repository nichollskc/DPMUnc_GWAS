source("process_v2_functions.R")
library(mclust)

psm_plots <- function(dataset, datasets, name, focus_dataset=NULL) {
    obsData = read.table(paste0("data/", dataset, "/beta.tsv"),
                         header=1, row.names=1)

    allocs=lapply(paste0(datasets, "/clusterAllocations.tsv"), fread)## read the allocations
    # This line is essential for some reason
    allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),]))
    #bigalloc = do.call(rbind, allocs)
    bigalloc=lapply(allocs, function(x) x[-c(1:(nrow(x)/2)),]) %>% do.call("rbind",.) ## combine, discarding first 50%
    bigpsm=calc_psm(bigalloc,burn=0) ## make a psm, don't discard any burn in because already discarded
    print(isSymmetric(bigpsm))
    print(max(bigpsm))
    print(min(bigpsm))
    print(diag(bigpsm))
    rownames(bigpsm) = rownames(obsData)
    colnames(bigpsm) = rownames(obsData)

    calls=maxpear(bigpsm) ## calls

    # MClust solution
    BIC <- mclustBIC(obsData)
    mclust_solution <- Mclust(obsData, x=BIC)

    annotations = get_ann_colors(calls$cl,obsData)
    annotations$ann$mclust = mclust_solution$classification
    annotations$colors$cluster = annotations$colors$call
    annotations$colors$mclust = annotations$colors$call[seq(1, 11, by=2)]
    annotations$myositis = NULL
    psm_heatmap = pheatmap(bigpsm,
                           show_rownames = FALSE,
                           show_colnames = TRUE,
                           annotation_row = annotations$ann,
                           annotation_col = annotations$ann,
                           fontsize_col=6,
                           annotation_legend=FALSE,
                           color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                            name = "Blues")))(100),
                           annotation_colors = annotations$colors,
                           width=14,
                           height=16,
                           filename=paste0("plots/psm_heatmap_", name, ".png"))

    if (! is.null(focus_dataset) ) {
        print("Plotting focus dataset")
        cluster = calls$cl[focus_dataset]
        print("Dataset is in cluster")
        print(cluster)
        rows_in_cluster = rownames(obsData)[calls$cl == cluster]
        print(rows_in_cluster)
        psm_heatmap = pheatmap(bigpsm[rows_in_cluster, rows_in_cluster],
                               show_rownames = FALSE,
                               show_colnames = TRUE,
                               annotation_row = annotations$ann,
                               annotation_col = annotations$ann,
                               fontsize_col=3,
                               annotation_legend=FALSE,
                               color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                                name = "Blues")))(100),
                               annotation_colors = annotations$colors,
                               filename=paste0("plots/psm_heatmap_focus_", name, ".png"))
    }

    generate_balanced_colours <- function(obj) {
        paletteLength = 100
        customColours <- colorRampPalette(c("firebrick3", "white", "navy"))(paletteLength)
        customBreaks <- c(seq(min(obj), 0, length.out=ceiling(paletteLength/2) + 1),
                          seq(max(obj)/paletteLength, max(obj), length.out=floor(paletteLength/2)))
        return(list("colours"=customColours, "breaks"=customBreaks))
    }
    customColours = generate_balanced_colours(obsData)

    obs_heatmap = pheatmap(t(obsData),
                           clustering_method="average",
                           cluster_col = psm_heatmap$tree_row,
                           annotation_legend = FALSE,
                           annotation_colors = annotations$colors,
                           annotation_col = annotations$ann,
                           color = customColours$colours,
                           fontsize_col = 3,
                           breaks = customColours$breaks,
                           filename=paste0("plots/obs_heatmap_", name, ".png"))
}

burren_heatmaps <- function(directory) {
  dir.create(paste0(directory, "/plots/"), recursive = TRUE)

  dataset = "with_subtypes_001"
  datasets = paste0("./trimmed_results/", dataset, "/seed", 1:5) 
  psm_plots(dataset, datasets, name=dataset)

}

burren_heatmaps("./")
