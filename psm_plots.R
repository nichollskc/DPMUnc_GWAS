source("process_v2_functions.R")
library(R.cache)
library(clue)
library(mclust)
library(ggplot2)
library(grid)
library(gridExtra)

adjust_labels_B_to_match_A <- function(calls_A, calls_B) {
    K_A = length(unique(calls_A))
    K_B = length(unique(calls_B))

    if (K_A < max(calls_A)) {
        print("WARNING: assumptions about cluster labels violated")
    }
    if (K_B < max(calls_B)) {
        print("WARNING: assumptions about cluster labels violated")
    }

    jaccard_mat = matrix(0,
                         nrow=max(K_A, K_B),
                         ncol=max(K_A, K_B))
    for (i in 1:K_A) {
        for (j in 1:K_B) {
            in_A = calls_A == i
                in_B = calls_B == j
                jacc = sum(in_A & in_B) / sum(in_A | in_B)
                jaccard_mat[i, j] = jacc
        }
    }

    new_labels_for_B = c(solve_LSAP(t(jaccard_mat), maximum=TRUE))[1:K_B]

    return(plyr::mapvalues(calls_B, from=1:K_B, to=new_labels_for_B))
}

palette <- c("#77AADD", "#000000", "#9E0142", "#D53E4F", "#F46D43", "#FEE08B", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "#FDAE61", "#E6F598", "#771155", "#AA4488", "#CC99BB", "#114477", "#774411", "#EEEEEE", "#117777", "#117744", "#44AA77", "#88CCAA", "#777711", "#44AAAA", "#AAAA44", "#77CCCC", "#DDDD77", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

get_ann_colors=function(calls, mclust_calls, obsData, verbose=TRUE) {
  # from spectral, plus some extras
  #palette= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788") %>%
#    matrix(.,7,3,byrow=TRUE) %>%
#    as.vector()
  counts = data.frame(table(calls))
  cluster_labels = paste0(LETTERS[1:length(counts$Freq)], " (", counts$Freq, ")")

  mclust_calls = adjust_labels_B_to_match_A(calls, mclust_calls)
  mclust_counts = data.frame(table(mclust_calls))
  mclust_cluster_labels = paste0(letters[1:length(mclust_counts$Freq)], " (", mclust_counts$Freq, ")")

  ann=data.frame(row.names=rownames(obsData),
                 cluster=factor(calls, labels=cluster_labels),
                 mclust=factor(mclust_calls, labels=mclust_cluster_labels))
  ncalls=length(unique(calls))
  ncalls_mclust=length(unique(mclust_calls))
  ann_colors=list(cluster = structure(palette[1:ncalls], names=cluster_labels),
                  mclust = structure(palette[1:ncalls_mclust], names=mclust_cluster_labels))
  list(ann=ann,colors=ann_colors)
}

raw_calc_psms <- function(dataset, datasets) {
    allocs=lapply(paste0(datasets, "/clusterAllocations.tsv"), fread)## read the allocations
    # This line is essential for some reason
    allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),]))
    #bigalloc = do.call(rbind, allocs)
    bigalloc=lapply(allocs, function(x) x[-c(1:(nrow(x)/2)),]) %>% do.call("rbind",.) ## combine, discarding first 50%
    bigpsm=calc_psm(bigalloc,burn=0) ## make a psm, don't discard any burn in because already discarded

    psms = lapply(allocs, function(x) calc_psm(x, burn=0))
    return(list(bigpsm=bigpsm, psms=psms))
}

calc_psms <- addMemoization(raw_calc_psms)

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
    psm_heatmap = pheatmap(bigpsm,
                           show_rownames = TRUE,
                           show_colnames = FALSE,
                           clustering_method = "complete",
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
                           cluster_row= psm_heatmap$tree_row,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           color = customColours$colours,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0("plots/obs_heatmap_", name, ".png"))

    customColours = generate_balanced_colours(obsVars)
    obs_vars_heatmap = pheatmap(obsVars,
                           clustering_method="complete",
                           cluster_row= psm_heatmap$tree_row,
                           annotation_colors = annotations$colors,
                           annotation_row = annotations$ann,
                           color = customColours$colours,
                           fontsize_col = 8,
                           fontsize_row = 6,
                           width=9,
                           height = 14,
                           breaks = customColours$breaks,
                           filename=paste0("plots/obs_vars_heatmap_", name, ".png"))

    heatmaps = lapply(psms, function(x) pheatmap(x,
                                                 legend=FALSE,
                                                 color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                                                  name = "Blues")))(100),
                                                 cluster_col = psm_heatmap$tree_row,
                                                 cluster_row = psm_heatmap$tree_row,
                                                 treeheight_col=0,
                                                 border_color = NULL,
                                                 treeheight_row=0))
    heatmap_grid = grid.arrange(grobs=lapply(heatmaps, function(x) x[[4]]))
    ggsave("plots/psm_heatmap_grid.png", heatmap_grid)

    individual_calls = do.call(cbind, lapply(psms, function(x) adjust_labels_B_to_match_A(calls$cl, maxpear(x)$cl)))
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
                             cluster_row = psm_heatmap$tree_row,
                             fontsize_row = 6,
                             width=8,
                             height=14,
                             filename="plots/calls_heatmap.png")
}

burren_heatmaps <- function(directory) {
  dir.create(paste0(directory, "/plots/"), recursive = TRUE)

  dataset = "with_subtypes_001"
  datasets = paste0("./trimmed_results/", dataset, "/seed", 1001:1010) 
  psm_plots(dataset, datasets, name=dataset)

}

burren_heatmaps("./")
