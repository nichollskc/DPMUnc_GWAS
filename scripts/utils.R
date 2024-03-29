library(clue)
library(data.table)
library(dplyr)
library(grid)
library(gridExtra)
library(magrittr)
library(mcclust)
library(mclust)
library(R.cache)

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
}

# From process_v2_functions.R
raw_calc_psm=function(x,burn=0.5) {
  n=nrow(x)
  if(burn>0)
    x=x[ (burn*n) : n, , drop=FALSE]
  if(any(is.na(x)))
    x=x[ apply(!is.na(x),1,all), ]
  unq=unique(as.vector(x))
  ## print(unq)
  m=matrix(0,ncol(x),ncol(x))
  for(k in unq) {
    xk=matrix(as.numeric(x==k),nrow(x),ncol(x))
    ## m=m + t(xk) %*% xk
    m=m + crossprod(xk)
  }
  psm=m/nrow(x)
  psm
}

calc_psm <- addMemoization(raw_calc_psm)

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

calc_psms <- function(datasets) {
    allocs=lapply(paste0(datasets, "/clusterAllocations.csv"), fread)## read the allocations
    # This line is essential for some reason
    allocs %<>% lapply(., function(x) as.matrix(x[1:nrow(x),]))
#    trimmed_allocs = lapply(allocs, function(x) x[-c(1:(nrow(x)/2)),])
    bigalloc = do.call(rbind, allocs)
    #bigalloc = trimmed_allocs %>% do.call("rbind",.) ## combine, discarding first 50%
    bigpsm=calc_psm(bigalloc,burn=0) ## make a psm, don't discard any burn in because already discarded

#    psms = lapply(trimmed_allocs, function(x) calc_psm(x, burn=0))
    psms = lapply(allocs, function(x) calc_psm(x, burn=0))
    return(list(bigpsm=bigpsm, psms=psms))
}


raw_calc_mclust <- function(obsData) {
    BIC <- mclustBIC(obsData)
    Mclust(obsData, x=BIC)
}

calc_mclust <- addMemoization(raw_calc_mclust)

plot_ari_for_datasets <- function(combined_calls, plots_dir, width=1400, height=1500) {
    print(colnames(combined_calls))
    all_ari = matrix(0, nrow=ncol(combined_calls), ncol=ncol(combined_calls))
    all_ari = data.frame(all_ari)
    rownames(all_ari) = colnames(combined_calls)
    colnames(all_ari) = colnames(combined_calls)
    for (var1 in colnames(all_ari)) {
        for (var2 in rownames(all_ari)) {
            print(paste(var1, var2))
            all_ari[var2, var1] = adjustedRandIndex(combined_calls[, var1], combined_calls[, var2])
        }
    }
    write.csv(all_ari, paste0(plots_dir, "/ari_all_calls.csv"))

    K = apply(combined_calls, MARGIN=2, FUN=function(x) length(unique(x)))
    pheatmap(combined_calls, file=paste0(plots_dir, "/all_calls_with_mclust.png"))
    K_df = data.frame(K)
    K_df$K = as.character(K)

    dist_from_ari <- function(ari_mat) {
        # Set any negative values to 0 and do 1-ARI to get disssimilarity matrix
        all_pos = apply(ari_mat, function(x) max(0, x), MARGIN=1:2)
        distances = apply(ari_mat, function(x) 1-x, MARGIN=1:2)
        as.dist(distances)
    }

    hclust.comp <- hclust(dist_from_ari(all_ari), method="complete")
    ari_hmp = pheatmap(all_ari,
                       annotation_row=K_df,
                       annotation_colors=list("K"=structure(palette[1:length(unique(K_df$K))],
                                                            names=sort(as.numeric(unique(K_df$K))))),
                       display_numbers=TRUE,
                       color=colorRampPalette((RColorBrewer::brewer.pal(n = 7,
                                                                        name = "Blues")))(100),
                       breaks=seq(0, 1, length.out=101),
                       number_color=ifelse(all_ari[hclust.comp$order, hclust.comp$order] < 0.5, "#000000", "#FFFFFF"),
                       cluster_cols=hclust.comp,
                       cluster_row=hclust.comp,
                       cellwidth=22,
                       cellheight=20,
                       legend=FALSE)
    save_pheatmap_png(ari_hmp, paste0(plots_dir, "/ari_all_calls.png"), width=width, height=height)
}

restrict_to_present_samples <- function(sample_info, expression) {
    # Requires expression to have colnames given by sample_info$neat_sample_id
    rownames(sample_info) = sample_info$neat_sample_id
    sample_info = sample_info[colnames(expression)[-c(1)], ]
    rownames(sample_info) = NULL
    return (sample_info)
}
