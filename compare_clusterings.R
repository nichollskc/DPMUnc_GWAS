library(mclust)
library(pheatmap)
library(ggplot2)
load("psm_results.Rda")
source("utils.R")

base_calls = psm_results$with_subtypes_noGA_001$calls$cl
is_novar_dataset = sapply(psm_results, function(x) grepl("novar", x$name))
print(is_novar_dataset)


noGA_calls = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$calls$cl[names(base_calls)])))
noGA_mclustcalls = do.call(cbind, lapply(psm_results[!is_novar_dataset], function(x) adjust_labels_B_to_match_A(base_calls, x$mclust_calls[names(base_calls)])))
colnames(noGA_calls) = lapply(psm_results, function(x) x$name)
colnames(noGA_mclustcalls) = lapply(psm_results[!is_novar_dataset], function(x) paste0("mclust_", x$name))
combined_calls = data.frame(cbind(noGA_calls, noGA_mclustcalls))
print(head(combined_calls))
plot_ari_for_datasets(combined_calls, "noGA")

map_to_call_counts <- function(calls) {
    call_labels = paste0(LETTERS[1:max(calls)], " (", table(calls), ")")
    return(plyr::mapvalues(calls, from=1:max(calls), to=call_labels))
}

cell_labels = apply(combined_calls, MARGIN=2, map_to_call_counts)
obsVars = read.table(paste0("data/", "with_subtypes_noGA_001", "/var.tsv"), header=1, row.names=1)
ann = data.frame("var.delta.PC1" = obsVars[, 1])
rownames(ann) = rownames(obsVars)
calls_heatmap = pheatmap(combined_calls,
                         display_numbers = cell_labels,
                         fontsize_number = 4,
                         number_color = "#CCCCCC",
                         show_rownames = TRUE,
                         annotation_row = ann,
                         show_colnames = TRUE,
                         color = palette[1:max(combined_calls)],
                         breaks = seq(0.5, max(combined_calls) + 0.5, by=1),
                         cluster_col = FALSE,
                         cluster_row = psm_results$with_subtypes_noGA_001$hclust.comp,
                         fontsize_row = 6,
                         width=8,
                         height=14,
                         filename=paste0("plots/all_calls_with_subtypes_noGA.png"))

base_calls = psm_results$with_subtypes_001$calls$cl
psm_results_withGA = psm_results[sapply(psm_results, function(x) !grepl("noGA", x$name))]
is_novar_dataset = sapply(psm_results_withGA, function(x) grepl("novar", x$name))
print(is_novar_dataset)

calls = do.call(cbind, lapply(psm_results_withGA, function(x) adjust_labels_B_to_match_A(base_calls, x$calls$cl)))
mclustcalls = do.call(cbind, lapply(psm_results_withGA[!is_novar_dataset], function(x) adjust_labels_B_to_match_A(base_calls, x$mclust_calls)))
colnames(calls) = lapply(psm_results_withGA, function(x) x$name)
colnames(mclustcalls) = lapply(psm_results_withGA[!is_novar_dataset], function(x) paste0("mclust_", x$name))
combined_calls = data.frame(cbind(calls, mclustcalls))
print(head(combined_calls))
plot_ari_for_datasets(combined_calls, "withGA")

cell_labels = apply(combined_calls, MARGIN=2, map_to_call_counts)
obsVars = read.table(paste0("data/", "with_subtypes_001", "/var.tsv"), header=1, row.names=1)
ann = data.frame("var.delta.PC1" = obsVars[, 1])
rownames(ann) = rownames(obsVars)
calls_heatmap = pheatmap(combined_calls,
                         display_numbers = cell_labels,
                         fontsize_number = 4,
                         number_color = "#CCCCCC",
                         show_rownames = TRUE,
                         annotation_row = ann,
                         show_colnames = TRUE,
                         color = palette[1:max(combined_calls)],
                         breaks = seq(0.5, max(combined_calls) + 0.5, by=1),
                         cluster_col = FALSE,
                         cluster_row = psm_results$with_subtypes_001$hclust.comp,
                         fontsize_row = 6,
                         width=8,
                         height=14,
                         filename=paste0("plots/all_calls_with_subtypes.png"))
