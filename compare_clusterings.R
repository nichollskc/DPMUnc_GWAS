library(mclust)
library(pheatmap)
library(ggplot2)
load("psm_results.Rda")
source("utils.R")

base_calls = psm_results$with_subtypes_noGA_001$calls$cl
noGA_calls = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$calls$cl[names(base_calls)])))
noGA_mclustcalls = do.call(cbind, lapply(psm_results, function(x) adjust_labels_B_to_match_A(base_calls, x$mclust_calls[names(base_calls)])))
colnames(noGA_calls) = lapply(psm_results, function(x) x$name)
colnames(noGA_mclustcalls) = lapply(psm_results, function(x) paste0("mclust_", x$name))
combined_calls = data.frame(cbind(noGA_calls, noGA_mclustcalls))
plot_ari_for_datasets(combined_calls, "noGA")

base_calls = psm_results$with_subtypes_001$calls$cl
psm_results_withGA = psm_results[sapply(psm_results, function(x) !grepl("noGA", x$name))]
calls = do.call(cbind, lapply(psm_results_withGA, function(x) adjust_labels_B_to_match_A(base_calls, x$calls$cl)))
mclustcalls = do.call(cbind, lapply(psm_results_withGA, function(x) adjust_labels_B_to_match_A(base_calls, x$mclust_calls)))
colnames(calls) = lapply(psm_results_withGA, function(x) x$name)
colnames(mclustcalls) = lapply(psm_results_withGA, function(x) paste0("mclust_", x$name))
combined_calls = data.frame(cbind(calls, mclustcalls))
plot_ari_for_datasets(combined_calls, "withGA")
