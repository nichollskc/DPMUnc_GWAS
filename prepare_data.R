data = readRDS("Kath_clustering_data_v3.RDS")
data$Var.Delta = as.numeric(data$Var.Delta)

save_subset <- function(rows_to_keep, dir_name) {
    delta = reshape(data[rows_to_keep, c("PC", "Delta", "Label")], idvar="Label", timevar = "PC", direction = "wide")
    delta_var = reshape(data[rows_to_keep, c("PC", "Var.Delta", "Label")], idvar="Label", timevar = "PC", direction = "wide")
    
    dir.create(paste0("data/", dir_name))
    write.table(delta, file=paste0("data/", dir_name, "/beta.tsv"), sep='\t', row.names=FALSE)
    write.table(delta_var, file=paste0("data/", dir_name, "/var.tsv"), sep='\t', row.names=FALSE)
}

is_myositis = data$First_Author %in% c("Miller", "Rothwell") | data$Label == "Dermatopolymyositis / FinnGen" | data$Label == "Polymyositis / FinnGen"
excluding_combined = (data$Label != "Myositis (European) / Miller")

save_subset(excluding_combined, "myositis_plus_imd_76")

is_myositis_not_combined = is_myositis & excluding_combined
save_subset(is_myositis_not_combined, "myositis_11")

save_subset(!is_myositis, "imd_65")
