library(dplyr)

save_subset <- function(rows_to_keep, dir_name) {                                                                          
  delta = reshape(df[rows_to_keep, c("PC", "delta", "label")], idvar="label", timevar = "PC", direction = "wide")
  delta_var = reshape(df[rows_to_keep, c("PC", "var.delta", "label")], idvar="label", timevar = "PC", direction = "wide")

  dir.create(paste0("data/", dir_name), recursive = TRUE)
  write.table(delta, file=paste0("data/", dir_name, "/beta.tsv"), sep='\t', row.names=FALSE)
  write.table(delta_var, file=paste0("data/", dir_name, "/var.tsv"), sep='\t', row.names=FALSE)

  dir.create(paste0("data/", dir_name, "_novar"), recursive = TRUE)
  write.table(delta, file=paste0("data/", dir_name, "_novar/beta.tsv"), sep='\t', row.names=FALSE)
  no_var = delta_var
  no_var[, 2:14] = no_var[, 2:14] * 1e-12
  write.table(no_var, file=paste0("data/", dir_name, "_novar/var.tsv"), sep='\t', row.names=FALSE)
}

df <- read.table("13073_2020_797_MOESM6_ESM.csv", skip=1, sep=',', header=TRUE) %>%
  mutate(fdr.overall.below.thresh = fdr.overall < 0.01,
         label = paste(category.label, trait.label))

# Look for combined traits - we will ignore these and use only their subtypes
is_combined_trait = grepl("comb", df$trait.label, ignore.case=TRUE)

# Traits which are subtypes of these combined traits
is_from_combined_trait_category = df$category.label %in% df[is_combined_trait, "category.label"]
is_subtype = (is_from_combined_trait_category & !is_combined_trait)

# Traits listed in Figure 5 of Astle - probably a sensible subset
# We can also define this subset as below, by excluding any traits that are combinations of other traits
from_astle_fig_5 = c("plt","mpv","pdw","pct","rbc","mcv","hct","mch","mchc","hgb",
                     "rdw","ret","irf","mono","neut","eo","baso","gran","lymph")

# discard any that are percentages and sums
# myeloid and wbc are sums, just with nice names
# hlr (high light scatter reticulocytes) is related to irf (immature reticulocyte fraction) so keep just one
is_complex_astle_trait = df$category.label == "ASTLE" & grepl("_p$|_p_|_sum|myeloid|wbc|hlr", df$trait.label)

is_simple_below_thresh = (df$fdr.overall.below.thresh & !is_complex_astle_trait & !is_combined_trait)
keep_dataset <- is_subtype | is_simple_below_thresh
save_subset(keep_dataset, "with_subtypes_001")

is_geneatlas = grepl("Geneatlas", df$label)
save_subset(keep_dataset & !is_geneatlas, "with_subtypes_noGA_001")

is_subtypes_noGA_noblood_001 = (!is_geneatlas & df$category.label != "ASTLE" & df$category.label != "Cytokines" & !is_combined_trait & df$fdr.overall.below.thresh) | is_subtype
save_subset(is_subtypes_noGA_noblood_001, "with_subtypes_noblood_noGA_001")
