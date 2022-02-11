source("~/rds/hpc-work/UncertainClustering/kcn25/basis-clustering/R/dpMix_cat.R")

read_matrix <- function(filename) {
    table <- read.table(filename, sep='\t', header=TRUE, row.names=1)
    mat <- as.matrix(table)
    colnames(mat) <- NULL
    rownames(mat) <- NULL
    return (mat)
}

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
    dataset = args[1]
    seed = 1
} else if (length(args)==2) {
    dataset = args[1]
    seed = as.integer(args[2])
}

obsData = read_matrix(paste0("data/", dataset, "/beta.tsv"))
obsVars = read_matrix(paste0("data/", dataset, "/var.tsv"))

directory = paste0("results/", dataset, "/seed", seed)
dir.create(directory, recursive=TRUE)
dpMix(obsData, obsVars, saveFileDir=directory, unique_id=seed, nIts=100000000, thinningFreq=1000)
