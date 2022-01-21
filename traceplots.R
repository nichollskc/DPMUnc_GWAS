library(bayesplot)
library(ggplot2)
library(abind)
library(magrittr)

readClusterParams = function(filepath, nDim) {
  values <- matrix(, nrow=0, ncol=nDim)
  
  con = file(filepath, "r")
  continue = TRUE
  while ( continue ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      continue = FALSE
    } else {
      if (grepl("tsv", filepath)) {
        all_entries <- strsplit(line, "   ")[[1]][-1]
      } else {
        all_entries <- strsplit(line, ",")[[1]]
      }
      new_values <- matrix(as.numeric(all_entries), ncol=nDim)
      values <- rbind(values, new_values)
    }
  }
  
  close(con)
  df = data.frame(values)
  df$file = filepath
  df$folder = strsplit(filepath, "/")[[1]][2]
  df$iter = strsplit(filepath, "/")[[1]][3]
  df
}

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

readNumeric = function(filepath) {
  values <- matrix(, nrow=0, ncol=1)
  
  con = file(filepath, "r")
  continue = TRUE
  while ( continue ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      continue = FALSE
    } else {
      new_values <- matrix(as.numeric(line), ncol=1)
      values <- rbind(values, new_values)
    }
  }
  
  close(con)
  df = data.frame(values)
  df$file = filepath
  df$folder = strsplit(filepath, "/")[[1]][2]
  df$iter = strsplit(filepath, "/")[[1]][3]
  df
}

traceplot <- function(datasets, name, directory) {
  print(paste0("Reading in data for " , name))
  K = read_numeric(datasets, "K.txt")
  alpha = read_numeric(datasets, "alpha.tsv")
  pLatent = read_numeric(datasets, "pLatentsGivenClusters.tsv")
  res=abind(alpha, K, pLatent, along=3)
  dimnames(res)=list(Iteration=NULL,
                     Chain=NULL,
                     Parameter=c("alpha","K","pLatent"))
  print(paste0("Making plot for " , name))
  g = mcmc_trace(res)
  print(class(g))
  print(str(g))
  print(paste0("Saving plot for " , name))
  ggsave(paste0(directory, "/plots/trace_", name,  ".png"), g)
}

burren_traceplots <- function(directory) {
  dir.create(paste0(directory, "/plots/"), recursive = TRUE)
  datasets = paste0(directory, "/trimmed_results/with_subtypes_001/seed", 1:5, "/")
  traceplot(datasets, "with_subtypes_001", dir=directory)
}

burren_traceplots("./")
