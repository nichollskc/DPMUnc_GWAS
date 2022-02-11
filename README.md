# Setup

Conda environment on HPC - contains all the necessary plotting packages

```
conda activate basis_clustering
```

Download Supplementary Table 6 from Burren et al. 2020 (Genome Medicine)

```
wget https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7687775/bin/13073_2020_797_MOESM6_ESM.csv
```

# Preprocessing

Generate actual data matrices - different subsets of the data which end up in `data/`

```
Rscript prepare_data.R
```

# Running

Run dpMix on the datasets for various seeds, uses `running_dpMix.R` script and the `submit_dpmix.sbatch` script to submit the job using sbatch. Results end up in `results/<subset>/seed<n>`:

```
./run_all.sh
```

To save time reading in data, we can trim the results at this stage rather than after loading. We save these trimmed versions in `trimmed_results/`. In particular we discard the first half of the run and thin the samples even more, so we take only every 10th line of the file (and thus only every 1000th sample of the chain):

```
./trim_datasets.sh
```

The trimming currently isn't done for cluster parameter files, but can easily be adapted to do this too.

# Plotting

Main plotting scripts:

```
Rscript psm_plots.R
Rscript traceplots.R
Rscript quantile_traceplots.R
```

The quantile traceplots are faster and show a good overview of the distributions over time as long as quantiles are well separated from each other.

