for seed in {1001..1010}
do
    sbatch submit_dpmix.sbatch "with_subtypes_001" $seed
done
