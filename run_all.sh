for seed in {1..5}
do
    sbatch submit_dpmix.sbatch "with_subtypes_001" $seed
done
