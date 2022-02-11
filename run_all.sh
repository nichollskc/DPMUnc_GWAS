for seed in {1001..1010}
do
    sbatch submit_dpmix.sbatch "with_subtypes_001" $seed
    sbatch submit_dpmix.sbatch "with_subtypes_001_novar" $seed
#    sbatch submit_dpmix.sbatch "with_subtypes_noGA_001" $seed
    sbatch submit_dpmix.sbatch "with_subtypes_noGA_001_novar" $seed
done
