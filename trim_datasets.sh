$file="with_subtypes_001"
for seed in {1..5}
do
    # Make dir
    OLDDIR="results/${file}/seed${seed}"
    DIR="trimmed_results/${file}/seed${seed}"
    mkdir -p $DIR

    # Trim all files
    awk 'NR > 2500000 && NR % 10 == 0' $OLDDIR/alpha.tsv > $DIR/alpha.tsv
    awk 'NR > 2500000 && NR % 10 == 0' $OLDDIR/clusterAllocations.tsv > $DIR/clusterAllocations.tsv
    awk 'NR > 2500000 && NR % 10 == 0' $OLDDIR/pLatentsGivenClusters.tsv > $DIR/pLatentsGivenClusters.tsv
    # Also find K from cluster allocations, using fact that we always use cluster labels 1..K and so max value in row is K
    awk 'NR > 2500000 && NR % 10 == 0 {max=0; for (i=0;i<=NF;i++) { if($i > max|| i == 0) { max = $i;} } print max}' $OLDDIR/clusterAllocations.tsv > $DIR/K.txt
done
