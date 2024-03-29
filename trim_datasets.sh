#file="with_subtypes_001"
#for seed in {1001..1010}
#do
#    # Make dir
#    OLDDIR="results/${file}/seed${seed}"
#    DIR="trimmed_results/${file}/seed${seed}"
#    mkdir -p $DIR
#
#    # Trim all files
#    LINE_CHOICE="NR <= 650000"
##    awk "${LINE_CHOICE}" $OLDDIR/alpha.tsv > $DIR/alpha.tsv
#    awk "${LINE_CHOICE}" $OLDDIR/clusterAllocations.tsv > $DIR/clusterAllocations.tsv
##    awk "${LINE_CHOICE}" $OLDDIR/pLatentsGivenClusters.tsv > $DIR/pLatentsGivenClusters.tsv
#    # Also find K from cluster allocations, using fact that we always use cluster labels 1..K and so max value in row is K
##    awk "${LINE_CHOICE} {max=0; for (i=0;i<=NF;i++) { if(\$i > max|| i == 0) { max = \$i;} } print max}" $OLDDIR/clusterAllocations.tsv > $DIR/K.txt
#done

file="with_subtypes_001"
for seed in {1001..1010}
do
    OLDDIR="results/${file}/seed${seed}"
    DIR="trimmed_results/${file}/seed${seed}"
    mkdir -p $DIR

    #alpha.tsv  clusterAllocations.tsv  clusterMeans.tsv  clusterVars.tsv  K.txt  latentObservations.tsv  pLatentsGivenClusters.tsv

    LINE_CHOICE="NR % 10 == 0 && NR <= 650000"
    awk "${LINE_CHOICE}" $OLDDIR/alpha.tsv > $DIR/alpha.tsv
    awk "${LINE_CHOICE}" $OLDDIR/clusterAllocations.tsv > $DIR/clusterAllocations.tsv
    awk "${LINE_CHOICE}" $OLDDIR/clusterMeans.tsv > $DIR/clusterMeans.tsv
    awk "${LINE_CHOICE}" $OLDDIR/clusterVars.tsv > $DIR/clusterVars.tsv
    awk "${LINE_CHOICE}" $OLDDIR/latentObservations.tsv > $DIR/latentObservations.tsv
    awk "${LINE_CHOICE}" $OLDDIR/pLatentsGivenClusters.tsv > $DIR/pLatentsGivenClusters.tsv

    # Also find K from cluster allocations, using fact that we always use cluster labels 1..K and so max value in row is K
    awk "${LINE_CHOICE} {max=0; for (i=0;i<=NF;i++) { if(\$i > max|| i == 0) { max = \$i;} } print max}" $OLDDIR/clusterAllocations.tsv > $DIR/K.txt
done
