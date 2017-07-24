DIR=/mnt/ws/sysbio/roygroup/home/sroy/old-local-home/data/cichilddata/motifs/
for SPECIES in `cut -f1 motifcollection.txt`
do
	TEMPFNAME=${SPECIES/$DIR/}
	TEMPDIR=${TEMPFNAME/allfimo.txt}
	mkdir $TEMPDIR
	echo "mkdir $TEMPDIR"
	echo "egrep "MA0073.1" $SPECIES > $TEMPDIR/allfimo.txt "
	egrep "MA0073.1" $SPECIES > $TEMPDIR/allfimo.txt 
done
