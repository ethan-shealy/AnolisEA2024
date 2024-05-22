## From bismark .cov files

for i in /scratch/es88065/cupredict/DSS/*.bismark.cov.gz
do
	cat $i | gunzip | awk '{ print $1"\t"$2"\t"$5+$6"\t"$5 }' > ${i/.bismark.cov.gz/.dss}
done


## From methylkit CpG files

for i in /scratch/es88065/cupredict/methylkit/S*_CpG.txt
do
	tail -n +1 $i | awk '{ print $2"\t"$3"\t"$5"\t"int($5*($6/100)+0.5) }' > ${i/_CpG.txt/.dss}
done