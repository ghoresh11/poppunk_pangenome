



## run pastml
for filename in *.tab; do
	bsub -J $filename -o ${filename}.o -e ${filename}.e  -R"select[mem>500] rusage[mem=500]" -M500 pastml --tree tree_for_treeseg.nwk --data $filename  -o ${filename}.out  --prediction_method ALL 
done


