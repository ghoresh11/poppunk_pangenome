
for i in $(seq 2 51); do
if [ $i -eq 50 ]; then
continue
fi
job_name=${i}_variation
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>500] rusage[mem=500]" -M500 python 3_get_gene_variation.py -i ${i}
done


for i in $(seq 1 49); do
for j in $(seq $((${i}+1)) 51); do
if [ $j -eq 50 ]; then
continue
fi
job_name=${i}_${j}_variation
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>500] rusage[mem=500]" -M500 python 4_get_pairwise_variation.py -i ${i} -j ${j}
done
done
