
for i in $(seq 1 51); do
if [ $i -eq 50 ]; then
continue
fi
if [ $i -eq 2 ]; then
continue
fi
job_name=${i}_variation
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>1000] rusage[mem=1000]" -M1000 python get_gene_variation.py -i ${i}
done


for i in $(seq 1 1); do
for j in $(seq $((${i}+1)) 2); do
if [ $j -eq 50 ]; then
continue
fi
job_name=${i}_${j}_variation
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>1000] rusage[mem=1000]" -M1000 python get_pairwise_variation.py -i ${i} -j ${j}
done
done
