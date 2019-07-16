

for i in $(seq 4 49); do
for j in $(seq $((${i}+1)) 51); do
if [ $j -eq 50 ]; then
continue
fi
job_name=${i}_${j}_merge
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>2000] rusage[mem=2000]" -M2000 -n2 -R"span[hosts=1]" python trial_merge.py --a ${i}/ --b ${j}/ --name_a ${i} --name_b ${j}
done
done


## cluster 1-3 need more memory
for i in $(seq 2 3); do
for j in $(seq $((${i}+1)) 51); do
if [ $j -eq 50 ]; then
continue
fi
job_name=${i}_${j}_merge
bsub -J ${job_name} -G team216 -o ${job_name}.o -e ${job_name}.e -R"select[mem>5000] rusage[mem=5000]" -M5000 -n2 -R"span[hosts=1]" python trial_merge.py --a ${i}/ --b ${j}/ --name_a ${i} --name_b ${j}
done
done
