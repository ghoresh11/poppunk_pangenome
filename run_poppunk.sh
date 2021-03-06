

assemblies_file=assemblies.txt
output_file=poppunk_db_parallel
threads=16
job_name=poppunk
metadata_file=/lustre/scratch118/infgen/team216/gh11/e_coli_collections/FILTERED_MD_FINAL_ALL.tab

## Step 1: activate the python 3 environemnt
conda activate python36

## Step 0: generate a file with all the locations of the assemblies, then run the script below
## and calculate the randindex between this clustering to ST and to MASH
python get_assembly_locs.py > assemblies.txt



## Step 2: Create a DB from all the assemblies
## --- this is the slow step ----
bsub -q parallel -J ${job_name} -R"select[mem>30000] rusage[mem=30000]" -M30000  -G team216 -o ${job_name}.o -e ${job_name}.e -n16 -R"span[hosts=1]" poppunk --create-db --r-files  ${assemblies_file} --output ${output_file} --threads ${threads} --plot-fit 5 --k-step 3  --max-k 31 --min-k 18

## I need to choose an informative kmer to estimate the line. The bigger the genomes, the higher k need to be

## next step: fit the model:
## --- for this number of genomes it was very fast ---

## trying different Ks should change the output of fit model
for K in $(seq 11 3 20 );
do echo $K;
#job_name=${output_file}_${K}
#bsub -q parallel -J ${job_name} -R"select[mem>30000] rusage[mem=30000]" -M30000  -G team216 -o ${job_name}.o -e ${job_name}.e -n16 -R"span[hosts=1]" poppunk --fit-model --distances ${output_file}/${output_file}.dists  --output ${K}_fit --full-db --ref-db ${output_file} --K ${K}
done

## otherwise try to use db-scan:
job_name=poppunk_dbscan  ## I couldn't get this to work, there's an error
bsub -q parallel -J ${job_name} -R"select[mem>40000] rusage[mem=40000]" -M40000  -G team216 -o ${job_name}.o -e ${job_name}.e -n16 -R"span[hosts=1]" poppunk --fit-model --distances ${output_file}/${output_file}.dists  --output dbscan_fit --full-db --ref-db ${output_file} --dbscan


## to refine the model and improve the fit:
## I chose K=11
# didn't manage to get the model refinement to work
# job_name=refine
# bsub -q parallel -J ${job_name} -R"select[mem>30000] rusage[mem=30000]" -M30000  -G team216 -o ${job_name}.o -e ${job_name}.e -n16 -R"span[hosts=1]" poppunk --refine-model --distances ${output_file}/${output_file}.dists --output refine_model --full-db --ref-db ${output_file} --threads ${threads}

## OTHER OUTPUTS:

## calculate rand indexes to find how much the poppunk clusters match the MASH/STs
clusters=mash_clusters.csv,st_cluster.csv
for K in $(seq 5 3 20 );
do
clusters=${clusters},clusters/${K}_fit_clusters.csv
done
python calculate_rand_indices.py --input $clusters


## 1. Extract the pairwise distances from the output files
runjob5 python extract_distances.py --distances ${output_file}/${output_file}.dists --output ${output_file}/${output_file}.dists.out


### after extracting the distances, calculate how big each cluster is and what the core and accessory extract_distances
# are within the cluster in order to run a bunch of roary jobs
## it should take about 10 minutes to run!
conda activate python27
runjob5 python 1_calc_cluster_dists.py --clusters_file ${output_file}/${output_file}_clusters.csv --dists_file ${output_file}/${output_file}.dists.out --metadata_file ${metadata_file}
runjob10 -o all_sizes.o -e all_sizes.e  1_python calc_cluster_dists.py --clusters_file ${output_file}/${output_file}_clusters.csv --dists_file ${output_file}/${output_file}.dists.out --metadata_file ${metadata_file} --out dists_analysis_all --min_cluster_size 1


### run roary
python run_roary.py

##post process the roary outputs by splitting some of the clusters based on the alignment length
python run_postprocessing.py


job_name=gene_extraction
bsub -q parallel -J ${job_name} -R"select[mem>3000] rusage[mem=3000]" -M3000  -G team216 -o ${job_name}.o -e ${job_name}.e -n16 -R"span[hosts=1]" python extract_gene_sequences.py


## At this stage: run-eggnog and also need to add here to compare to known databases so
## that there's an extra column - in addition to the eggnog column, saying
## "AMR", "Phage" or "Vir", to add more specific annotation

### After running roary, this is in 3.1 to extract all the properties of all the genes
## get the properties of all the genes
## most clusters only need 500MB, cluster 2 compeleted with 1GB and 1 with 4GB
## this submits with the correct memory allocation for the different clusters
job_name=props
for i in $(seq 1 1); do
if [ $i -eq 1 ]; then
mem=6000
elif [ $i -eq 2 ]; then
mem=1000
else
mem=500
fi
bsub -J ${job_name}_${i} -R"select[mem>${mem}] rusage[mem=${mem}]" -M${mem}  -G team216 -o ${job_name}_${i}.o -e ${job_name}_${i}.e  python 1_get_gene_props.py --cluster ${i};
done

## after extracting genes, blast them against Phage, Virulence and AMR
job_name=blast_dbs
for i in $(seq 1 39);
#do bsub -q parallel -J ${job_name}_${i} -R"select[mem>5000] rusage[mem=5000]" -M5000  -G team216 -o ${job_name}_${i}.o -e ${job_name}_${i}.e -n8 -R"span[hosts=1]" python blast_against_DBs.py --cluster ${i};
do python blast_against_DBs.py --cluster ${i};
done
