#human:/share/database/openData/GRCh38_GENCODE/transcripts_index_salmon
#mus:/datapool/yanzeqin/database/GRCm39_GENCODE/salmon
input=/datapool/yanzeqin/project/ZDW/RNA/merge/MI
output=${input}/03.salmon
index=/share/database/openData/GRCh38_GENCODE/transcripts_index_salmon

mkdir -p ${output}

qstat|grep Eqw|awk '{print $1}' |xargs qdel -j
ls ${input}/03.salmon/*/quant.sf|awk -F "03.salmon" '{print $2}' |awk -F "/" '{print $2}' > sample_finish
cat samples_info|grep -vwf sample_finish |while read l1 l2 l3 l4
do
	
	mkdir -p ${output}/$l1
	outdir=${output}/$l1
	fastq1=${l3}
	fastq2=${l4}
	
	run_salmon=${outdir}/run_${l1}_salmon.sh
	cd $outdir
	echo -e "/share/Data01/yanzeqin/software/snakemake/conda/envs/salmon/bin/salmon quant -i ${index} --gcBias -l A --thread 4 -1 $fastq1 -2 $fastq2 -o $outdir " > $run_salmon
	qsub -clear -cwd -q all_el6.q  -l vf=20G,num_proc=4 $run_salmon
done
