



RESULT_DIR = config["output"]["relative"]

SAMPLES = {}



with open(config["units"], 'rb') as sinfo:
	for line in sinfo:
		parts = line.decode('utf-8').split()
		sample = parts[0]
		SAMPLES[sample] = [parts[1],parts[2]]

print(SAMPLES)

rule all:
	input:expand([RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_BQSR.bam"],sample=SAMPLES.keys())



rule fastp:
	input:
		rawfq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
		rawfq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
	output:
		read1 =  RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_fastp1.fq.gz",
		read2 =  RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_fastp2.fq.gz",
		fastp_html = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}_report.html",
		fastp_json = RESULT_DIR + "/{sample}/" + config["output"]["fastp"] +  "/{sample}.fastp.json"
	params:
		v1 = "{sample}",
	log: RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_fastp.log"
	shell:
		r'''
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/fastp \
		-i {input.rawfq1} -o {output.read1} -I {input.rawfq2} -O {output.read2} \
		--json {output.fastp_json} --html {output.fastp_html} \
		--report_title {params.v1} \
		--detect_adapter_for_pe --compression 6 --cut_front --cut_tail	 2> {log}
		'''

rule BWA:
	input:
		read1= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_fastp1.fq.gz",
		read2= RESULT_DIR + "/{sample}/" + config["output"]["fastp"] + "/{sample}_fastp2.fq.gz"
	output:
		bam1= temp(RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_Sorted.bam")
	params:
		rg=r"@RG\tID:{sample}"+r"\tSM:{sample}"+r"\tLB:{sample}"+r"\tPU:{sample}"+r"\tPL:ILLUMINA",
		GRCh38= config["databases"]["GRCh38"],
		v1 = "{sample}"
	log: RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}.log"
	shell:
		r'''
		/share/Data01/yanzeqin/software/snakemake/conda/bin/bwa \
		mem -M -t 8 -R "{params.rg}" {params.GRCh38} {input.read1} {input.read2}|/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools \
		sort -m 4G -l 6 -O BAM -T {params.v1} --threads 6 -o {output.bam1} 2> {log}
		'''


rule MarkDuplicates:
	input:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["bwa"] + "/{sample}_Sorted.bam"
	output:
		bam2= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_MD.bam",
		bam3= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_MD_Sorted.bam",
		metrics= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_dup_metrics.txt"
	params:
		v1 = "{sample}",
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_MarkDuplicates.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx100G -Xms16G" MarkDuplicates -I {input.bam1} -O {output.bam2} -M {output.metrics} --TMP_DIR {params.TempDir} 2> {log}
		echo -e "Gatk MarkDuplicates done!" >> {log}
		/share/Data01/yanzeqin/software/snakemake/conda/envs/fastp/bin/samtools sort -l 6 -O bam -T {params.v1} --threads 8 -o {output.bam3} {output.bam2} 2>> {log}
		rm -rf {params.TempDir}
		'''

rule BQSR_T:
	input:bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_MD_Sorted.bam"
	output:table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}.table" 
	params:
		GRCh38= config["databases"]["GRCh38"],
		dbSNP= config["databases"]["dbSNP"],
		G1000= config["databases"]["G1000"],
		ClinVar= config["databases"]["ClinVar"],
		COSMIC1= config["databases"]["COSMIC1"],
		COSMIC2= config["databases"]["COSMIC2"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_BQSRtable.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" BaseRecalibrator -I {input.bam} -R {params.GRCh38} \
		--known-sites {params.dbSNP} --known-sites {params.G1000} --known-sites {params.ClinVar} --known-sites {params.COSMIC1} --known-sites {params.COSMIC2} \
		-O {output.table} --bqsr-baq-gap-open-penalty 30 --tmp-dir {params.TempDir} 2> {log}
		rm -rf {params.TempDir}
		'''

rule ApplyBQSR_T:
	input:
		bam= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_MD_Sorted.bam",
		table= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}.table"	
	output:
		bam1= RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_BQSR.bam"
	params:
		GRCh38= config["databases"]["GRCh38"],
		TempDir= RESULT_DIR + "/{sample}/" + "TempDir"
	log: RESULT_DIR + "/{sample}/" + config["output"]["dup"] + "/{sample}_ApplyBQSR.log"
	shell:
		'''
		mkdir -p {params.TempDir}
		/share/Data01/pengguoyu/bin/gatk --java-options "-Xmx80G -Xms16G" ApplyBQSR -R {params.GRCh38} -I {input.bam} --bqsr-recal-file {input.table} -O {output.bam1} --tmp-dir {params.TempDir} 2> {log}
		rm -rf {params.TempDir}
		'''
