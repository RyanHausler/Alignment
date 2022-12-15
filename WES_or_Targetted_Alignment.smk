import shutil
import glob
import os
import errno
import yaml
from collections import defaultdict

configfile : "yaml_files/config.yaml"

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno==errno.EEXIST and os.path.isdir(path):
			pass

FASTQs = glob.glob("FASTQ/*R1.fastq.gz")
SAMPLES=[]
for sample in FASTQs:
    x=os.path.basename(sample).replace("_R1.fastq.gz", "")
    SAMPLES.append(x)

### ### ### RULES ### ### ###

rule all:
	input:
		expand("bam_input/final/{sample}/{reference}/{sample}.ready.bam",sample=SAMPLES,reference=config['reference']['key'])

rule aln_pe:
	input:
		R1 = "FASTQ/{sample}_R1.fastq.gz",
		R2 = "FASTQ/{sample}_R2.fastq.gz"
	output:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.bam"
	params:
		fasta=config['reference']['fasta']
	threads:
		16
	shell:
		"bwa mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools addreplacerg -r 'ID:1.{wildcards.sample}' -r 'PU:1.{wildcards.sample}.1' -r 'PL:illumina' -r 'SM:{wildcards.sample}' -@ {threads} - | samtools sort -@ {threads} -o {output}"

rule index:
	input:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.bam"
	output:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.bam.bai"
	shell: 
		"samtools index {input}"

rule MarkDuplicates:
	input:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.bam"
	output:
		bam="bam_input/work/{sample}/{reference}/mDup.bam",
		metrics="bam_input/final/{sample}/metrics/{reference}/mark_duplicates.table"
	params:
		memory="10240m",
		picard=config['software']['picard']
	shell:
		"java -Xmx{params.memory} -jar {params.picard} MarkDuplicates I={input} O={output.bam} M={output.metrics} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"

rule RealignerTargetCreator:
	input:
		"bam_input/work/{sample}/{reference}/mDup.bam"
	output:
		"bam_input/work/{sample}/{reference}/IndelRealigner.intervals"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
		thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
		thousandg_gold_standard=config['reference']['thousandG_gold_standard']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T RealignerTargetCreator -I {input} -o {output} -known {params.thousandg_phase1} -known {params.thousandg_gold_standard}"

rule IndelRealigner:
	input:
		bam="bam_input/work/{sample}/{reference}/mDup.bam",
		targets="bam_input/work/{sample}/{reference}/IndelRealigner.intervals"
	output:
		"bam_input/work/{sample}/{reference}/realign.bam"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
		thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
		thousandg_gold_standard=config['reference']['thousandG_gold_standard']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T IndelRealigner -I {input.bam} -o {output} -targetIntervals {input.targets} -known {params.thousandg_phase1} -known {params.thousandg_gold_standard}"

rule FirstPass_BaseRecalibrator:#update resources
	input:
		"bam_input/work/{sample}/{reference}/realign.bam"
	output:
		"bam_input/final/{sample}/metrics/{reference}/recal_data.table"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
        thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
        thousandg_gold_standard=config['reference']['thousandG_gold_standard'],
		dbsnp_138_b37=config['reference']['dbsnp_138_b37']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T BaseRecalibrator -I {input} -o {output} -knownSites {params.dbsnp_138_b37} -knownSites {params.thousandg_phase1} -knownSites {params.thousandg_gold_standard}"

rule SecondPass_BaseRecalibrator:
	input:
		"bam_input/work/{sample}/{reference}/realign.bam",
		"bam_input/final/{sample}/metrics/{reference}/recal_data.table"
	output:
		"bam_input/final/{sample}/metrics/{reference}/post_recal_data.table"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
        thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
        thousandg_gold_standard=config['reference']['thousandG_gold_standard'],
		dbsnp_138_b37=config['reference']['dbsnp_138_b37']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T BaseRecalibrator -I {input[0]} -BQSR {input[1]} -o {output} -knownSites {params.dbsnp_138_b37} -knownSites {params.thousandg_phase1} -knownSites {params.thousandg_gold_standard}}"

rule AnalyzeCovariates:
	input:
		before="bam_input/final/{sample}/metrics/{reference}/recal_data.table",
		after="bam_input/final/{sample}/metrics/{reference}/post_recal_data.table"
	output:
		csv="bam_input/final/{sample}/metrics/BQSR.csv",
		pdf="bam_input/final/{sample}/metrics/BQSR.pdf"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T AnalyzeCovariates -before {input.before} -after {input.after} -csv {output.csv} -plots {output.pdf}"

rule PrintReads:
	input:
		bam="bam_input/work/{sample}/{reference}/realign.bam",
		bqsr="bam_input/final/{sample}/metrics/{reference}/recal_data.table"
	output:
		"bam_input/work/{sample}/{reference}/recal.bam"
	params:
		memory="1020m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T PrintReads -I {input.bam} -BQSR {input.bqsr} -o {output}"

rule ValidateSamFile:
	input:
		"bam_input/work/{sample}/{reference}/recal.bam"
	output:
		"bam_input/work/{sample}/{reference}/validation_data.table"
	params:
		memory="10240m",
		picard=config['software']['picard']
	shell:
		"java -Xmx{params.memory} -jar {params.picard} ValidateSamFile I={input} O={output} MODE=SUMMARY"

rule validation_pass:
	input:
		"bam_input/work/{sample}/{reference}/validation_data.table"
	output:
		"bam_input/final/{sample}/metrics/{reference}/validation_data.table"
	run:
		with open(input[0],'r') as file:
			lines=file.read().splitlines()
		if all(not x.startswith('ERROR') for x in lines):
			shutil.copyfile(input[0],output[0])
		else:
			for x in lines:
				if x.startswith('ERROR'):
					print(x)

rule ready_bam:
	input:
		bam="bam_input/work/{sample}/{reference}/recal.bam",
		table="bam_input/final/{sample}/metrics/{reference}/validation_data.table"
	output:
		"bam_input/final/{sample}/{reference}/{sample}.ready.bam"
	shell:
		"""
		rsync -v {input.bam} {output}
		samtools index {output}
		"""
