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

FASTQs = glob.glob("Fastq_Combined/*R1.fastq.gz")
SAMPLES=[]
for sample in FASTQs:
    x=os.path.basename(sample).replace("_R1.fastq.gz", "")
    SAMPLES.append(x)

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9' ,'10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
SAMPLES = ['109303-128-Tumor-R-LFS']

### ### ### RULES ### ### ###

rule all:
	input:
		expand("bam_input/final/{sample}/{reference}/{sample}.ready.bam", sample=SAMPLES, reference=config['reference']['key'])

# Run BWA-mem on paired end FASTQ file
rule aln_pe:
	input:
		R1 = "Fastq_Combined/{sample}_R1.fastq.gz",
		R2 = "Fastq_Combined/{sample}_R2.fastq.gz"
	output:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.bam"
	params:
		fasta=config['reference']['fasta']
	threads:
		16
	shell:
		"bwa mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools addreplacerg -r 'ID:1.{wildcards.sample}' -r 'PU:1.{wildcards.sample}.1' -r 'PL:illumina' -r 'SM:{wildcards.sample}' -@ {threads} - | samtools sort -@ {threads} -o {output}"

# Split bam file by chromosome
rule split:
    input:
        "bam_input/work/{sample}/{reference}/{sample}_mapped.bam"
    output:
        "bam_input/work/{sample}/{reference}/{sample}_mapped.REF_{chrom}.bam"
    shell:
        "/project/kmaxwell/software/bamtools/build/src/bamtools split -in {input} -reference"

# Index split bam files
rule index:
	input:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.REF_{chrom}.bam"
	output:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.REF_{chrom}.bam.bai"
	shell: 
		"samtools index {input}"

# Mark duplicates in bam files
rule MarkDuplicates:
	input:
		"bam_input/work/{sample}/{reference}/{sample}_mapped.REF_{chrom}.bam"
	output:
		bam="bam_input/work/{sample}/{reference}/mDup.REF_{chrom}.bam",
		metrics="bam_input/final/{sample}/metrics/{reference}/mark_duplicates.table.REF_{chrom}"
	params:
		memory="10240m",
		picard=config['software']['picard']
	shell:
		"java -Xmx{params.memory} -jar {params.picard} MarkDuplicates I={input} O={output.bam} M={output.metrics} CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"

# Find intervals to be realigned
rule RealignerTargetCreator:
	input:
		"bam_input/work/{sample}/{reference}/mDup.REF_{chrom}.bam"
	output:
		"bam_input/work/{sample}/{reference}/IndelRealigner.REF_{chrom}.intervals"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
		thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
		thousandg_gold_standard=config['reference']['thousandG_gold_standard']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T RealignerTargetCreator -I {input} -o {output} -known {params.thousandg_phase1} -known {params.thousandg_gold_standard}"

# Realign problem areas		
rule IndelRealigner:
	input:
		bam="bam_input/work/{sample}/{reference}/mDup.REF_{chrom}.bam",
		targets="bam_input/work/{sample}/{reference}/IndelRealigner.REF_{chrom}.intervals"
	output:
		"bam_input/work/{sample}/{reference}/realign.REF_{chrom}.bam"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
		thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
		thousandg_gold_standard=config['reference']['thousandG_gold_standard']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T IndelRealigner -I {input.bam} -o {output} -targetIntervals {input.targets} -known {params.thousandg_phase1} -known {params.thousandg_gold_standard}"

# first pass of base quality score recalibration
rule FirstPass_BaseRecalibrator:
	input:
		"bam_input/work/{sample}/{reference}/realign.REF_{chrom}.bam"
	output:
		"bam_input/final/{sample}/metrics/{reference}/recal_data.REF_{chrom}.table"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta'],
        thousandg_phase1=config['reference']['thousandG_phase1_indels_b37'],
        thousandg_gold_standard=config['reference']['thousandG_gold_standard'],
		dbsnp_138_b37=config['reference']['dbsnp_138_b37']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T BaseRecalibrator -I {input} -o {output} -knownSites {params.dbsnp_138_b37} -knownSites {params.thousandg_phase1} -knownSites {params.thousandg_gold_standard}"

# Second pass of base quality score recalibration
rule SecondPass_BaseRecalibrator:
	input:
		"bam_input/work/{sample}/{reference}/realign.REF_{chrom}.bam",
		"bam_input/final/{sample}/metrics/{reference}/recal_data.REF_{chrom}.table"
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

# Evaluate and compare base quality score recalibration tables
rule AnalyzeCovariates:
	input:
		before="bam_input/final/{sample}/metrics/{reference}/recal_data.REF_{chrom}.table",
		after="bam_input/final/{sample}/metrics/{reference}/post_recal_data.REF_{chrom}.table"
	output:
		csv="bam_input/final/{sample}/metrics/BQSR.csv",
		pdf="bam_input/final/{sample}/metrics/BQSR.pdf"
	params:
		memory="10240m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T AnalyzeCovariates -before {input.before} -after {input.after} -csv {output.csv} -plots {output.pdf}"

# Write reads that have a passing BQSR to BAM file
rule PrintReads:
	input:
		bam="bam_input/work/{sample}/{reference}/realign.REF_{chrom}.bam",
		bqsr="bam_input/final/{sample}/metrics/{reference}/recal_data.REF_{chrom}.table"
	output:
		"bam_input/work/{sample}/{reference}/recal.REF_{chrom}.bam"
	params:
		memory="1020m",
        gatk=config['software']['gatk'],
		reference=config['reference']['fasta']
	shell:
		"java -Xmx{params.memory} -jar {params.gatk} -R {params.reference} -T PrintReads -I {input.bam} -BQSR {input.bqsr} -o {output}"

# Combine split chromosomes to single bam file
rule combine_bam:
	input:
		expand("bam_input/work/{{sample}}/{reference}/recal.REF_{chrom}.bam", chrom = chromosomes, sample=SAMPLES, reference=config['reference']['key'])
	output:
		"bam_input/final/{sample}/{reference}/{sample}.ready.bam"
	shell:
		"""
		samtools merge {output} {input}
		samtools index {output}
		"""(base) 
