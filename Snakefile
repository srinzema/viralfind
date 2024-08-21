from scripts import utils
from pathlib import Path
import pandas as pd


configfile: "config.yaml"

GENOME_DIR = Path(config["genome_dir"])
FIRST_ASSEBLY_DIR = GENOME_DIR / config["first_assembly"]
SECOND_ASSEMBLY_DIR = GENOME_DIR / config["second_assembly"] 

samples = utils.load_samples(config["samplesheet"], config["fastq_dir"])
paired_end = samples[samples["read2"].notnull()]
single_end = samples[samples["read2"].isna()]

# These following functions exist to link wildcards to existing files.

def get_unmapped_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"run/unmapped/{wildcards.sample}.unmapped.fastq.gz"]
    return [
        f"run/unmapped/{wildcards.sample}_R1.unmapped.fastq.gz", 
        f"run/unmapped/{wildcards.sample}_R2.unmapped.fastq.gz"
    ]


def get_trimmed_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"run/trimmed/{wildcards.sample}.trimmed.fastq.gz"]
    return [
        f"run/trimmed/{wildcards.sample}_R1.trimmed.fastq.gz",
        f"run/trimmed/{wildcards.sample}_R2.trimmed.fastq.gz"
    ]


def original_read_1(wildcards):
    "Gets the first fastq for FASTP."
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read1.iloc[0]


def original_read_2(wildcards):
    "Gets the second fastq for FASTP, not used for paired end samples."
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read2.iloc[0]


def get_accession(wildcards):
    row = assemblies.loc[assemblies["species"] == wildcards.assembly]
    accession = row["name"].iloc[0]
    return accession


rule all:
    """
    This rule depends on the raw_counts.tsv table, 
    which is depends on all other rules.
    """
    input:
        "run/counts/raw_counts.tsv"


rule summarize_counts:
    "Merge the counts of all samples into a single count table."
    input: 
        counts = expand("run/counts/per_sample/{sample}.tsv", sample=samples.alias),
        samplesheet = config["samplesheet"]
    params: "-m" if config["merge_replicates"] else ""
    log: "run/logs/summarize_counts.log"
    output: "run/counts/raw_counts.tsv"
    shell: f"python3 {config['root']}/scripts/merge_counts.py -o {{output}} -i {{input.counts}} -s {{input.samplesheet}} -f {{params}} > {{log}} 2>&1"


rule featurecounts:
    "Run feature counts on the BAM resulting from the second pass."
    input: 
        bam = "run/alignment/second_pass/{sample}_Aligned.out.bam",
        gtf = SECOND_ASSEMBLY_DIR / (config["second_assembly"] + ".annotation.gtf")
    output: "run/counts/per_sample/{sample}.tsv"
    params: paired_end = lambda wildcards: "" if wildcards.sample in single_end.alias else "-p"
    log: "run/logs/featurecounts/{sample}.log"
    threads: 4
    shell:
        """
        featureCounts -a {input.gtf} {input.bam} {params.paired_end} -T {threads} -o {output} > {log} 2>&1
        """

# Below are all rules related to mapping

rule bam2fastq_pe:
    "Converts the unmapped reads into two paired end fastq files."
    input: "run/unmapped/{sample}.unmapped.sorted.bam"
    output: 
        read1 = "run/unmapped/{sample}_R1.unmapped.fastq.gz",
        read2 = "run/unmapped/{sample}_R2.unmapped.fastq.gz"
    log: "run/logs/bam2fastq/{sample}.log"
    threads: 8
    shell:
        """
        samtools fastq -@ {threads} {input} \
        -1 {output.read1} -2 {output.read2} \
        -0 /dev/null -s /dev/null -n > {log} 2>&1
        """


rule bam2fastq_se:
    "Converts the unmapped reads into a single end fastq file."
    input: "run/unmapped/{sample}.unmapped.sorted.bam"
    output: "run/unmapped/{sample}.unmapped.fastq.gz"
    log: "run/logs/bam2fastq/{sample}.log"
    threads: 8
    shell:
        """
        samtools fastq -@ {threads} {input} \
        -o {output} -0 /dev/null -s /dev/null -n > {log} 2>&1
        """


rule sort_unmapped:
    "Sorts all unmapped reads."
    input: "run/unmapped/{sample}.unmapped.bam"
    output: temp("run/unmapped/{sample}.unmapped.sorted.bam")
    log: "run/logs/sort_unmapped/{sample}.log"
    shell: "samtools sort -n {input} -o {output} > {log} 2>&1"


rule extract_unmapped:
    "Extract all unmapped reads from the first pass into temporary BAM file."
    input: "run/alignment/first_pass/{sample}_Aligned.out.bam"
    output: temp("run/unmapped/{sample}.unmapped.bam")
    log: "run/logs/extract_unmapped/{sample}.log"
    shell: "samtools view -f 4 {input} 1> {output} 2> {log}"


rule star_sp:
    """
    Second pass of STAR, maps all unmapped reads against a custom metagenome.
    """
    input:
        index = SECOND_ASSEMBLY_DIR / "index",
        reads = get_unmapped_reads
    output:
        bam = "run/alignment/second_pass/{sample}_Aligned.out.bam"
    params: "run/alignment/second_pass/{sample}_"
    log: "run/logs/alignmend/second_pass/{sample}.log"
    threads: 16
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} \
        --readFilesIn {input.reads} --readFilesCommand zcat \
        --outSAMunmapped Within --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params} > {log} 2>&1
        """


rule star_fp:
    """
    First pass of STAR, runs against the first assembly.
    All unmapped reads are kept in the resulting BAM file.
    """
    input:
        index = FIRST_ASSEBLY_DIR / "index",
        reads = get_trimmed_reads
    output: 
        bam = "run/alignment/first_pass/{sample}_Aligned.out.bam"
    params: "run/alignment/first_pass/{sample}_"
    log: "run/logs/alignment/first_pass/{sample}.log"
    threads: 16
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.index} \
        --readFilesIn {input.reads} --readFilesCommand zcat \
        --outSAMunmapped Within --outSAMtype BAM Unsorted \
        --outFileNamePrefix {params} > {log} 2>&1
        """


rule star_build_index:
    "Build index rule"
    input:
        fa = f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa",
        annotation_gtf = f"{GENOME_DIR}/{{assembly}}/{{assembly}}.annotation.gtf",
    output: directory(f"{GENOME_DIR}/{{assembly}}/index")
    log: "run/logs/star_build_index/" + "{assembly}.log"
    threads: 16
    params: genome_dir = SECOND_ASSEMBLY_DIR
    shell:
        """
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} \
        --genomeFastaFiles {input.fa} --sjdbGTFfile {input.annotation_gtf} > {log} 2>&1
        """

# Fastp rules

rule fastp_pe:
    input:
        read1 = original_read_1,
        read2 = original_read_2,
    output:
        json="run/qc/fastp/{sample}.json",
        html="run/qc/fastp/{sample}.html",
        out1="run/trimmed/{sample}_R1.trimmed.fastq.gz",
        out2="run/trimmed/{sample}_R2.trimmed.fastq.gz",
    log: "run/logs/fastp/{sample}.log"
    threads: 2
    shell:
        "fastp -w {threads} --in1 {input.read1} --in2 {input.read2} --out1 {output.out1} --out2 {output.out2} -h {output.html} -j {output.json} > {log} 2>&1"


rule fastp_se:
    input: original_read_1
    output:
        json = "run/qc/fastp/{sample}.json",
        html = "run/qc/fastp/{sample}.html",
        out = "run/trimmed/{sample}.trimmed.fastq.gz"
    log: "run/logs/fastp/{sample}.log"
    # wildcard_constraints: sample = "(?!.*_R\d)"
    threads: 2
    shell:
        "fastp -w {threads} --in1 {input} --out1 {output.out} -h {output.html} -j {output.json} > {log} 2>&1"


# Rules for constructing a metagenome to map all unmapped reads against

rule construct_metagenome:
    input: config["assembly_file"]
    output:
        outdir = directory(SECOND_ASSEMBLY_DIR),
        fa = SECOND_ASSEMBLY_DIR / f"{config['second_assembly']}.fa",
        annotation_gtf = SECOND_ASSEMBLY_DIR / f"{config['second_assembly']}.annotation.gtf"
    log: f"run/logs/construct_metagenome/{config['second_assembly']}.log"
    threads: 8
    params: 
        outdir = GENOME_DIR,
        outname = config["second_assembly"]
    shell: 
        f"""
        python3 {config['root']}/scripts/construct_metagenome.py \
        {{input}} \
        {{params.outdir}} \
        {{params.outname}} \
        --cores {{threads}} > {{log}} 2>&1
        """