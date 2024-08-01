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
assemblies = utils.load_assemblies("assemblies.tsv")


def get_unmapped_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"run/unmapped/{wildcards.sample}.unmapped.fastq.gz"]
    return [f"run/unmapped/{wildcards.sample}_R1.unmapped.fastq.gz", f"run/unmapped/{wildcards.sample}_R2.unmapped.fastq.gz"]


def get_trimmed_reads(wildcards):
    sample_info = samples.loc[samples["alias"] == wildcards.sample].iloc[0]
    if pd.isna(sample_info["read2"]): # SE sample
        return [f"run/trimmed/{wildcards.sample}.trimmed.fastq.gz"]
    return [f"run/trimmed/{wildcards.sample}_R1.trimmed.fastq.gz", f"run/trimmed/{wildcards.sample}_R2.trimmed.fastq.gz"]


def original_read_1(wildcards):
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read1.iloc[0]


def original_read_2(wildcards):
    alias = wildcards.sample
    sample_info = samples.loc[samples["alias"] == alias]
    return sample_info.read2.iloc[0]


rule all:
    input:
        "run/counts/raw_counts.tsv"


rule summarize_counts:
    input: 
        counts = expand("run/counts/per_sample/{sample}.tsv", sample=samples.alias),
        samplesheet = config["samplesheet"]
    params: "-m" if config["merge_replicates"] else ""
    log: "run/logs/summarize_counts.log"
    output: "run/counts/raw_counts.tsv"
    shell: "python3 scripts/merge_counts.py -o {output} -i {input.counts} -s {input.samplesheet} -f {params} > {log} 2>&1"


rule featurecounts:
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


rule bam2fastq_pe:
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
    input: "run/unmapped/{sample}.unmapped.bam"
    output: temp("run/unmapped/{sample}.unmapped.sorted.bam")
    log: "run/logs/sort_unmapped/{sample}.log"
    shell: "samtools sort -n {input} -o {output} > {log} 2>&1"


rule extract_unmapped:
    input: "run/alignment/first_pass/{sample}_Aligned.out.bam"
    output: temp("run/unmapped/{sample}.unmapped.bam")
    log: "run/logs/extract_unmapped/{sample}.log"
    shell: "samtools view -f 4 {input} 1> {output} 2> {log}"


rule star_sp:
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

# Misschien outFilterMismatchNmax=<10 ?
rule star_fp:
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
                # --outFilterMultimapNmax 1 \


rule star_build_index:
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


rule fastp_pe:
    input:
        read1 = original_read_1,
        read2 = original_read_2,
        # read1 = Path(config["fastq_dir"]) / "{sample}_R1.fastq.gz",
        # read2 = Path(config["fastq_dir"]) / "{sample}_R2.fastq.gz"
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


rule construct_metagenome:
    input: 
        directories = expand(f"{GENOME_DIR}/{{assembly}}", assembly=assemblies.species),
        annotation_bed = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.annotation.bed", assembly=assemblies.species),
        annotation_gtf = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.annotation.gtf", assembly=assemblies.species),
        fa = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa", assembly=assemblies.species),
        fa_fai = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa.fai", assembly=assemblies.species),
        fa_sizes = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa.sizes", assembly=assemblies.species),
        gaps_bed = expand(f"{GENOME_DIR}/{{assembly}}/{{assembly}}.gaps.bed", assembly=assemblies.species),
    output:
        annotation_bed = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.annotation.bed",
        annotation_gtf = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.annotation.gtf",
        fa = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.fa",
        fa_fai = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.fa.fai",
        fa_sizes = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.fa.sizes",
        gaps_bed = f"{SECOND_ASSEMBLY_DIR}/{config['second_assembly']}.gaps.bed",
    shell:
        """
        bash scripts/merge_assembly.sh {output.annotation_bed} {input.annotation_bed};
        bash scripts/merge_assembly.sh {output.annotation_gtf} {input.annotation_gtf};
        bash scripts/merge_assembly.sh {output.fa} {input.fa};
        bash scripts/merge_assembly.sh {output.fa_fai} {input.fa_fai};
        bash scripts/merge_assembly.sh {output.fa_sizes} {input.fa_sizes};
        bash scripts/merge_assembly.sh {output.gaps_bed} {input.gaps_bed};
        """


def get_accession(wildcards):
    row = assemblies.loc[assemblies["species"] == wildcards.assembly]
    accession = row["name"].iloc[0]
    return accession


rule download_assembly:
    output: 
        temp(directory(f"{GENOME_DIR}/{{assembly}}/")),
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.annotation.gtf",
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.annotation.bed",
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa",
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa.fai",
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.fa.sizes",
        f"{GENOME_DIR}/{{assembly}}/{{assembly}}.gaps.bed",
    params: 
        accession = get_accession,
        assembly = lambda w: w.assembly,
    log: "run/logs/download_assembly/{{assembly}}.log"
    wildcard_constraints: assembly = "(?!\/index)"
    threads: 1
    shell: f"genomepy install {{params.accession}} -g {GENOME_DIR} -l {{params.assembly}} -t {{threads}} -a > {{log}} 2>&1"
