#! /usr/bin/env bash

#Snakefile for analysis of ChIP Seq data from PAO1.
    #Constituitive on and off mutants are E and N, respectively.
    #P is non-chip control.
    #Data from Dr. MJ Schurr.


FASTA = "dbases/PAgenome.fa"
FNA = "dbases/Pseudomonas_aeruginosa_PAO1_107.fna"
CHROMS = "PA_length.txt"
GENES = "DBPA.bed"
WIN = "src/get_tss_windows.sh"

CHIPS, = glob_wildcards("{chip}.fastq")


#12: generate visualization of changes, as list of genes


#11: determine differences between each E, N, and P
    statistical cut-off 


#10: annotate genes of peaks


#9: Get fasta format from peaks


#8: call peaks
    Macs2


#7: pull up peaks in pdf files
rule all:
    input:
        expand("plots/{chip}.pdf", chip = CHIPS)

#6: plot data from the coverage
rule plot_data:
    input:
        "bedfiles/{chip}/coverage.bed"
    output:
        "plots/{chip}.pdf"
    shell:
        "src/plot_data.R {input} {output}"

#5: get total coverage of peaks?
rule get_coverage:
    input:
        windows = "bedfiles/tss_windows.bed",
        bedgraph = "bowtie/{chip}.bedgraph"
    output:
        coverage = "bedfiles/{chip}/coverage.bed"
    shell:
        """
        bedtools map -a {input.windows} \
            -b {input.bedgraph} \
            -c 4 -o mean -null 0 > {output.coverage}
        """

#4:
rule prepare_windows:
    input:
        genes = GENES,
        chroms = CHROMS,
    output:
        "bedfiles/tss_windows.bed"
    shell:
        "src/get_tss_windows.sh {input.genes} {input.chroms} {output}"

#3: create bedgraphs from bowtie map
rule make_bedgraphs:
    input:
        "bowtie/{chip}.bam"
    output:
        "bowtie/{chip}.bedgraph"
    shell:
        """
        bedtools genomecov \
            -ibam {input} \
            -bg \
            -g {CHROMS} > {output}
        """

#2: map with the bowtie index
rule bowtie_mapping:
    input:
        fq = "{chip}.fastq",
        idx = "dbases/bowtie_idx/PAO1.1.bt2"
    output: "bowtie/{chip}.bam"
    params: idx = "dbases/bowtie_idx/PAO1"
    shell:
        """
        bowtie2 \
            -x {params.idx} \
            -U {input.fq} \
            -S {output}.tmp

        samtools sort {output}.tmp > {output}
        samtools index {output}
        """

#1: make the index
rule bowtie_index:
    input: FNA
    output: "dbases/bowtie_idx/PAO1.1.bt2"
    params:
        output_name = "dbases/bowtie_idx/PAO1"
    shell:
        """
        bowtie2-build {input} {params.output_name}
        """
