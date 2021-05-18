import os
from os.path import expanduser

configfile: 'config.yaml'

reads = ['1', '2']
fastqc = 'fastqc'
multiqc = 'multiqc'
trimmed = 'trimmed'
ref = 'ref_annot'
index = 'genome_index'
aligned = 'aligned'
samtools = 'samtools_sorted'
fcounts = 'fcounts'
deseq = 'deseq'


rule all:
    input:
        expand(os.path.join(fastqc, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc, 'multiqc_report.html'),
        expand(os.path.join(trimmed, '{sample}_R1_001.fastq'), sample=config['samples']),
        expand(os.path.join(trimmed, '{sample}_R2_001.fastq'), sample=config['samples']),
        expand(os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc, trimmed, 'multiqc_report.html'),
        os.path.join(index, 'genomeParameters.txt'),
        os.path.join(index, 'Genome'),
        expand(os.path.join(aligned, '{sample}', 'Aligned.sortedByCoord.out.bam'), sample=config['samples']),
        expand(os.path.join(samtools, '{sample}.bam'), sample=config['samples']),
        os.path.join(fcounts, 'fcount_result.txt'),
        expand(os.path.join(deseq, '{prep}_DE_volcano.pdf'), prep=config['prep']),
        expand(os.path.join(deseq, '{prep}_DE_results.csv'), prep=config['prep'])        
        

rule fastqc:
    input:
        os.path.join('data', '{sample}_R{read}_001.fastq.gz')
    output:
        os.path.join(fastqc, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} -o {fastqc}"

rule multiqc:
    input:
        expand(os.path.join(fastqc, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    output:
        os.path.join(multiqc, 'multiqc_report.html')
    shell:
        "multiqc {fastqc} -o {multiqc}"

rule trimming:
    input:
        in1 = os.path.join('data', '{sample}_R1_001.fastq.gz'),
        in2 = os.path.join('data', '{sample}_R2_001.fastq.gz')
    output:
        out1 = os.path.join(trimmed, '{sample}_R1_001.fastq'),
        out2 = os.path.join(trimmed, '{sample}_R2_001.fastq')
    shell:
        "bbduk.sh in1={input.in1} out1={output.out1} in2={input.in2} out2={output.out2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

rule fastqc_trimmed:
    input:
        os.path.join(trimmed, '{sample}_R{read}_001.fastq')
    output:
        os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html')
    shell:
        "fastqc {input} -o {fastqc}/{trimmed}"

rule multiqc_trimmed:
    input:
        expand(os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    output:
        os.path.join(multiqc, trimmed, 'multiqc_report.html')
    shell:
        "multiqc {fastqc} -o {multiqc}/{trimmed}"

rule generate_index:
    input:
        fa = os.path.join(ref, 'chr19_20Mb.fa'),
        gtf = os.path.join(ref, 'chr19_20Mb.gtf')
    output:
        os.path.join(index, 'genomeParameters.txt'),
        os.path.join(index, 'Genome')
    shell:
        "STAR --genomeDir {index} --runMode genomeGenerate --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 11 --outFileNamePrefix star_files"

rule align:
    input:
        os.path.join(index, 'genomeParameters.txt'),
        os.path.join(index, 'Genome'),
        in1 = os.path.join(trimmed, '{sample}_R1_001.fastq'),
        in2 = os.path.join(trimmed, '{sample}_R2_001.fastq'),
        gtf = os.path.join(ref, 'chr19_20Mb.gtf')
    output:
        os.path.join(aligned, '{sample}', 'Aligned.sortedByCoord.out.bam')
    shell:
        "STAR --readFilesIn {input.in1} {input.in2} --genomeDir {index} --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} --outFileNamePrefix {aligned}/{wildcards.sample}/"



rule sort:
    input:
        os.path.join(aligned, '{sample}', 'Aligned.sortedByCoord.out.bam')
    output:
        os.path.join(samtools, '{sample}.bam')
    shell:
        "samtools sort -n -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"

def fcounts_inputs(wildcards):
    files = expand(os.path.join(samtools, '{sample}.bam'), sample=config['samples'])
    return files

rule feature_counts:
    input:
        fcounts_inputs
    output:
        os.path.join(fcounts, 'fcount_result.txt')
    shell:
        "featureCounts {input} -p -t exon -g gene_id -a {ref}/chr19_20Mb.gtf -o {output} -s 1"

rule deseq:
    input:
        os.path.join(fcounts, 'fcount_result.txt'),
    params:
        deseq,
        config['prep']
    output:
        os.path.join(deseq, '{prep}_DE_volcano.pdf'),
        os.path.join(deseq, '{prep}_DE_results.csv')
    script:
        "deseq2.R"

