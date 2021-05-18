import os
from os.path import expanduser

configfile: 'config.yaml'

reads = ['1', '2']
fastqc = 'fastqc'
multiqc = 'multiqc'
trimmed = 'trimmed'

rule all:
    input:
        expand(os.path.join(fastqc, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc, 'multiqc_report.html'),
        expand(os.path.join(trimmed, '{sample}_R1_001.fastq'), sample=config['samples']),
        expand(os.path.join(trimmed, '{sample}_R2_001.fastq'), sample=config['samples']),
        expand(os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads),
        os.path.join(multiqc, trimmed, 'multiqc_report.html')
        

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

