import os

configfile: 'config.yaml'

reads = ['1', '2']
fastqc = 'fastqc'
multiqc = 'multiqc'
trimmed = 'trimmed'

rule all:
    input:
        os.path.join(multiqc, trimmed, 'multiqc_report.html'),
        os.path.join(multiqc, 'multiqc_report.html')

rule fastqc:
    input:
        expand(os.path.join('data', '{sample}_R{read}_001.fastq.gz'), sample=config['samples'], read=reads)
    output:
        expand(os.path.join(fastqc, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
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
        expand(os.path.join('data', '{sample}.fastq.gz'), sample=config['names'].split())
    output:
        expand(os.path.join(trimmed, '{sample}.fastq'), sample=config['names'].split())
    shell:"""
        bbduk.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[2]} in2={input[3]} out1={output[2]} out2={output[3]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[4]} in2={input[5]} out1={output[4]} out2={output[5]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[6]} in2={input[7]} out1={output[6]} out2={output[7]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[8]} in2={input[9]} out1={output[8]} out2={output[9]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[10]} in2={input[11]} out1={output[10]} out2={output[11]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[12]} in2={input[13]} out1={output[12]} out2={output[13]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        bbduk.sh in1={input[14]} in2={input[15]} out1={output[14]} out2={output[15]} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10
        """

rule fastqc_trimmed:
    input:
        expand(os.path.join(trimmed, '{sample}_R{read}_001.fastq'), sample=config['samples'], read=reads)
    output:
        expand(os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    shell:
        "fastqc {input} -o {fastqc}/{trimmed}"

rule multiqc_trimmed:
    input:
        expand(os.path.join(fastqc, trimmed, '{sample}_R{read}_001_fastqc.html'), sample=config['samples'], read=reads)
    output:
        os.path.join(multiqc, trimmed, 'multiqc_report.html')
    shell:
        "multiqc {fastqc} -o {multiqc}/{trimmed}"