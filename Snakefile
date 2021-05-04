configfile: "config.yaml"
rule all:
    input:
        'multiqc_report.html'

rule green_fastqc:
    input:
        # fastq = expand('data/{sample}.fastq.gz, sample=config["samples"])'
        expand('data/{sample}.fastq.gz', sample=config["samples"])
    output:
        expand('green_fastqc/{sample}_fastqc.html', sample=config["samples"]),
        expand('green_fastqc/{sample}_fastqc.zip', sample=config["samples"])
    shell:
        'fastqc {input} -o green_fastqc'

rule green_MultiQC:
    input:
        expand('green_fastqc/{sample}_fastqc.zip', sample=config["samples"])
    output:
        'green_multiqc_report.html'
    shell:
        'multiqc {input} -n {output}'

# rule trimming:
#     input:
#         in1 = expand('data/{stem}_R1_001.fastq.gz', stem=config["stems"]),
#         in2 = expand('data/{stem}_R2_001.fastq.gz', stem=config["stems"])
#     output:
#         out1 = expand('trimmed/{stem}_R1_001.fq', stem=config["stems"]),
#         out2 = expand('trimmed/{stem}_R2_001.fq', stem=config["stems"])
#     shell:
#         'bbduk.sh in1={input.in1} in2={input.in2} out1={output.out1}' +
#         ' out2={output.out2} ref=adapters.fa ktrim=r k=23 mink=11 hdist=1' +
#         ' tpe tbo qtrim=r trimq=10'

rule trimming:
    input:
        in1 = expand('data/{stem}.fastq.gz', stem=config["names"].split())
    output:
        out1 = expand('trimmed/{stem}.fq', stem=config["names"].split())
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

rule fastqc:
    input:
        expand('trimmed/{stem}.fq', stem=config["names"].split())
    output:
        expand('fastqc/{sample}_fastqc.html', sample=config["samples"]),
        expand('fastqc/{sample}_fastqc.zip', sample=config["samples"])
    shell:
        'fastqc {input} -o fastqc'

rule MultiQC:
    input:
        expand('fastqc/{sample}_fastqc.zip', sample=config["samples"])
    output:
        'multiqc_report.html'
    shell:
        'multiqc {input} -n {output}'