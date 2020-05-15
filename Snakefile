inp_dir = config['input']
out_dir = config['output']
stems = config['stems'].split()
R = ['1', '2']

rule all:
    input:
        out_dir + "/multiqc_report_raw.html",
        out_dir + "/multiqc_report_raw_data.zip",
        expand(inp_dir + "/{stem}_R1_trimmed.fastq.gz", stem=stems),
        expand(inp_dir + "/{stem}_R2_trimmed.fastq.gz", stem=stems),
        out_dir + "/multiqc_report_trimmed.html",
        out_dir + "/multiqc_report_trimmed_data.zip"


rule fastqc:
    input:
        inp_dir + "/{stem}.fastq.gz"
    output:
        out_dir + "/fastqc/{stem}_fastqc.zip",
        out_dir + "/fastqc/{stem}_fastqc.html"
    params:
        out_d = out_dir + "/fastqc/"
    shell:
        "fastqc {input} -o {params.out_d}"

rule multiqc_raw:
    input:
        expand(out_dir + "/fastqc/{stem}_L001_R{R}_001_fastqc.zip", stem = stems, R=R)
    output:
        out_dir + "/multiqc_report_raw.html",
        out_dir + "/multiqc_report_raw_data.zip"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_report_raw -z"

rule multiqc_trimmed:
    input:
        expand(out_dir + "/fastqc/{stem}_R{R}_trimmed_fastqc.zip", stem = stems, R=R)
    output:
        out_dir + "/multiqc_report_trimmed.html",
        out_dir + "/multiqc_report_trimmed_data.zip"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_report_trimmed -z"

rule trim_fastq:
    input:
        inp_dir + "/{stem}_L001_R1_001.fastq.gz",
        inp_dir + "/{stem}_L001_R2_001.fastq.gz"
    output:
        inp_dir + "/{stem}_R1_trimmed.fastq.gz",
        inp_dir + "/{stem}_R2_trimmed.fastq.gz"
    params:
        adapters = config['ADAPT']
    shell:
        "bbduk.sh ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 " +
        "in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]}"

