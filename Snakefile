inp_dir = config['input']
out_dir = config['output']
stems = config['stems'].split()
R = ['1', '2']
kits = ['KAPA', 'Collibri']

def get_stranded(wildcards):
    stem = wildcards.stem
    if 'KAPA' in stem:
        return '-s 2'
    elif 'Collibri' in stem:
        return '-s 1'

rule all:
    input:
        out_dir + "/multiqc_report_raw.html",
        out_dir + "/multiqc_report_raw_data.zip",
        expand(inp_dir + "/{stem}_R1_trimmed.fastq.gz", stem=stems),
        expand(inp_dir + "/{stem}_R2_trimmed.fastq.gz", stem=stems),
        out_dir + "/multiqc_report_trimmed.html",
        out_dir + "/multiqc_report_trimmed_data.zip",
        expand(out_dir + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam", stem=stems),
        expand(out_dir + "/STRD/{stem}_str.log", stem=stems),
        expand(out_dir + "/counts/counts_{stem}.log", stem=stems),
        out_dir + "/multiqc_STAR.html",
        out_dir + "/multiqc_featureCounts.html",
        expand("{kit}_DE_genes.csv", kit = kits),
        expand("{kit}_DOWN_DE_genes.csv", kit = kits),
        expand("{kit}_UP_DE_genes.csv", kit = kits),
        expand("{kit}_vulcano.png", kit = kits),
        "venn.png"

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

rule multiqc_STAR:
    input:
        expand(out_dir + "/BAMs/MAPPED_{stem}_Log.final.out", stem = stems)
    output:
        out_dir + "/multiqc_STAR.html",
        out_dir + "/multiqc_STAR_data.zip"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_STAR -z"

rule featureCounts_multiqc:
    input:
        expand(out_dir + "/counts/counts_{stem}.log.summary", stem = stems)
    output:
        out_dir + "/multiqc_featureCounts.html",
        out_dir + "/multiqc_featureCounts_data.zip"
    shell:
        "multiqc {input} -o " + out_dir + " -n multiqc_featureCounts -z"   

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
rule STAR_index:
    output:
        directory(out_dir + "/STAR_genome")
    threads: 2
    params:
        ref = config['ref'],
        gtf = config['gtf']
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} " +
        "--genomeFastaFiles {params.ref} --sjdbGTFfile {params.gtf}"

rule STAR:
    input:
        inp_dir + "/{stem}_R1_trimmed.fastq.gz",
        inp_dir + "/{stem}_R2_trimmed.fastq.gz",
        out_dir + "/STAR_genome"
    output:
        out_dir + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam",
        out_dir + "/BAMs/MAPPED_{stem}_Log.final.out"
    threads: 2
    params:
        pref = out_dir + "/BAMs/MAPPED_{stem}_"
    shell:
        "STAR --runThreadN {threads} --genomeDir {input[2]} " +
        "--readFilesIn {input[0]} {input[1]} --readFilesCommand zcat " +
        "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.pref}"

rule sort:
    input:
        out_dir + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam"
    output:
        out_dir + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    shell:
        "samtools sort -n -o {output} {input}"

rule get_strand:
    input:
        out_dir + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    output:
        out_dir + "/STRD/{stem}_str.log"
    params:
        bed = config['bed']
    shell:
        "infer_experiment.py -r {params.bed} -i {input} > {output}"

rule featureCounts:
    input:
        out_dir + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    output:
        out_dir + "/counts/counts_{stem}.log",
        out_dir + "/counts/counts_{stem}.log.summary"
    params:
        strd = get_stranded,
        gtf = config['gtf']
    shell:
        "featureCounts -p {params.strd} -t exon -g gene_id -a {params.gtf} -o {output[0]} {input}"

rule deseq2:
    input:
        expand(out_dir + "/counts/counts_{stem}.log", stem=stems)
    output:
        expand("{kit}_DE_genes.csv", kit = kits),
        expand("{kit}_DOWN_DE_genes.csv", kit = kits),
        expand("{kit}_UP_DE_genes.csv", kit = kits),
        expand("{kit}_vulcano.png", kit = kits),
        "venn.png"
    shell:
        "Rscript analyse_DE.r {input}"

