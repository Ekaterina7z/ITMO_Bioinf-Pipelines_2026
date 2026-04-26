nextflow.enable.dsl=2

params.sra_id = "DRR030302"
params.reads = null
params.reference = null
params.outdir = "results"

workflow {

    reads_ch = Channel.empty()

    if (params.reads) {
        reads_ch = Channel.fromFilePairs(params.reads)
    } else {
        reads_ch = DOWNLOAD_SRA(params.sra_id)
    }

    raw_qc = FASTQC_RAW(reads_ch)

    trimmed_reads = TRIM_READS(reads_ch)

    trimmed_qc = FASTQC_TRIM(trimmed_reads)

    ref_ch = Channel.empty()

    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference)
    } else {
        ref_ch = ASSEMBLE(trimmed_reads)
    }

    mapped = MAP_READS(trimmed_reads, ref_ch)

    coverage = COVERAGE(mapped)
}

process DOWNLOAD_SRA {

    tag "$sra_id"

    input:
    val sra_id

    output:
    tuple path("${sra_id}_1.fastq"), path("${sra_id}_2.fastq")

    script:
    """
    prefetch $sra_id
    fasterq-dump $sra_id --split-files
    """
}

process FASTQC_RAW {

    tag "raw_qc"

    input:
    tuple path(r1), path(r2)

    output:
    path "*_fastqc.html"

    script:
    """
    fastqc $r1 $r2
    """
}

process TRIM_READS {

    tag "trimming"

    input:
    tuple path(r1), path(r2)

    output:
    tuple path("trim_R1.fastq"), path("trim_R2.fastq")

    script:
    """
    fastp \
        -i $r1 \
        -I $r2 \
        -o trim_R1.fastq \
        -O trim_R2.fastq \
        --html fastp_report.html
    """
}

process FASTQC_TRIM {

    tag "trim_qc"

    input:
    tuple path(r1), path(r2)

    output:
    path "*_fastqc.html"

    script:
    """
    fastqc $r1 $r2
    """
}

process ASSEMBLE {

    tag "spades"

    input:
    tuple path(r1), path(r2)

    output:
    path "contigs.fasta"

    script:
    """
    spades.py \
      -1 $r1 \
      -2 $r2 \
      -o spades_out

    cp spades_out/contigs.fasta .
    """
}

process MAP_READS {

    tag "mapping"

    input:
    tuple path(r1), path(r2)
    path ref

    output:
    path "aligned.bam"

    script:
    """
    bwa index $ref

    bwa mem $ref $r1 $r2 | \
        samtools sort -o aligned.bam

    samtools index aligned.bam
    """
}

process COVERAGE {

    tag "coverage"

    input:
    path bam

    output:
    path "coverage.txt"

    script:
    """
    samtools depth $bam > coverage.txt
    """
}

