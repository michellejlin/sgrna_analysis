// Use Trimmomatic to trim files, above Q20, minlen of 75
// Trim fastqs
process Trimming { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        tuple val(base), file(R1), file(R2) // from input_read_ch
        file ADAPTERS
        val MINLEN
    output: 
        tuple val(base), file("${base}.trimmed.fastq.gz") //into Trim_out_ch3

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}

    cat *paired.fastq.gz > ${base}.trimmed.fastq.gz
    """
}

// Use Trimmomatic to trim files, above Q20, minlen of 75
// Single end version
process Trimming_SE { 
    container "quay.io/biocontainers/trimmomatic:0.35--6"

    // Retry on fail at most three times 
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 //from input_read_ch
        file ADAPTERS
        val MINLEN
    output: 
        tuple env(base),file("*.trimmed.fastq.gz") //into Trim_out_ch3_SE

    publishDir "${params.OUTDIR}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    base=\$(echo ${R1} | awk -F'_R1' '{print \$1}')
    trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}
    """
}

// Counting sgRNAs
process CountSubgenomicRNAs {
    container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz") //from Trim_out_ch3
        file SGRNAS 
    output:
        file("*stats*")
        tuple val(base),file("*_sgrnas.fastq.gz")
    
    publishDir "${params.OUTDIR}sgRNAs", mode: 'copy'

    script:
    """
    #!/bin/bash
    /usr/local/bin/bbduk.sh in=${base}.trimmed.fastq.gz outm=${base}_sgrnas.fastq.gz ref=${SGRNAS} stats=${base}_sgrnas_stats.txt refstats=${base}_sgrnas_refstats.txt k=40 qhdist=1 -Xmx6g
    """
}

// Mapping sgRNAs to double-check
process MapSubgenomics {
    container "quay.io/biocontainers/bbmap:38.86--h1296035_0"

    // Retry on fail at most three times
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}_sgrnas.fastq.gz") //from CountSubgenomics
        file FULL_SGRNAS 
    output:
        file("${base}_sgrnas_mapped.bam")
    
    publishDir "${params.OUTDIR}sgRNAs", mode: 'copy', pattern: '*.bam'

    script:
    """
    #!/bin/bash
    /usr/local/bin/bbmap.sh in=${base}_sgrnas.fastq.gz outm=${base}_sgrnas_mapped.bam ref=${FULL_SGRNAS} -Xmx6g 2>&1
    """
}
