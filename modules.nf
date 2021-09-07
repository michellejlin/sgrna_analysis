#!/usr/bin/env nextflow

// Adapted from TAYLOR (https://github.com/greninger-lab/covid_swift_pipeline)

// Using the Nextflow DSL-2 to account for the logic flow of this workflow
nextflow.preview.dsl=2

// Print help message
def helpMessage() {
    log.info"""
    Usage: 
    An example command for running the pipeline is as follows:
    nextflow run michellejlin/sgrna_analysis -resume -with-docker ubuntu:18.04 --INPUT example/ --OUTDIR example/output/
    
    Parameters:
        --INPUT         Input folder where all fastqs are located.
                        ./ can be used for current directory.
                        Fastqs should all be gzipped. This can be done with the command gzip *.fastq. [REQUIRED]
        --OUTDIR        Output directory. [REQUIRED]
        --SINGLE_END    Optional flag for single end reads. By default, this pipeline does 
                        paired-end reads.
        --MIN_LEN       Set minimum length for Trimmomatic. Default is 75.
        -with-docker ubuntu:18.04   [REQUIRED]
        -resume [RECOMMENDED]
        -profile        Specify which profile to run. For AWS, run with -profile cloud_big. For large memory-intensive runs on AWS, run with -profile cloud_bigger.
        
    """.stripIndent()
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*          SET UP CONFIGURATION VARIABLES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Initializing flags
params.INPUT = false
params.OUTDIR= false
params.SINGLE_END = false
params.MIN_LEN = 75

// Checking for argument validity
// Throw error if --INPUT not set
if (params.INPUT == false) {
    println( "Must provide an input directory with --INPUT") 
    exit(1)
}
// Make sure INPUT ends with trailing slash
if (!params.INPUT.endsWith("/")){
    println("Make sure your input directory ends with trailing slash.")
   exit(1)
}
// Throw error if --OUTDIR not set
if (params.OUTDIR == false) {
    println( "Must provide an output directory with --OUTDIR") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.OUTDIR.endsWith("/")){
   println("Make sure your output directory ends with trailing slash.")
   exit(1)
}

// Setting up files 
REFERENCE_FASTA = file("${baseDir}/NC_045512.2.fasta")
REFERENCE_FASTA_FAI = file("${baseDir}/NC_045512.2.fasta.fai")
ADAPTERS = file("${baseDir}/All_adapters.fa")
SGRNAS = file("${baseDir}/sgRNAs_60.fasta")
FULL_SGRNAS=file("${baseDir}/sgRNAs.fasta")

// Import processes 
include { Trimming } from './modules.nf'
include { Trimming_SE } from './modules.nf' 
include { CountSubgenomicRNAs } from './modules.nf'
include { MapSubgenomics } from './modules.nf'

// Import reads depending on single end vs. paired end
if(params.SINGLE_END == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.INPUT}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.INPUT} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.INPUT}*_R1.fastq.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    // Paired end first few steps
    if(params.SINGLE_END == false) {
        Trimming (
            input_read_ch, 
            ADAPTERS,
            params.MIN_LEN
        )
        CountSubgenomicRNAs (
            Trimming.out[2],
            SGRNAS
        )
        MapSubgenomics (
            CountSubgenomicRNAs.out[1],
            FULL_SGRNAS
        )
    } else {
    // Single end first few steps
        Trimming_SE (
            input_read_ch,
            ADAPTERS,
            params.MIN_LEN
        )
        CountSubgenomicRNAs (
           Trimming_SE.out[2],
           SGRNAS
        )
        MapSubgenomics (
           CountSubgenomicRNAs.out[1],
           FULL_SGRNAS
        )
    }
}
