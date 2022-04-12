/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowExample.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
    params.input,
    params.genome_fasta,
    params.known_indels,
    params.known_snps,
    params.known_snps_tbi,
    params.known_indels_tbi
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}

// genome ref params
ch_genome_fasta = params.genome_fasta ? file(params.genome_fasta) : []
ch_known_indels = params.known_indels ? Channel.fromPath(params.known_indels).collect() : []
ch_known_indels_tbi = params.known_indels_tbi ? Channel.fromPath(params.known_indels_tbi).collect()   : []
ch_known_snps = params.known_snps ? Channel.fromPath(params.known_snps).collect()   : []
ch_known_snps_tbi = params.known_snps_tbi ? Channel.fromPath(params.known_snps_tbi).collect()   : []
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { PICARD_COLLECTALIGNMENTSUMMARYMETRICS } from '../modules/local/picard/collectalignmentsummarymetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../modules/local/picard/collectinsertsizemetrics/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { PREPARE_GENOME     } from '../subworkflows/local/prepare_genome'
include { BWA_ALIGN          } from '../subworkflows/local/bwa_align'
include { GATK4_BASE_RECALIB } from '../subworkflows/local/gatk4_base_recalib'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { GATK4_MARKDUPLICATES        } from '../modules/nf-core/modules/gatk4/markduplicates/main'
include { SAMTOOLS_DEPTH              } from '../modules/nf-core/modules/samtools/depth/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow EXAMPLE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            meta.read_group = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:illumina"
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // prepare reference
    PREPARE_GENOME ()
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)


    CUTADAPT(
        CAT_FASTQ.out.reads
    )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)


    BWA_ALIGN(
        CUTADAPT.out.reads,
        PREPARE_GENOME.out.bwa_index,
        PREPARE_GENOME.out.fasta
    )

    ch_versions = ch_versions.mix(BWA_ALIGN.out.versions)

    GATK4_MARKDUPLICATES(
        BWA_ALIGN.out.bam,
    )

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    ch_bam_bai = GATK4_MARKDUPLICATES.out.bam.tap{log1}.map{ meta, bam -> [meta, bam, [], []]}.tap{log2}

    GATK4_BASE_RECALIB(
        ch_bam_bai,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.dict,
        ch_known_indels.concat(ch_known_snps),
        ch_known_indels_tbi.concat(ch_known_snps_tbi)
    )

    ch_versions = ch_versions.mix(GATK4_BASE_RECALIB.out.versions)

    PICARD_COLLECTINSERTSIZEMETRICS(
        GATK4_BASE_RECALIB.out.bam
    )

    ch_versions = ch_versions.mix(PICARD_COLLECTINSERTSIZEMETRICS.out.versions)

    PICARD_COLLECTALIGNMENTSUMMARYMETRICS(
        GATK4_BASE_RECALIB.out.bam,
        PREPARE_GENOME.out.fasta,
        PREPARE_GENOME.out.dict
    )

    ch_versions = ch_versions.mix(PICARD_COLLECTALIGNMENTSUMMARYMETRICS.out.versions)

    SAMTOOLS_DEPTH(
        GATK4_BASE_RECALIB.out.bam
    )

    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
