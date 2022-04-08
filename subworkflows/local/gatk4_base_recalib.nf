/*
 * Alignment with Cellranger
 */

include { GATK4_BASERECALIBRATOR } from '../../modules/nf-core/modules/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR        } from '../../modules/nf-core/modules/gatk4/applybqsr/main'
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/modules/samtools/index/main'

// Define workflow to subset and index a genome region fasta file
workflow GATK4_BASE_RECALIB {
    take:
        ch_input           // meta, bam, bai, intervals
        reference_fasta
        reference_fai
        reference_dict
        known_sites
        known_sites_tbi

    main:
        ch_versions = Channel.empty()

        GATK4_BASERECALIBRATOR(
            ch_input,
            reference_fasta,
            reference_fai,
            reference_dict,
            known_sites,
            known_sites_tbi
        )

        ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions)

        ch_input_and_tbl = ch_input.join(GATK4_BASERECALIBRATOR.out.table).map{
            meta, bam, bai, intervals, table ->
            [meta, bam, bai, table, intervals]
        }.tap{log2}
        log2.view { "Log 2 base recalib: $it" }

        GATK4_APPLYBQSR(
            ch_input_and_tbl,
            reference_fasta,
            reference_fai,
            reference_dict
        )

        ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions)

        SAMTOOLS_INDEX(GATK4_APPLYBQSR.out.bam)

    emit:
        versions = ch_versions
        bam = GATK4_APPLYBQSR.out.bam
        bai = SAMTOOLS_INDEX.out.bai
}
