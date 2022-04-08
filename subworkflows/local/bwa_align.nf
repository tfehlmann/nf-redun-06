include { BWAMEM2_MEM } from "../../modules/nf-core/modules/bwamem2/mem/main"
include { GATK4_SETNMMDANDUQTAGS } from "../../modules/local/gatk-sort-and-setnmmdanduqtags/main"

workflow BWA_ALIGN {
    take:
        ch_fastq
        reference_index
        reference_fasta

    main:
        ch_versions = Channel.empty()

        BWAMEM2_MEM(
            ch_fastq,
            reference_index,
            false
        )
        ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

        GATK4_SETNMMDANDUQTAGS(
            BWAMEM2_MEM.out.bam,
            reference_fasta
        )
        ch_versions = ch_versions.mix(GATK4_SETNMMDANDUQTAGS.out.versions)

    emit:
        versions = ch_versions
        bam      = GATK4_SETNMMDANDUQTAGS.out.bam
}
