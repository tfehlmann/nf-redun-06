/*
 * Alignment with Cellranger
 */

include {CELLRANGER_COUNT} from "../../modules/nf-core/modules/cellranger/count/main.nf"

// Define workflow to subset and index a genome region fasta file
workflow CELLRANGER_ALIGN {
    take:
        cellranger_reference
        ch_fastq

    main:
        ch_versions = Channel.empty()

        // Obtain read counts
        CELLRANGER_COUNT (
             ch_fastq.map{ meta, reads -> [meta + ["gem": meta.id, "samples": [meta.id]], reads]},
             cellranger_reference
        )
        ch_versions = ch_versions.mix(CELLRANGER_COUNT.out.versions)

    emit:
        ch_versions
        cellranger_out  = CELLRANGER_COUNT.out.outs
}