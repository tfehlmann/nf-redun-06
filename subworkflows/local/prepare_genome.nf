//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/modules/gunzip/main'

include { BWAMEM2_INDEX                     } from '../../modules/nf-core/modules/bwamem2/index/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/modules/custom/getchromsizes/main'

include { GATK4_CREATESEQUENCEDICTIONARY    } from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'



workflow PREPARE_GENOME {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.genome_fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.genome_fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.genome_fasta)
    }


    //
    // Create chromosome sizes file and fai
    //
    CUSTOM_GETCHROMSIZES (ch_fasta)
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)

    GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
    ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict 
    ch_versions   = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)

    ch_bwa_index = BWAMEM2_INDEX ( ch_fasta ).index
    ch_versions   = ch_versions.mix(BWAMEM2_INDEX.out.versions)



    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    fai              = ch_fai              //    path: genome.fai
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    bwa_index        = ch_bwa_index       //     path: genome.[amb,0123,ann,...]
    dict             = ch_dict
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
