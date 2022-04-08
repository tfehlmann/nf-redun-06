process PICARD_COLLECTALIGNMENTSUMMARYMETRICS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta
    path dict

    output:
    tuple val(meta), path("*.alignment_metrics.txt"), emit: metrics
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 1
    if (!task.memory) {
        log.info '[Picard CollectAlignmentSummaryMetrics] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        CollectAlignmentSummaryMetrics \\
        -Xmx${avail_mem}g \\
        --INPUT $bam \\
        --REFERENCE_SEQUENCE $fasta \\
        --OUTPUT ${prefix}.alignment_metrics.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectAlignmentSummaryMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
