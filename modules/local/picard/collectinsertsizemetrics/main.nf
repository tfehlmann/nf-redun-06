process PICARD_COLLECTINSERTSIZEMETRICS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::picard=2.26.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:2.26.10--hdfd78af_0' :
        'quay.io/biocontainers/picard:2.26.10--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.insert_size_metrics.txt"), emit: metrics
    tuple val(meta), path("*.insert_size_histogram.pdf"), emit: histogram
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = 1
    if (!task.memory) {
        log.info '[Picard CollectInsertSizeMetrics] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    picard \\
        CollectInsertSizeMetrics \\
        -Xmx${avail_mem}g \\
        --INPUT $bam \\
        --OUTPUT ${prefix}.insert_size_metrics.txt \\
        --Histogram_FILE ${prefix}.insert_size_histogram.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard CollectInsertSizeMetrics --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
