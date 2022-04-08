process GATK4_SETNMMDANDUQTAGS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.5.0--hdfd78af_0' :
        'quay.io/biocontainers/gatk4:4.2.5.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(input)
    path  fasta

    output:
    tuple val(meta), path("*.bam"),  emit: bam, optional: true
    tuple val(meta), path("*.bai"),  emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"           ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()

    def avail_mem = 14
    if (!task.memory) {
        log.info '[GATK SetNmMdAndUqTags] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    gatk --java-options "-Xmx${avail_mem}g" SortSam \\
        --INPUT $input \\
        --OUTPUT /dev/stdout \\
        --SORT_ORDER coordinate \\
        --CREATE_INDEX false \\
        --CREATE_MD5_FILE false \\
        | \\
    gatk --java-options "-Xmx${avail_mem}g" SetNmMdAndUqTags \\
            --INPUT /dev/stdin \\
            --OUTPUT ${prefix}.tagged.${file_type} \\
            --CREATE_INDEX true \\
            --CREATE_MD5_FILE true \\
            --REFERENCE_SEQUENCE $fasta \\
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
