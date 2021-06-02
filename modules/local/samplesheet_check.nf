// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path samplesheet
    val skip_star
    val skip_classinput
    val skip_glm

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/rnaseq/bin/
    def star_flag = skip_star ? '--skip-star' : ''
    def classinput_flag = skip_classinput ? '--skip-classinput' : ''
    def glm_flag = skip_glm ? '--skip-glm' : ''
    """
    check_samplesheet.py \\
        ${star_flag} \\
        ${classinput_flag} \\
        ${glm_flag} \\
        $samplesheet \\
        samplesheet.valid.csv
    """
}
