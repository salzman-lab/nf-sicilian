// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_ALIGN {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process),  meta:[:], publish_by_meta:[]) }

    // Note: 2.7X indices incompatible with AWS iGenomes.
    conda (params.enable_conda ? 'bioconda::star=2.7.5a' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/star:2.7.5a--0'
    } else {
        container 'quay.io/biocontainers/star:2.7.5a--0'
    }

    input:
    tuple val(sample_id), path(reads)
    path  index
    path  gtf

    output:
    tuple val(sample_id), path('*d.out.bam')       , emit: bam
    tuple val(sample_id), path('*Log.final.out')   , emit: log_final
    tuple val(sample_id), path('*Log.out')         , emit: log_out
    tuple val(sample_id), path('*Log.progress.out'), emit: log_progress
    path  '*.version.txt'                     , emit: version

    tuple val(sample_id), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(sample_id), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(sample_id), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(sample_id), path('*fastq.gz')               , optional:true, emit: fastq
    tuple val(sample_id), path('*SJ.out.tab')             , optional:true, emit: sj_out_tab
    tuple val(sample_id), path('*Chimeric.out.junction')  , optional:true, emit: chimeric_out_junction
    tuple val(sample_id), path('*ReadsPerGene.out.tab')   , optional:true, emit: reads_per_gene


    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${sample_id}${options.suffix}" : "${sample_id}"
    def ignore_gtf = params.star_ignore_sjdbgtf ? '' : "--sjdbGTFfile $gtf"
    def seq_center = params.seq_center ? "--outSAMattrRGline ID:$prefix 'CN:$params.seq_center' 'SM:$prefix'" : "--outSAMattrRGline ID:$prefix 'SM:$prefix'"
    def out_sam_type = (options.args.contains('--outSAMtype')) ? '' : '--outSAMtype BAM Unsorted'
    def mv_unsorted_bam = (options.args.contains('--outSAMtype BAM Unsorted SortedByCoordinate')) ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''
    def reads_v2 = params.tenx ? "${reads[1]}" : "${reads}"
    if (params.tenx) {
        if (params.skip_umitools) {
            // Skipping umi tools, so providing already-extracted R2 reads. No R1 at all --> take the first one
            reads_v2 = "${reads[0]}"
        } else {
            // Did UMI tools, and only R2 has sequence to align since R1 is empty after UMI tools extract
            reads_v2 = "${reads[1]}"
        }
    } else {
        reads_v2 = "${reads}"
    }
    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads_v2  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $out_sam_type \\
        $ignore_gtf \\
        $seq_center \\
        $options.args

    $mv_unsorted_bam

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi

    STAR --version | sed -e "s/STAR_//g" > ${software}.version.txt
    """
}
