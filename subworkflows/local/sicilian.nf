
/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

params.classinput_options          = [:]
params.glm_options  = [:]
params.annsplices_options = [:]
params.consolidate_options = [:]
params.process_ci_10x_options = [:]
params.postprocess_options = [:]

// Per-sample SICILIAN functions
include { SICILIAN_CLASSINPUT } from '../../modules/local/sicilian/classinput.nf'        addParams( options: params.classinput_options )
include { SICILIAN_GLM        } from '../../modules/local/sicilian/glm.nf'               addParams( options: params.glm_options )
include { SICILIAN_ANNSPLICES } from '../../modules/local/sicilian/annsplices.nf'        addParams( options: params.annsplices_options )

// Postprocessing of SICILIAN output, to consolidate output from mutliple samples into one summary file
include { SICILIAN_CONSOLIDATE    } from '../../modules/local/sicilian/consolidate.nf'  addParams( options: params.consolidate_options )
include { SICILIAN_PROCESS_CI_10X } from '../../modules/local/sicilian/processci10x.nf' addParams( options: params.process_ci_10x_options )
include { SICILIAN_POSTPROCESS    } from '../../modules/local/sicilian/postprocess.nf'  addParams( options: params.postprocess_options )


workflow SICILIAN {
    take:
    bam                 // channel: [ val(meta), bam ]
    sj_out_tab          // channel: [ val(meta), *SJ.out.tab ]
    reads_per_gene      // channel: [ val(meta), *ReadsPerGene.out.tab ]
    chimeric_junction   // channel: [ val(meta), *Chimeric.out.junction ]
    gtf                 // channel: /path/to/genome.gtf
    domain              // channel: /path/to/domain.txt
    annotator           // channel: /path/to/annotator.pkl
    exon_bounds         // channel: /path/to/exon_bounds.pkl
    splices             // channel: /path/to/splices.pkl
    class_input         // channel: [ val(meta), class_input_file ]
    glm_output          // channel: [ val(meta), glm_output_file ]


    main: 
    ch_software_versions = Channel.empty()

    if (!params.skip_classinput) {
        SICILIAN_CLASSINPUT (
            bam,
            gtf,
            annotator,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_CLASSINPUT.out.version.ifEmpty(null))
        ch_class_input = SICILIAN_CLASSINPUT.out.class_input
    } else {
        ch_class_input = class_input
    }

    if (!params.skip_glm) {
        ch_classinput_sjouttab_chimericjunctions_readspergene = ch_class_input.join(
            sj_out_tab
        ).join( 
            chimeric_junction
        ).join(
            reads_per_gene
        )

        SICILIAN_GLM (
            gtf,
            domain,
            exon_bounds,
            splices,
            ch_classinput_sjouttab_chimericjunctions_readspergene
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_GLM.out.version.ifEmpty(null))
        ch_glm_output = SICILIAN_GLM.out.glm_output

        SICILIAN_ANNSPLICES(
            SICILIAN_GLM.out.sicilian_called_splices,
            exon_bounds,
            splices,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_GLM.out.version.ifEmpty(null))
    } else {
        ch_glm_output = glm_output
    }

    SICILIAN_CONSOLIDATE(
        // Take the 2nd (1-index) item, which is the file only, and not the sample id
        ch_glm_output.collect{ it[1] }.dump(tag: 'ch_glm_output_collected')
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_CONSOLIDATE.out.version.ifEmpty(null))

    ch_class_input
        .dump( tag: 'class_input' )
        .multiMap { 
            meta, class_input ->
                concatenation_ids: meta.concatenation_id
                class_inputs: class_input
        }
        .set { class_input_multimap }
    ch_process_ci_concatenation_ids = class_input_multimap.concatenation_ids
        .collect()
        .dump( tag: 'class_input_concatenation_ids_collected' )
    ch_process_ci_class_inputs = class_input_multimap.class_inputs
        .collect()
        .dump( tag: 'class_input_files_collected' )

    if (params.tenx || (!params.smartseq2)) {
            SICILIAN_PROCESS_CI_10X (
            ch_process_ci_concatenation_ids,
            ch_process_ci_class_inputs,
            gtf,
            exon_bounds,
            splices,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_PROCESS_CI_10X.out.version.ifEmpty(null))
        ch_sicilian_junctions_tsv = SICILIAN_PROCESS_CI_10X.out.sicilian_junctions_tsv
    } else {
        // Set an empty sicilian junctions file
        ch_sicilian_junctions_tsv = Channel.empty()
    }


    SICILIAN_POSTPROCESS (
        ch_sicilian_junctions_tsv,
        SICILIAN_CONSOLIDATE.out.glm_consolidated
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_POSTPROCESS.out.version.ifEmpty(null))


    emit:
    version     = ch_software_versions         // path: *.version.txt
}
