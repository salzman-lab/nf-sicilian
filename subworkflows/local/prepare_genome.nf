//
// Uncompress and prepare reference genome files
//

params.genome_options       = [:]
params.index_options        = [:]
params.gffread_options      = [:]
params.star_index_options   = [:]
params.sicilian_createannotator_options   = [:]


include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF } from '../../modules/nf-core/software/gunzip/main'                            addParams( options: params.genome_options          )
include { UNTAR as UNTAR_STAR_INDEX   } from '../../modules/nf-core/software/untar/main'                addParams( options: params.star_index_options      )
include { GFFREAD                     } from '../../modules/nf-core/software/gffread/main'              addParams( options: params.gffread_options         )
include { STAR_GENOMEGENERATE         } from '../../modules/nf-core/software/star/genomegenerate/main'  addParams( options: params.star_index_options      )
include { SICILIAN_CREATEANNOTATOR    } from '../../modules/local/sicilian_createannotator.nf'          addParams( params.sicilian_createannotator_options )

workflow PREPARE_GENOME {

    main:

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gffread_version = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_gffread_version = GFFREAD.out.version
    }

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index   = Channel.empty()
    ch_star_version = Channel.empty()
    if (params.star_index) {
        if (params.star_index.endsWith('.tar.gz') || params.star_index.endswith('.tgz')) {
            ch_star_index = UNTAR_STAR_INDEX ( params.star_index ).untar
        } else {
            ch_star_index = file(params.star_index)
        }
    } else {
        ch_star_index   = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
        ch_star_version = STAR_GENOMEGENERATE.out.version
    }


    //
    if (!params.annotator) {
        SICILIAN_CREATEANNOTATOR ( ch_gtf )
        ch_sicilian_gene_names  = SICILIAN_CREATEANNOTATOR.out.gene_names
        ch_sicilian_splices     = SICILIAN_CREATEANNOTATOR.out.splices
        ch_sicilian_exon_bounds = SICILIAN_CREATEANNOTATOR.out.exon_bounds
    } else {
        ch_sicilian_gene_names  = file(params.annotator)
        ch_sicilian_splices     = file(params.splicesites)
        ch_sicilian_exon_bounds = file(params.exon_bounds)
    }

    emit:
    fasta            = ch_fasta            // path: genome.fasta
    gtf              = ch_gtf              // path: genome.gtf
    star_index       = ch_star_index       // path: star/index/
    sicilian_gene_names = ch_sicilian_gene_names
    sicilian_splices = ch_sicilian_splices
    sicilian_exon_bounds = ch_sicilian_exon_bounds
    star_version     = ch_star_version     // path: *.version.txt
    gffread_version  = ch_gffread_version  // path: *.version.txt
}