# nf-core/sicilian: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - Read quality control
* [MultiQC](#multiqc) - Aggregate report describing results from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
* [Alignment](#alignment)
  * [STAR](#star) - Fast spliced aware genome alignment and transcriptome quantification
* [Workflow reporting and genomes](#workflow-reporting-and-genomes)
  * [Reference genome files](#reference-genome-files) - Saving reference genome indices/files
  * [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## UMI-tools

[UMI-tools](https://github.com/CGATOxford/UMI-tools) contains tools for dealing with Unique Molecular Identifiers (UMIs)/Random Molecular Tags (RMTs) and single cell RNA-Seq cell barcodes. Two commands are used here: `whitelist` and `extract`.

For further reading and documentation see the [UMI-tools help pages](https://umi-tools.readthedocs.io/en/latest/).

### Whitelist

This command Builds a whitelist of the 'real' cell barcodes.
This is useful for droplet-based single cell RNA-Seq where the identity of the true cell barcodes is unknown. Whitelist can then be used to filter with extract.

### Extract

Flexible removal of UMI sequences from fastq reads.
UMIs are removed and appended to the read name. Any other barcode, for example a library barcode, is left on the read. Can also filter reads by quality or against a whitelist (see above)

**Output files:**

* `umitools/whitelist`
  * `*_whitelist.txt`: Whitelist of true cell barcodes
* `umitools/logs/`
  * `*.log`: Logged output of `umi_tools whitelist`
* `umitools/cell_thresholds/`
  * `*.tsv`: TODO
* `umitools/plots/`
  * `*_cell_barcode_counts.png`: TODO
  * `*_cell_barcode_knee.png`: TODO

<!-- ## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability.

For more information about how to use MultiQC reports, see [https://multiqc.info](https://multiqc.info).

**Output files:**

* `multiqc/`
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats. -->

## STAR

<details markdown="1">
<summary>Output files</summary>

* `star/`
  * `*.Aligned.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the reference genome will be placed in this directory.
  * `*.Aligned.toTranscriptome.out.bam`: If `--save_align_intermeds` is specified the original BAM file containing read alignments to the transcriptome will be placed in this directory.
* `star/log/`
  * `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  * `*.Log.final.out`: STAR alignment report containing the mapping results summary.
  * `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
* `star/unmapped/`
  * `*.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

</details>

[STAR](https://github.com/alexdobin/STAR) is a read aligner designed for splice aware mapping typical of RNA sequencing data. STAR stands for *S*pliced *T*ranscripts *A*lignment to a *R*eference, and has been shown to have high accuracy and outperforms other aligners by more than a factor of 50 in mapping speed, but it is memory intensive.

## SICILIAN

<details markdown="1">
<summary>Output files</summary>

* `sicilian/`
  * `*.Aligned.out.bam`: The original BAM file containing read alignments to the reference genome will be placed in this directory.
* `star/log/`
  * `*.SJ.out.tab`: File containing filtered splice junctions detected after mapping the reads.
  * `*.Log.final.out`: STAR alignment report containing the mapping results summary.
  * `*.Log.out` and `*.Log.progress.out`: STAR log files containing detailed information about the run. Typically only useful for debugging purposes.
* `star/unmapped/`
  * `*.fastq.gz`: If `--save_unaligned` is specified, FastQ files containing unmapped reads will be placed in this directory.

[SICILIAN](https://github.com/salzmanlab/SICILIAN/) (SIngle Cell precIse spLice estImAtioN) is a statistical method for identifying RNA splice junctions using alignments reported from a spliced aligner. SICILIAN is currently implemented for the STAR aligner, and will be adapted to more spliced aligner in the near future.

### Reference genome files

<details markdown="1">
<summary>Output files</summary>

* `genome/`
  * `*.fa`, `*.gtf`, `*.gff`, `*.bed`, `.tsv`: If the `--save_reference` parameter is provided then all of the genome reference files will be placed in this directory.
* `genome/index/`
  * `star/`: Directory containing STAR indices.
  * `sicilian/`: Directory containing SICILIAN reference files
    * `*exon_bounds.pkl`:  is an optional input for SICILIAN and is used to determine whether or not the splice sites in a junction are annotated exon boundaries
    * `*gene_names.pkl`: is a required input for SICILIAN and is used to add gene names to junction ids
    * `*splices.pkl`: is an optional input for SICILIAN and is used to determine whether or not the splice site is annotated in the annotation file

</details>

A number of genome-specific files are generated by the pipeline because they are required for the downstream processing of the results. If the `--save_reference` parameter is provided then these will be saved in the `genome/` directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build new indices so that you can save them somewhere locally. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space.

## Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

**Output files:**

* `pipeline_info/`
  * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Documentation for interpretation of results in HTML format: `results_description.html`.
