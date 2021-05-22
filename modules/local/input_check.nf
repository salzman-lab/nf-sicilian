// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(row) {
    sample_id = row[0]
    fastqs = row[1]

    def meta = [:]
    meta.id           = row.sample_id
    // meta.single_end   = row.single_end.toBoolean()
    // meta.strandedness = row.strandedness

    def array = []
    if (!file(fastqs[0]).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (fastqs.size() == 1) {
        array = [ meta, [ file(fastqs[0]) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(fastqs[0]), file(fastqs[1]) ] ]
    }
    return array
}