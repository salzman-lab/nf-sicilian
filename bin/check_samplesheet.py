#!/usr/bin/env python

import os
import sys
import errno
import argparse

FASTQ_COLUMNS = "fastq_1,fastq_2"
FASTQ_COLUMNS_LIST = FASTQ_COLUMNS.split(",")
STAR_COLUMNS = "bam,sj_out_tab,reads_per_gene,chimeric_junction"
STAR_COLUMNS_LIST = STAR_COLUMNS.split(",")
CLASSINPUT_COL = "class_input"
GLM_COL = "glm_output"

SKIP_STAR_COLUMNS = STAR_COLUMNS_LIST
SKIP_CLASSINPUT_COLUMNS = STAR_COLUMNS_LIST + [CLASSINPUT_COL]
SKIP_GLM_COLUMNS = SKIP_CLASSINPUT_COLUMNS + [GLM_COL]

# Cribbed from https://github.com/nf-core/rnaseq/blob/master/bin/check_samplesheet.py
def parse_args(args=None):
    Description = (
        "Reformat czbiohub/nf-sicilian samplesheet file and check its contents."
    )
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument(
        "--skip-star",
        action="store_true",
        help=f"Skipping alignment, so not looking for only '{FASTQ_COLUMNS}' columns. "
        f"Looking for '{FASTQ_COLUMNS},{STAR_COLUMNS}' columns",
    )
    parser.add_argument(
        "--skip-classinput",
        action="store_true",
        help=f"Skipping alignment, so not looking for only '{FASTQ_COLUMNS}' columns. "
        f"Looking for '{FASTQ_COLUMNS},{STAR_COLUMNS},{CLASSINPUT_COL}' columns",
    )
    parser.add_argument(
        "--skip-glm",
        action="store_true",
        help=f"Skipping alignment, so not looking for '{FASTQ_COLUMNS}' columns. "
        f"Looking for '{FASTQ_COLUMNS},{STAR_COLUMNS},{CLASSINPUT_COL},{GLM_COL}' columns",
    )
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = f"ERROR: Please check samplesheet -> {error}"
    if context != "" and context_str != "":
        error_str = f"ERROR: Please check samplesheet -> {error}\n{context.strip()}: '{context_str.strip()}'"
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out, skip_star, skip_classinput, skip_glm):
    """
    This function checks that the samplesheet follows the following structure:

    sample_id,fastq_1,fastq_2,single_end,strandedness,concatenation_id
    SAMPLE_PE_RUN1,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,false,forward,SAMPLE_PE
    SAMPLE_PE_RUN1,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz,false,forward,SAMPLE_PE
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,,false,forward,SAMPLE_SE

    For an example see:
    TODO: Create new example
    https://github.com/nf-core/test-datasets/blob/rnaseq/samplesheet/v3.1/samplesheet_test.csv
    """

    sample_mapping_dict = {}

    additional_cols = FASTQ_COLUMNS_LIST
    if skip_star:
        additional_cols += SKIP_STAR_COLUMNS
    if skip_classinput:
        additional_cols += [CLASSINPUT_COL]
    if skip_glm:
        additional_cols += [GLM_COL]

    do_alignment = not (skip_star or skip_classinput or skip_glm)

    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = [
            "sample_id",
            "strandedness",
            "concatenation_id",
        ] + additional_cols
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        header_set = set(header)
        HEADER_SET = set(HEADER)
        if len(header_set) < len(header):
            seen = {}
            duplicates = []
            for x in header:
                if x not in seen:
                    seen[x] = 1
                else:
                    if seen[x] > 1:
                        duplicates.append(x)
                    seen[x] += 1
            print(
                "ERROR: Please check samplesheet header:"
                f"\n-> Duplicate columns found: {','.join(duplicates)}"
            )
        if len(header_set) < len(HEADER_SET):
            expected_columns = sorted(list(HEADER_SET.difference(header_set)))
            unexpected_columns = sorted(list(header_set.difference(HEADER_SET)))
            if expected_columns:
                expected_columns_str = f"\n-> Expected {','.join(expected_columns)} columns but did not see them"
            else:
                expected_columns_str = ""
            if unexpected_columns:
                unexpected_columns_str = f"\n-> Saw {','.join(unexpected_columns)} columns but but did not expect them, ignoring"
            else:
                unexpected_columns_str = ""
            print(
                f"ERROR: Please check samplesheet header:"
                + expected_columns_str
                + unexpected_columns_str
                + f"\nParameters: skip_star={skip_star}, skip_classinput={skip_classinput}, skip_glm={skip_glm}"
            )
            if expected_columns:
                sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            ## Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    f"Invalid number of columns (minimum = {len(HEADER)})!",
                    "Line",
                    line,
                )

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    f"Invalid number of populated columns (minimum = {MIN_COLS})!",
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample_id, strandedness, concatenation_id, *remaining_cols = lspl[
                : len(HEADER)
            ]
            if sample_id:
                if sample_id.find(" ") != -1:
                    print_error("Sample entry contains spaces!", "Line", line)
            else:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for col in remaining_cols:
                if col:
                    if col.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if do_alignment:
                        # If doing alignment, then the remaining columns can only be fastqs
                        if not col.endswith(".fastq.gz") and not col.endswith(".fq.gz"):
                            print_error(
                                "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                                "Line",
                                line,
                            )

            ## Check strandedness
            strandednesses = ["unstranded", "forward", "reverse"]
            if strandedness:
                if strandedness not in strandednesses:
                    print_error(
                        f"Strandedness must be one of '{', '.join(strandednesses)}'!",
                        "Line",
                        line,
                    )
            else:
                print_error(
                    f"Strandedness has not been specified! Must be one of {', '.join(strandednesses)}.",
                    "Line",
                    line,
                )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2, strandedness]
            if (
                sample_id and do_alignment and remaining_cols[0] and remaining_cols[1]
            ):  ## Paired-end short reads
                sample_info = [
                    "0",
                    strandedness,
                    concatenation_id,
                    remaining_cols[0],
                    remaining_cols[1],
                ]
            elif (
                sample_id and do_alignment and remaining_cols[0]
            ):  ## Single-end short reads
                sample_info = [
                    "1",
                    strandedness,
                    concatenation_id,
                    remaining_cols[0],
                    "",
                ]
            elif sample_id:
                sample_info = ["1", strandedness, concatenation_id] + remaining_cols
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, strandedness ]]}
            if sample_id not in sample_mapping_dict:
                sample_mapping_dict[sample_id] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample_id]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample_id].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample_id",
                        "single_end",
                        "strandedness",
                        "concatenation_id",
                    ]
                    + additional_cols
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                if not all(
                    x[0] == sample_mapping_dict[sample][0][0]
                    for x in sample_mapping_dict[sample]
                ):
                    print_error(
                        f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                ## Check that multiple runs of the same sample are of the same strandedness
                if not all(
                    x[-1] == sample_mapping_dict[sample][0][-1]
                    for x in sample_mapping_dict[sample]
                ):
                    print_error(
                        f"Multiple runs of a sample must have the same strandedness!",
                        "Sample",
                        sample,
                    )

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join([f"{sample}_T{idx+1}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(
        args.FILE_IN, args.FILE_OUT, args.skip_star, args.skip_classinput, args.skip_glm
    )


if __name__ == "__main__":
    sys.exit(main())
