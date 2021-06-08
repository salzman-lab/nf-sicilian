#!/usr/bin/env python
import argparse
import logging
import os

import boto3
import pandas as pd

# Create a logger
logging.basicConfig(format="%(name)s - %(asctime)s %(levelname)s: %(message)s")
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


# cribbed from https://github.com/chanzuckerberg/s3mi/blob/master/scripts/s3mi
def s3_bucket_and_key(s3_uri, require_prefix=False):
    """
    Return the bucket and key name as a two-element list, key here could be
    file or folder under the bucket
    """
    prefix = "s3://"
    if require_prefix:
        assert s3_uri.startswith(prefix), "The path must start with s3://"
    elif not s3_uri.startswith(prefix):
        prefix = ""

    return s3_uri[len(prefix) :].split("/", 1)


def make_samplesheet_from_s3_path(s3_path, suffix="fastq.gz", strandedness="forward"):
    bucket_name, key = s3_bucket_and_key(s3_path)
    bucket = s3.Bucket(bucket_name)
    file_objects = bucket.objects.filter(Prefix=key).all()
    logger.info(f"Finding files ending with {suffix} in {s3_path}")

    fastqs = pd.Series(
        [
            os.path.join(f"s3://{bucket_name}", file_obj.key)
            for file_obj in file_objects
            if file_obj.key.endswith(suffix)
        ],
        name="fastq",
    ).to_frame()

    # Remove all index (I1, I2) reads
    fastqs = fastqs.loc[~fastqs["fastq"].str.contains("_I[12]_")]

    # Get the R1/R2 read numbers
    fastqs["read_number"] = fastqs["fastq"].str.extract("_(R[12])_")
    fastqs["fastq_number"] = fastqs.read_number.str.replace("R", "fastq_")

    fastqs["sample_id"] = fastqs.apply(
        # Split on `_S\d+` is a biohub specific thing
        lambda x: os.path.basename(x["fastq"]).split("_" + x["read_number"] + "_")[0],
        axis=1,
    )

    # Split on `_S\d+` is a biohub specific thing
    # This converts the per-lane sample IDs to a single concatenation ID
    # MACA_24m_M_BLADDER_58_S5_L001 --> MACA_24m_M_BLADDER_58
    # MACA_24m_M_BLADDER_58_S5_L002 --> MACA_24m_M_BLADDER_58
    fastqs["concatenation_id"] = fastqs["sample_id"].str.split("_S\d+").str[0]

    samplesheet = fastqs.pivot(
        index="sample_id", columns="fastq_number", values="fastq"
    )
    # Convert the sample id to the concatenation id
    sample_id_to_concatenation_id = dict(
        zip(fastqs["sample_id"], fastqs["concatenation_id"])
    )
    samplesheet["concatenation_id"] = samplesheet.index.map(
        sample_id_to_concatenation_id
    )
    samplesheet = samplesheet.reset_index()
    # samplesheet = samplesheet.rename(columns={"lane": "concatenation_id"})
    # samplesheet["single_end"] = samplesheet["fastq_2"].isnull()

    # Remove the column name so it doesn't get written
    samplesheet.columns.name = None

    # Add empty required columns
    samplesheet["strandedness"] = strandedness

    # Reorder columns to match
    samplesheet = samplesheet[
        [
            "sample_id",
            "strandedness",
            "concatenation_id",
            "fastq_1",
            "fastq_2",
        ]
    ]
    return samplesheet


def main():
    parser = argparse.ArgumentParser(
        description="""Create an input csv of sample R1, R2 fastqs. Must have AWS credentials set up already."""
    )
    parser.add_argument("s3_path", type=str, help="S3 path to parse")
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="samplesheet.csv",
        type=str,
        help="Comma-separated variable file of samplesheet data",
    )
    method = parser.add_mutually_exclusive_group(required=True)
    method.add_argument(
        "--tenx",
        action="store_true",
        help="Specify 10x input, which is stranded",
    )
    method.add_argument(
        "--smartseq2",
        action="store_true",
        help="Specify SmartSeq2 input, which is unstranded",
    )
    parser.add_argument(
        "-s",
        "--suffix",
        default="fastq.gz",
        type=str,
        help="Suffix of files on s3 to look for",
    )

    args = parser.parse_args()

    global s3
    s3 = boto3.resource("s3")

    strandedness = "forward" if args.tenx else "unstranded"
    samplesheet = make_samplesheet_from_s3_path(args.s3_path, args.suffix, strandedness)
    logger.info(f"Writing samplesheet to {args.output}")
    samplesheet.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
