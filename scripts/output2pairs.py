#!/usr/bin/env python3
"""
Script to parse the output of an Argos pipeline run
and convert it into a pairs samplesheet for use with various pipeline CLI operators

Usage
$ ./output2pairs.py /path/to/argos_output /path/to/argos_output/sample_pairing.txt pairs.tsv

Example
$ ./output2pairs.py /juno/work/ci/helix_filters_01/fixtures/demo/ /juno/work/ci/helix_filters_01/fixtures/demo/inputs/demo_sample_pairing.txt pairs.tsv
"""
import sys
import os
from collections import OrderedDict

def main():
    """
    Main control function for the script
    """
    args = sys.argv[1:]
    argos_output_dir = args[0]
    sample_pairing_file = args[1] # os.path.join(argos_output_dir, "sample_pairing.txt")
    try:
        output_file = args[2]
        fout = open(output_file, "w")
    except IndexError:
        fout = sys.stdout

    bam_dir = os.path.join(argos_output_dir, "bam")
    maf_dir = os.path.join(argos_output_dir, "maf")
    snp_pileup_dir = os.path.join(argos_output_dir, "snp-pileup") # "snp_pileup" ; this one comes from Facets output

    samples = []

    # load the input samples and find their files
    with open(sample_pairing_file) as f:
        for line in f:
            parts = line.split()
            normal_id = parts[0]
            tumor_id = parts[1]
            pair_id = "{}.{}".format(tumor_id, normal_id)
            pair_maf = os.path.join(maf_dir, "{}.muts.maf".format(pair_id))
            tumor_bam_file = os.path.join(bam_dir, "{}.rg.md.abra.printreads.bam".format(tumor_id))
            normal_bam_file = os.path.join(bam_dir, "{}.rg.md.abra.printreads.bam".format(normal_id))

            sample = OrderedDict([
                ("tumor_id", tumor_id),
                ("normal_id", normal_id),
                ("pair_id", pair_id),
                ("pair_maf", pair_maf),
                ("tumor_bam", tumor_bam_file),
                ("normal_bam", normal_bam_file),
                ])
            if os.path.exists(snp_pileup_dir):
                sample["snp_pileup"] = os.path.join(snp_pileup_dir, "{}.snp_pileup.gz".format(pair_id))
            samples.append(sample)

    # print headers
    header = '\t'.join(samples[0].keys()) + '\n'
    fout.write(header)
    # print samples
    for sample in samples:
        line = '\t'.join(sample.values()) + '\n'
        fout.write(line)
    fout.close()

if __name__ == '__main__':
    main()
