#!/usr/bin/env python3
# parse the output of an Argos run and convert it into input format for pluto CLI operator
# need to gather up the tumor/normal files and output a formatted .tsv
import sys
import os
from collections import OrderedDict

def main():
    """
    """
    args = sys.argv[1:]
    sample_pairing_file = args[0]
    output_dir = args[1]
    bam_dir = os.path.join(output_dir, "bam")
    maf_dir = os.path.join(output_dir, "maf")

    samples = []

    # load the input samples and find their files
    with open(sample_pairing_file) as f:
        for line in f:
            parts = line.split()
            normal_id = parts[0]
            tumor_id = parts[1]
            pair_id = "{}.{}".format(tumor_id, normal_id)
            pair_maf = os.path.join(maf_dir, "{}.muts.maf".format(pair_id))
            # we will want these files later
            # tumor_bam_file = os.path.join(bam_dir, "{}.rg.md.abra.printreads.bam".format(tumor_id))
            # normal_bam_file = os.path.join(bam_dir, "{}.rg.md.abra.printreads.bam".format(normal_id))
            # "snp_pileup"
            sample = OrderedDict([
                ("tumor_id", tumor_id),
                ("normal_id", normal_id),
                ("pair_id", pair_id),
                ("pair_maf", pair_maf)
                ])
            samples.append(sample)

    # save the output ; write it to stdout so we can just pipe it to final file
    fout = sys.stdout
    # print headers
    header = '\t'.join(samples[0].keys()) + '\n'
    fout.write(header)
    # print samples
    for sample in samples:
        line = '\t'.join(sample.values()) + '\n'
        fout.write(line)

if __name__ == '__main__':
    main()
