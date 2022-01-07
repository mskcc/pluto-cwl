#!/usr/bin/env python

# script for updating the "path" entries in a CWL input JSON used with the pluto-cwl / helix_filters workflows
# script will search for each File's basename in the provided search dir
# certain JSON keys' Files are expected to be in specific sub-dirs under the search dir
# prints output JSON to stdout

# Example search_dir
# $ ls -1 ../../demo/
# bam
# maf
# snp-pileup

import os
import sys
import json
import copy

def find(search_dir: str, file_name:str) -> str:
    """
    Find a file with the same basename in the search_dir
    Returns path to last file found with matching basename
    """
    basename = os.path.basename(file_name)
    path = None
    for dirpath, dirnames, filenames in os.walk(search_dir):
        for filename in filenames:
            if filename == basename:
                path = os.path.join(dirpath, filename)
    if path is None:
        print(">>> ERROR: path not found for: ", search_dir, basename)
        raise
    return(path)

# load input args
args = sys.argv[1:]
input_json = args[0]
search_dir = args[1]

# need the complete path for CWL inputs
search_dir = os.path.abspath(search_dir)

# expected file locations
maf_dir = os.path.join(search_dir, "maf")
bam_dir = os.path.join(search_dir, "bam")
snp_pileup_dir  = os.path.join(search_dir, "snp_pileup")

# load input JSON
with open(input_json) as f:
    data = json.load(f)

# update each recognized key based on known criteria
for key in data:

    if key == 'pairs':
        for pair in data[key]:
            old_path = pair['pair_maf']['path']
            new_path = find(maf_dir, old_path)
            pair['pair_maf']['path'] = new_path

            old_path = pair['snp_pileup']['path']
            new_path = find(snp_pileup_dir, old_path)
            pair['snp_pileup']['path'] = new_path

    elif key == 'tumor_bam_files':
        for bam in data[key]:
            old_path = bam['path']
            new_path = find(bam_dir, old_path)
            bam['path'] = new_path

            old_path = bam['secondaryFiles'][0]['path']
            new_path = find(bam_dir, old_path)
            bam['secondaryFiles'][0]['path'] = new_path

    elif key == 'normal_bam_files':
        for bam in data[key]:
            old_path = bam['path']
            new_path = find(bam_dir, old_path)
            bam['path'] = new_path

            old_path = bam['secondaryFiles'][0]['path']
            new_path = find(bam_dir, old_path)
            bam['secondaryFiles'][0]['path'] = new_path

    elif key == 'mutation_svs_maf_files':
        for svs in data[key]:
            old_path = svs['path']
            new_path = find(maf_dir, old_path)
            svs['path'] = new_path

    elif key == 'mutation_svs_txt_files':
        for svs in data[key]:
            old_path = svs['path']
            new_path = find(maf_dir, old_path)
            svs['path'] = new_path

print(json.dumps(data, indent = 4))
