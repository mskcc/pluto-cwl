#!/bin/bash
# script for copying all the recognized input files from a CWL input.json file
# to subdirs in the local directory
# designed for usage with pluto-cwl / helix_filters_01 CWL input.json files
# also make sure the JSON file has indentation enabled otherwise this wont work
set -eu
input_file="${1}"

maf_dir="${PWD}/maf"
bam_dir="${PWD}/bam"
snp_dir="${PWD}/snp_pileup"

mkdir -p "${maf_dir}"
mkdir -p "${bam_dir}"
mkdir -p "${snp_dir}"

while read line; do
    case $line in
        *.maf|*.svs.pass.vep.portal.txt)
        cp -va "$line" "${maf_dir}/"
        ;;

        *.bam|*.bai)
        cp -va "$line" "${bam_dir}/"
        ;;

        *.dat.gz)
        echo "snp_pileup"
        cp -va "$line" "${snp_dir}/"
        ;;

        *)
        echo "no match $line"
        ;;
    esac
done < <(grep "path" "${input_file}" | cut -d ':' -f2 | tr -d '"' | tr -d ',')
