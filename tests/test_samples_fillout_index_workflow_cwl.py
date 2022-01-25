#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_index_workflow cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, md5_obj
sys.path.pop(0)

class TestSamplesFilloutIndex(PlutoTestCase):
    cwl_file = CWLFile('samples_fillout_index_workflow.cwl')

    def setUp(self):
        super().setUp()

        self.sample1_bam = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample1.bam")
        self.sample1_maf = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        self.sample4_bam = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample4.bam")
        self.sample4_maf = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.muts.maf")

        self.sample24_bam = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample24.bam")
        self.sample24_maf = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample24.Sample23.muts.maf")

    def test_run_fillout_workflow(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache for samples fillout workflow it breaks on split_vcf_to_mafs

        self.input = {
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
            "sample_ids": ["Sample1"],
            "bam_files": [
                { "class": "File", "path": self.sample1_bam }
            ],
            "maf_files": [
                { "class": "File", "path": self.sample1_maf }
            ],
            "unindexed_sample_ids": ["Sample4", "Sample24"],
            "unindexed_bam_files": [
                { "class": "File", "path": self.sample4_bam },
                { "class": "File", "path": self.sample24_bam }
            ],
            "unindexed_maf_files": [
                { "class": "File", "path": self.sample4_maf },
                { "class": "File", "path": self.sample24_maf }
            ],
            "fillout_output_fname": 'output.maf'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$d8d63a0aca2da20d2d15e26fcf64fd6295eda05e',
                # 'size': 25838928,
                'path':  os.path.join(output_dir,'output.maf')
                }
            }
        # NOTE: for some reason, this file keeps coming out with different annotations for 'splice_acceptor_variant' or `splice_donor_variant`
        # this keeps changing the byte size and checksum so need to remove those here for now
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']
        comments, mutations = self.load_mutations(output_file)

        self.assertEqual(len(mutations), 23742)

        # Need to remove these fields because they are inconsistent on the output maf file;
        for mut in mutations:
            mut.pop('all_effects')
            mut.pop('Consequence')
            mut.pop('Variant_Classification')

        hash = md5_obj(mutations)
        expected_hash = 'f153c68bc79a6f28e261ec04f51b2111'
        self.assertEqual(hash, expected_hash)


if __name__ == "__main__":
    unittest.main()
