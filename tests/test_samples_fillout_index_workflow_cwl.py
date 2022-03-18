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
        # self.preserve = True
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache for samples fillout workflow it breaks on split_vcf_to_mafs
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "Sample1-N",
                    "prefilter": True,
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": self.sample1_maf },
                    "bam_file": { "class": "File", "path": self.sample1_bam }
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "DMP-MATCHED-NORMAL",
                    "prefilter": False,
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": self.sample4_maf },
                    "bam_file": { "class": "File", "path": self.sample4_bam }
                },
                {
                    "sample_id": "Sample24",
                    "normal_id": "DMP-UNMATCHED-NORMAL",
                    "prefilter": False,
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": self.sample24_maf },
                    "bam_file": { "class": "File", "path": self.sample24_bam }
                },
            ],
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()
        output_file = os.path.join(output_dir,'output.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_file,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$d8d63a0aca2da20d2d15e26fcf64fd6295eda05e',
                # 'size': 25838928,
                'path':  output_file
                }
            }
        # NOTE: for some reason, this file keeps coming out with different annotations for 'splice_acceptor_variant' or `splice_donor_variant`
        # this keeps changing the byte size and checksum so need to remove those here for now
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        # instead of checksum and size, count the number of mutations and take a checksum on the mutation contents
        comments, mutations = self.load_mutations(output_file)
        self.assertEqual(len(mutations), 23742)

        # Need to remove these fields because they are inconsistent on the output maf file;
        for mut in mutations:
            mut.pop('all_effects')
            mut.pop('Consequence')
            mut.pop('Variant_Classification')

        hash = md5_obj(mutations)
        expected_hash = '40cb934e354cf985a9d0a03bc65747c8'
        self.assertEqual(hash, expected_hash)


if __name__ == "__main__":
    unittest.main()
