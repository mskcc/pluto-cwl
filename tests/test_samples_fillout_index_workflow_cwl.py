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
from pluto.serializer import OFile, ODir
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
            'output_file': OFile(name = 'output.maf', dir = output_dir),
        }

        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        # instead of checksum and size, count the number of mutations and take a checksum on the mutation contents
        comments, mutations = self.load_mutations(output_file, strip = True)
        self.assertEqual(len(mutations), 23742)
        hash = md5_obj(mutations)
        expected_hash = '53ee95f4e07084992362bba4917d1e82'
        self.assertEqual(hash, expected_hash)


if __name__ == "__main__":
    unittest.main()
