#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the samples_fillout_workflow.cwl Operator
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import TableReader, PlutoTestCase
from operators.samples_fillout import SamplesFillout
sys.path.pop(0)

class TestSamplesFilloutWorkflowOperator(PlutoTestCase):
    def setUp(self):
        super().setUp()
        # make a demo samplesheet from the demo files
        self.maf1 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        self.maf2 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.muts.maf")
        self.bam1 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample1.bam")
        self.bam2 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample4.bam")
        self.ref_fasta = self.DATA_SETS['demo']['REF_FASTA']
        lines = [
        ['sample_id', 'bam_file', 'maf_file'],
        ['Sample1', self.bam1, self.maf1],
        ['Sample4', self.bam2, self.maf2],
        ]
        self.samplesheet = self.write_table(self.tmpdir, filename = "samples.tsv", lines = lines)

    def test_operator_inputs(self):
        """
        Test that the CWL inputs are created correctly based on Operator args
        """
        self.maxDiff = None
        operator = SamplesFillout(samplesheet = self.samplesheet, ref_fasta = self.ref_fasta)
        input = operator.input

        # this is the old input format that works with CWL version 1.1
        # expected_input = {
        #     'ref_fasta': {'class': 'File', 'path': ref_fasta},
        #     'samples': [
        #         {
        #             'sample_id': 'Sample1',
        #             'maf_file': {'class': 'File', 'path': maf1},
        #             'bam_file': {'class': 'File', 'path': bam1}
        #             },
        #         {
        #             'sample_id': 'Sample2',
        #             'maf_file': {'class': 'File', 'path': maf2},
        #             'bam_file': {'class': 'File', 'path': bam2}
        #         }
        #     ]
        # }
        expected_input = {
        'sample_ids': ['Sample1', 'Sample4'],
        'bam_files': [{'class': 'File', 'path': self.bam1}, {'class': 'File', 'path': self.bam2}],
        'maf_files': [{'class': 'File', 'path': self.maf1}, {'class': 'File', 'path': self.maf2}],
        'ref_fasta': {'class': 'File', 'path': self.ref_fasta}
        }
        self.assertEqual(input, expected_input)

    def test_operator_runs1(self):
        """
        Test case that the Operator can run the CWL correctly
        """
        operator = SamplesFillout(samplesheet = self.samplesheet, ref_fasta = self.ref_fasta,
            dir=self.tmpdir, output_dir=self.tmpdir, verbose = False)
        output_json, output_dir, output_json_file = operator.run()
        expected_output = {
            'output_file':
                {
                'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$3fec625cac6d500a0bc074753ce686a5685229f4',
                'size': 17724,
                'path': os.path.join(output_dir, 'output.maf')
                }
            }
        self.maxDiff = None
        self.assertDictEqual(output_json, expected_output)

        path = output_json['output_file']['path']
        comments, mutations = self.load_mutations(path)
        self.assertEqual(len(mutations), 106)



if __name__ == "__main__":
    unittest.main()
