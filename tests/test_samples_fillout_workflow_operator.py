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
    def test_operator_inputs(self):
        """
        Test that the CWL inputs are created correctly based on Operator args
        """
        self.maxDiff = None
        maf1 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.maf")
        maf2 = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample2.maf")
        bam1 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")
        bam2 = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        ref_fasta = self.DATA_SETS['demo']['REF_FASTA']
        lines = [
        ['sample_id', 'bam_file', 'maf_file'],
        ['Sample1', bam1, maf1],
        ['Sample2', bam2, maf2],
        ]
        samplesheet = self.write_table(self.tmpdir, filename = "samples.tsv", lines = lines)

        operator = SamplesFillout(samplesheet = samplesheet, ref_fasta = ref_fasta)
        input = operator.input
        # print(input)
        expected_input = {
            'ref_fasta': {'class': 'File', 'path': ref_fasta},
            'samples': [
                {
                    'sample_id': 'Sample1',
                    'maf_file': {'class': 'File', 'path': maf1},
                    'bam_file': {'class': 'File', 'path': bam1}
                    },
                {
                    'sample_id': 'Sample2',
                    'maf_file': {'class': 'File', 'path': maf2},
                    'bam_file': {'class': 'File', 'path': bam2}
                }
            ]
        }
        self.assertEqual(input, expected_input)



if __name__ == "__main__":
    unittest.main()
