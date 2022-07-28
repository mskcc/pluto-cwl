#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

comments = [
['# comment 1'],
['# comment 2']
]
row1 = {
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
}
row2 = {
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'splice_region_variant,synonymous_variant' # exclude due to synonymous_variant
}
row3 = { # this one should pass filter
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}
row4 = {
't_af': '0.01', # exclude due to low AF
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}
row5 = {
't_af': '0.51',
't_depth': '90', # exclude due to low coverage
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}
row6 = { # this one should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'splice_region_variant'
}
row7 = { # this one should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295340', # good value; is_TERT_promoter = True
'Consequence': 'splice_region_variant'
}
row8 = { # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295339', # good value; is_TERT_promoter = True
'Consequence': 'splice_region_variant'
}
row9 = { # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295341', # bad value; is_TERT_promoter = False
'Consequence': 'splice_region_variant' # include anyway because its not synonymous_variant
}
row10 ={ # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295339', # good value; is_TERT_promoter = True
'Consequence': 'synonymous_variant' # include even though its synonymous_variant
}

class TestTMBCWL(PlutoTestCase):
    cwl_file = CWLFile('tmb.cwl')
    def test_tmb_workflow1(self):
        """
        Test case for the TMB analysis workflow
        """
        self.maxDiff = None

        maf_rows = [ row1, row2, row3, row4, row5, row6, row7, row8, row9, row10 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = comments)
        input_maf = self.write_table(self.tmpdir, filename = "input.maf", lines = maf_lines)
        # output_file = os.path.join(self.tmpdir, "output.txt")

        self.input = {
            "mutations_file": {
                  "class": "File",
                  "path": input_maf
                },
            "assay_coverage":  '1000',
            "sample_id": "Sample1",
            "normal_id": "Sample1-N"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'Sample1.Sample1-N.tmb.tsv'),
                'basename': 'Sample1.Sample1-N.tmb.tsv',
                'class': 'File',
                'checksum': 'sha1$7059b4e8d0e2fb18c6b3b7a6270ec945f4307857',
                'size': 38,
                'path':  os.path.join(output_dir,'Sample1.Sample1-N.tmb.tsv')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'Sample1.Sample1-N.tmb.tsv')
        with open(output_file) as fin:
            lines = [ l.strip().split() for l in fin ]
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['7000.0', 'Sample1']
        ]
        self.assertEqual(lines, expected_lines)

    def test_tmb_workflow2(self):
        # A pooled Normal gives a NA result output
        maf_rows = [ row1, row2, row3, row4, row5, row6, row7, row8, row9, row10 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = comments)
        input_maf = self.write_table(self.tmpdir, filename = "input.maf", lines = maf_lines)

        self.input = {
            "mutations_file": {
                  "class": "File",
                  "path": input_maf
                },
            "assay_coverage":  '1000',
            "sample_id": "Sample1",
            "normal_id": "Sample1-PooledNormal"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'Sample1.Sample1-PooledNormal.tmb.tsv'),
                'basename': 'Sample1.Sample1-PooledNormal.tmb.tsv',
                'class': 'File',
                'checksum': 'sha1$f822b3bbbe1ef2281f1caee3c5efec04c7740b41',
                'size': 34,
                'path':  os.path.join(output_dir,'Sample1.Sample1-PooledNormal.tmb.tsv')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'Sample1.Sample1-PooledNormal.tmb.tsv')
        with open(output_file) as fin:
            lines = [ l.strip().split() for l in fin ]
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['NA', 'Sample1']
        ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
