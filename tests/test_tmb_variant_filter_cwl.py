#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile
)
sys.path.pop(0)

# handle for errors arising from python3 -m unittest ...
try:
    import fixtures_tmb as fxt
except ModuleNotFoundError:
    sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
    import fixtures_tmb as fxt
    sys.path.pop(0)


class TestTMBVariantFilter(PlutoTestCase):
    cwl_file = CWLFile('tmb_variant_filter.cwl')

    def test_tmb_filter(self):
        """
        Test cases for filtering variants for TMB
        """
        maf_lines = self.dicts2lines(dict_list = fxt.rows1, comment_list = fxt.comments)
        input_maf = self.write_table(self.tmpdir, filename = "input.maf", lines = maf_lines)
        output_file = os.path.join(self.tmpdir, "output.txt")

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "output_filename":  'output.maf',
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$b019b5b3c6aba861371c135fe47520c969fab5ae',
                'size': 352,
                'path':  os.path.join(output_dir,'output.maf')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(os.path.join(output_dir,'output.maf'))

        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
        {'t_af': '0.50', 't_depth': '550', 'Hugo_Symbol': 'EGFR', 'Start_Position': '1', 'Consequence': 'missense_variant'},
          {'Consequence': 'missense_variant', 'Hugo_Symbol': 'EGFR', 'Start_Position': '1', 't_af': '0.51', 't_depth': '90'},
        {'t_af': '0.45', 't_depth': '590', 'Hugo_Symbol': 'EGFR', 'Start_Position': '1', 'Consequence': 'splice_region_variant'},
        {'t_af': '0.45', 't_depth': '590', 'Hugo_Symbol': 'TERT', 'Start_Position': '1295340', 'Consequence': 'splice_region_variant'},
        {'t_af': '0.45', 't_depth': '590', 'Hugo_Symbol': 'TERT', 'Start_Position': '1295339', 'Consequence': 'splice_region_variant'},
        {'t_af': '0.45', 't_depth': '590', 'Hugo_Symbol': 'TERT', 'Start_Position': '1295341', 'Consequence': 'splice_region_variant'},
        {'t_af': '0.45', 't_depth': '590', 'Hugo_Symbol': 'TERT', 'Start_Position': '1295339', 'Consequence': 'synonymous_variant'}
        ]

        self.assertEqual(comments, expected_comments)
        self.assertEqual(mutations, expected_mutations)

if __name__ == "__main__":
    unittest.main()
