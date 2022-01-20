#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestTMBVariantFilter(PlutoTestCase):
    cwl_file = CWLFile('tmb_variant_filter.cwl')

    def test_tmb_filter(self):
        """
        Test cases for filtering variants for TMB
        """
        self.maxDiff = None
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
        maf_rows = [ row1, row2, row3, row4, row5, row6, row7, row8, row9, row10 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = comments)
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
