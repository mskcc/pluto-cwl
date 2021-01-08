#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import unittest

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import TmpDirTestCase, load_mutations, run_cwl, write_table, dicts2lines
    from .settings import CWL_DIR

if __name__ == "__main__":
    from tools import TmpDirTestCase, load_mutations, run_cwl, write_table, dicts2lines
    from settings import CWL_DIR

cwl_file = os.path.join(CWL_DIR, 'tmb_variant_filter.cwl')

class TestTMBVariantFilter(TmpDirTestCase):
    def test_tmb_filter(self):
        comments = [
        ['# comment 1'],
        ['# comment 2']
        ]
        row1 = {
        't_af': '0.50',
        't_depth': '550',
        'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
        }
        row2 = {
        't_af': '0.50',
        't_depth': '550',
        'Consequence': 'splice_region_variant,synonymous_variant' # exclude due to synonymous_variant
        }
        row3 = { # this one should pass filter
        't_af': '0.50',
        't_depth': '550',
        'Consequence': 'missense_variant'
        }
        row4 = {
        't_af': '0.01', # exclude due to low AF
        't_depth': '550',
        'Consequence': 'missense_variant'
        }
        row5 = {
        't_af': '0.51',
        't_depth': '90', # exclude due to low coverage
        'Consequence': 'missense_variant'
        }
        row6 = { # this one should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Consequence': 'splice_region_variant'
        }
        maf_rows = [ row1, row2, row3, row4, row5, row6 ]
        maf_lines = dicts2lines(dict_list = maf_rows, comment_list = comments)
        input_maf = write_table(self.tmpdir, filename = "input.maf", lines = maf_lines)
        output_file = os.path.join(self.tmpdir, "output.txt")

        input_json = {
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "output_filename":  'output.maf',
            }
        output_json, output_dir = run_cwl(
            testcase = self,
            tmpdir = self.tmpdir,
            input_json = input_json,
            cwl_file = cwl_file,
            print_command = False,
            )

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$8a39725de0268c835923bc66d8363d66d6876676',
                'size': 106,
                'path':  os.path.join(output_dir,'output.maf')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        comments, mutations = load_mutations(output_json['output_file']['path'])
        
        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
        {'t_af': '0.50', 't_depth': '550', 'Consequence': 'missense_variant'},
        {'t_af': '0.45', 't_depth': '590', 'Consequence': 'splice_region_variant'}
        ]

        self.assertEqual(comments, expected_comments)
        self.assertEqual(mutations, expected_mutations)

if __name__ == "__main__":
    unittest.main()
