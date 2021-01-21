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

cwl_file = os.path.join(CWL_DIR, 'calc-tmb.cwl')

class TestCalcTMB(TmpDirTestCase):
    def test_calc_tmb_from_file(self):
        """
        Test case for calculating TMB by reading the number of variants from a file
        NOTE: variant entries are not parsed, only counted, so their contents do not matter here
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol', 'Chromosome'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        input_maf_file = write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_json = {
            "input_file": {
                  "class": "File",
                  "path": input_maf_file
                },
            "output_filename":  'output.txt',
            "genome_coverage": "1000",
            "normal_id": "Sample1-N"
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
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$2afff0129f070c317af9367958e0c9fd8713654d',
                'size': 12,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            result = next(fin).strip()
        expected_result = '0.000000005'
        self.assertEqual(result, expected_result)

    def test_calc_tmb_poolednormal1(self):
        """
        Test case for skipping TMB calculation if the normal is a pooled normal
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol', 'Chromosome'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        input_maf_file = write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_json = {
            "input_file": {
                  "class": "File",
                  "path": input_maf_file
                },
            "output_filename":  'output.txt',
            "genome_coverage": "1000",
            "normal_id": "Sample1PooledNormal"
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
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$7d9f637deb947080f063e9d70bdd6090968e1a7e',
                'size': 3,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            result = next(fin).strip()
        expected_result = 'NA'
        self.assertEqual(result, expected_result)

if __name__ == "__main__":
    unittest.main()
