#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the maf_col_filter.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile # run_command, load_mutations, run_cwl, write_table,
from pluto.settings import DATA_SETS
from pluto.serializer import OFile
sys.path.pop(0)

# handle for errors arising from python3 -m unittest ...
try:
    import fixtures_maf_col_filter_cwl as fxt
except ModuleNotFoundError:
    sys.path.insert(0, THIS_DIR)
    import fixtures_maf_col_filter_cwl as fxt
    sys.path.pop(0)


class TestMafColFilter(PlutoTestCase):
    cwl_file = CWLFile('maf_col_filter.cwl')

    def test_filter_maf_file_cols(self):
        """
        Filter columns in a tiny demo maf file
        """
        maf_lines = [
            ['# comment 1'], # keep the comments
            ['# comment 2'],
            ['Hugo_Symbol', 'foo_value'], # foo_value column should be removed in output
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        # run the script in a temporary directory
        input_maf_file = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_maf_file
                },
            "output_filename": "output.maf"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = 'output.maf', hash = "e55f7bdaa146f37b48d6c920ed27184e394ef1e6", size = 46, dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # validate number of lines output
        with open(expected_output['output_file']['path']) as fin:
            output_maf_lines = len(fin.readlines())
        self.assertEqual(output_maf_lines, 5)

        # validate file contents
        self.assertMutFileContains(
            filepath = expected_output['output_file']['path'],
            expected_comments = ['# comment 1', '# comment 2'],
            expected_mutations = [{'Hugo_Symbol': 'SUFU'}, {'Hugo_Symbol': 'GOT1'}],
            identical = True)

    def test_filter_maf_file_cols_full(self):
        """
        Test col filter on a full sized dataset


        ======================================================================
        FAIL: test_filter_maf_file_cols_full (__main__.TestMafColFilter)
        ----------------------------------------------------------------------
        Traceback (most recent call last):
          File "tests/test_maf_col_filter_cwl.py", line 99, in test_filter_maf_file_cols_full
            self.assertMutHeadersContain(filepath = os.path.join(output_dir, 'output.maf'), expected_headers = fxt.cols_to_keep)
          File "/juno/work/ci/kellys5/projects/pluto-dev/pluto-cwl4/pluto/tools.py", line 1246, in assertMutHeadersContain
            self.assertEqual(len(missingWanted), 0, message, *args, **kwargs)
        AssertionError: 1 != 0 : Expected columns {'Amino_Acid_Change'} missing from mutation file
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "output_filename": "output.maf"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = 'output.maf', hash = "a2f5b9f1533fd443b41561ca718ffca62ab45f36", size = 2710681, dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # validate number of lines output
        with open(os.path.join(output_dir, 'output.maf')) as fin:
            output_maf_lines = len(fin.readlines())
        self.assertEqual(output_maf_lines, 12518)

        # validate file contents
        self.assertNumMutations(os.path.join(output_dir, 'output.maf'), 12514)

        self.assertMutHeadersContain(filepath = os.path.join(output_dir, 'output.maf'), expected_headers = fxt.cols_to_keep)
        # for key in mutations[0].keys():
        #     self.assertTrue(key in fxt.cols_to_keep)
        #
        # # make sure there are fewer than or equal to the number of columns in new output as there are entries to keep
        # self.assertTrue( len(mutations[0].keys()) <= len(fxt.cols_to_keep) )


if __name__ == "__main__":
    unittest.main()
