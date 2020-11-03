#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the maf_filter.cwl
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command, load_mutations, run_cwl
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT, IMPACT_FILE

if __name__ == "__main__":
    from tools import run_command, load_mutations, run_cwl
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT, IMPACT_FILE

cwl_file = os.path.join(CWL_DIR, 'add_is_in_impact.cwl')

class TestMafFilter(unittest.TestCase):
    def test_filter_a_maf_file(self):
        """
        Test that a maf file with is_in_IMPACT column comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        impact_file = os.path.join(IMPACT_FILE)
        output_maf = 'output.maf'


        with open(input_maf) as fin:
            input_maf_lines = len(fin.readlines())
        self.assertEqual(input_maf_lines, 12518)

        with TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename":  output_maf,
                "IMPACT_filename": {
                      "class": "File",
                      "path": impact_file
                    },
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            with open(output_json['IMPACT_col_added_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 12518)

            # validate output mutation file contents
            input_comments,  input_mutations  = load_mutations(input_maf)
            output_comments, output_mutations = load_mutations(output_json['IMPACT_col_added_file']['path'])

            true_count=[row['is_in_impact'] for row in output_mutations].count('True')
            false_count=[row['is_in_impact'] for row in output_mutations].count('False')

            self.assertTrue(true_count == 8367)
            self.assertTrue(false_count == 4147)

            self.assertTrue(len(input_mutations[1])+1==len(output_mutations[1]))



if __name__ == "__main__":
    unittest.main()
