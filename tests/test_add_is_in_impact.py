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
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT

if __name__ == "__main__":
    from tools import run_command, load_mutations, run_cwl
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT

cwl_file = os.path.join(CWL_DIR, 'maf_filter.cwl')

class TestMafFilter(unittest.TestCase):
    def test_filter_a_maf_file(self):
        """
        Test that a filtered maf file comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        impact_file = os.path.join(, "") ##### <== add IMPACT file name here

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
                "IMPACT_filename": {
                      "class": "File",
                      "path": impact_file
                    },
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            with open(output_json['analysis_mutations_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 12518)

            # validate output mutation file contents
            comments, mutations = load_mutations(output_json['analysis_mutations_file']['path'])
            expected_comments, expected_mutations = load_mutations(os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "analyst_file.txt"))

            for mutation in expected_mutations:
                self.assertTrue(mutation in mutations)

            self.assertEqual(len(mutations), len(expected_mutations))

            comments, mutations = load_mutations(output_json['cbio_mutation_data_file']['path'])
            expected_comments, expected_mutations = load_mutations(os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Sample1", "portal_file.txt"))

            for mutation in expected_mutations:
                self.assertTrue(mutation in mutations)

            self.assertEqual(len(mutations), len(expected_mutations))

            expected_output = {
                'analysis_mutations_file': {
                    'location': 'file://' + os.path.join(output_dir, "Proj_08390_G.muts.maf"),
                    'basename': "Proj_08390_G.muts.maf",
                    'class': 'File',
                    'checksum': 'sha1$24421ab8d1a39a71f48eecbb0dd167d5d9f5c529',
                    'size': 28079,
                    'path': os.path.join(output_dir, "Proj_08390_G.muts.maf")
                    }
                }
            self.assertDictEqual(output_json, expected_output)


if __name__ == "__main__":
    unittest.main()
