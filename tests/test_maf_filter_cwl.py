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

        with open(input_maf) as fin:
            input_maf_lines = len(fin.readlines())

        self.assertEqual(input_maf_lines, 12518)

        with TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "maf_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "argos_version_string": ARGOS_VERSION_STRING,
                "is_impact": True,
                "analysis_mutations_filename": "Proj_08390_G.muts.maf",
                "cbio_mutation_data_filename": 'data_mutations_extended.txt'
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            with open(output_json['analysis_mutations_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 27)

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
                    },
                'cbio_mutation_data_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_mutations_extended.txt'),
                    'basename': 'data_mutations_extended.txt',
                    'class': 'File',
                    'checksum': 'sha1$6131494536ce956d741c820378e7e2ce1c714403',
                    'size': 4534,
                    'path': os.path.join(output_dir, 'data_mutations_extended.txt')
                    },
                'rejected_file': {
                    'basename': 'rejected.muts.maf',
                    'checksum': 'sha1$a06789623715703c5006db6876ecb58b8498f938',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'rejected.muts.maf'),
                    'path': os.path.join(output_dir, 'rejected.muts.maf'),
                    'size': 18627626
                    }
                }
            self.assertDictEqual(output_json, expected_output)

    def test_maf_filter_argos_3_2_0(self):
        """
        Test the maf filter script results when used with argos_version_string 3.2.0
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        with open(input_maf) as fin:
            input_maf_lines = len(fin.readlines())

        self.assertEqual(input_maf_lines, 12518)

        with TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "maf_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "argos_version_string": "3.2.0",
                "is_impact": True,
                "analysis_mutations_filename": "Proj_08390_G.muts.maf",
                "cbio_mutation_data_filename": 'data_mutations_extended.txt'
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

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
                    'checksum': 'sha1$fd78842c9410e7e622dee270ec9c0e7628811f18',
                    'size': 28081,
                    'path': os.path.join(output_dir, "Proj_08390_G.muts.maf")
                    },
                'cbio_mutation_data_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_mutations_extended.txt'),
                    'basename': 'data_mutations_extended.txt',
                    'class': 'File',
                    'checksum': 'sha1$47e716eabbfda3408b2d9a08b9bb432b2cb8fce8',
                    'size': 4536,
                    'path': os.path.join(output_dir, 'data_mutations_extended.txt')
                    },
                'rejected_file': {
                    'basename': 'rejected.muts.maf',
                    'checksum': 'sha1$a06789623715703c5006db6876ecb58b8498f938',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'rejected.muts.maf'),
                    'path': os.path.join(output_dir, 'rejected.muts.maf'),
                    'size': 18627626
                    }
                }

            with open(output_json['analysis_mutations_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 27)

            self.assertDictEqual(output_json, expected_output)

    def test_filter_maf_file_impact_false(self):
        """
        Test that a filtered maf file comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        with open(input_maf) as fin:
            input_maf_lines = len(fin.readlines())

        self.assertEqual(input_maf_lines, 12518)

        with TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "maf_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "argos_version_string": ARGOS_VERSION_STRING,
                "is_impact": False,
                "analysis_mutations_filename": "Proj_08390_G.muts.maf",
                "cbio_mutation_data_filename": 'data_mutations_extended.txt'
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            with open(output_json['analysis_mutations_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 23)

            expected_output = {
                'analysis_mutations_file': {
                    'location': 'file://' + os.path.join(output_dir, "Proj_08390_G.muts.maf"),
                    'basename': "Proj_08390_G.muts.maf",
                    'class': 'File',
                    'checksum': 'sha1$9fb9d43c71e546750ddec6aea2313dda28547b3a',
                    'size': 24524,
                    'path': os.path.join(output_dir, "Proj_08390_G.muts.maf")
                    },
                'cbio_mutation_data_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_mutations_extended.txt'),
                    'basename': 'data_mutations_extended.txt',
                    'class': 'File',
                    'checksum': 'sha1$15ca06249511c32c32e058c246a757ec8df11d83',
                    'size': 3931,
                    'path': os.path.join(output_dir, 'data_mutations_extended.txt')
                    },
                'rejected_file': {
                    'basename': 'rejected.muts.maf',
                    'checksum': 'sha1$e7441703699e82cef500d9557bfcbd3464ce8eab',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'rejected.muts.maf'),
                    'path': os.path.join(output_dir, 'rejected.muts.maf'),
                    'size': 18790398
                    }
                }
            self.assertDictEqual(output_json, expected_output)

    def test_large_maf_file(self):
        """
        Test that a giant maf file with tons of variants gets filtered as expected
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "Proj_08390_G.muts.maf")

        with TemporaryDirectory() as tmpdir:
            output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "maf_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "argos_version_string": "2.x",
                "is_impact": True,
                "analysis_mutations_filename": "Proj_08390_G.muts.maf",
                "cbio_mutation_data_filename": 'data_mutations_extended.txt'
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

            with open(output_json['analysis_mutations_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 1664)

            with open(output_json['cbio_mutation_data_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 1141)

            # validate output mutation file contents
            comments, mutations = load_mutations(output_json['analysis_mutations_file']['path'])
            expected_comments, expected_mutations = load_mutations(os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "analyst_file.txt"))

            for mutation in expected_mutations:
                self.assertTrue(mutation in mutations)

            self.assertEqual(len(mutations), len(expected_mutations))

            comments, mutations = load_mutations(output_json['cbio_mutation_data_file']['path'])
            expected_comments, expected_mutations = load_mutations(os.path.join(DATA_SETS['Proj_08390_G']['MAF_FILTER_DIR'], "Proj_08390_G", "portal_file.txt"))

            for mutation in expected_mutations:
                self.assertTrue(mutation in mutations)

            self.assertEqual(len(mutations), len(expected_mutations))


if __name__ == "__main__":
    unittest.main()
