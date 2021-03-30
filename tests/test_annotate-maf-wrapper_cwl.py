#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the annotate-maf-wrapper.cwl file
"""
import os
import sys
# import json
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
from pluto.settings import ENABLE_LARGE_TESTS
sys.path.pop(0)

class TestAnnotateMafWrapperCWL(PlutoTestCase):
    cwl_file = CWLFile('annotate-maf-wrapper.cwl')

    def test_annotate_demo1(self):
        """
        Test case for running Facets maf annotation on a smaller demo maf file
        """
        input_maf = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        input_rds = os.path.join(self.DATA_SETS['demo']['FACETS_DIR'], "Sample1_hisens.rds")
        self.input = {
            "maf_file": {
                "path": input_maf,
                "class": "File"
            },
            "facets_rds": {
                "path": input_rds,
                "class": "File"
            },
            "output_filename": "Sample1_hisens.ccf.maf"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed_txt': None,
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, 'Sample1_hisens.ccf.maf'),
                'basename': 'Sample1_hisens.ccf.maf',
                'class': 'File',
                'checksum': 'sha1$8cd487056bd86177d19d3dd0fe072747d31fd9b7',
                'size': 60230,
                'path': os.path.join(output_dir, 'Sample1_hisens.ccf.maf')
            },
            'stderr_txt': {
                'basename': 'annotate_maf_stderr.txt',
                'checksum': 'sha1$b87a1dfbeeaf0f3addd5e5efdddbd2e7bbc55f03',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir,'annotate_maf_stderr.txt'),
                'path': os.path.join(output_dir,'annotate_maf_stderr.txt'),
                'size': 37
            },
           'stdout_txt': {
                'basename': 'annotate_maf_stdout.txt',
                'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir,'annotate_maf_stdout.txt'),
                'path': os.path.join(output_dir,'annotate_maf_stdout.txt'),
                'size': 0
            }
        }
        self.maxDiff = None
        self.assertDictEqual(output_json, expected_output)

        path = os.path.join(output_dir, 'Sample1_hisens.ccf.maf')
        table_reader = TableReader(path)
        comments = table_reader.comment_lines
        fieldnames = table_reader.get_fieldnames()
        records = [ rec for rec in table_reader.read() ]
        self.assertEqual(len(records), 41)

        stdout_txt = os.path.join(output_dir,'annotate_maf_stdout.txt')
        stderr_txt = os.path.join(output_dir,'annotate_maf_stderr.txt')
        # with open(stdout_txt) as f:
        #     lines = [ l for l in f ]
        # print(lines)
        # with open(stderr_txt) as f:
        #     lines = [ l for l in f ]
        # self.assertEqual(lines, [])


    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_facets_wrapper(self):
        """
        Test case for running Facets maf annotation
        """
        input_maf = os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        input_rds = os.path.join(self.DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample1_hisens.rds")
        self.input = {
            "maf_file": {
                "path": input_maf,
                "class": "File"
            },
            "facets_rds": {
                "path": input_rds,
                "class": "File"
            },
            "output_filename": "Sample1_hisens.ccf.maf"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed_txt': None,
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, 'Sample1_hisens.ccf.maf'),
                'basename': 'Sample1_hisens.ccf.maf',
                'class': 'File',
                'checksum': 'sha1$7e478a8a44d27735f26e368989c672ed6ef5d52a',
                'size': 19217199,
                'path': os.path.join(output_dir, 'Sample1_hisens.ccf.maf')
            },
            'stderr_txt': {
                'basename': 'annotate_maf_stderr.txt',
                'checksum': 'sha1$2e672f99c23a2d827c1d33e06377870cdd9c8090',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir,'annotate_maf_stderr.txt'),
                'path': os.path.join(output_dir,'annotate_maf_stderr.txt'),
                'size': 105
            },
           'stdout_txt': {
                'basename': 'annotate_maf_stdout.txt',
                'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir,'annotate_maf_stdout.txt'),
                'path': os.path.join(output_dir,'annotate_maf_stdout.txt'),
                'size': 0
            }
        }
        self.maxDiff = None
        self.assertDictEqual(output_json, expected_output)
        # with open(output_json['stdout_txt']["path"]) as f:
        #     for line in f:
        #         print(line)
        #
        # with open(output_json['stderr_txt']["path"]) as f:
        #     for line in f:
        #         print(line)

if __name__ == "__main__":
    unittest.main()
