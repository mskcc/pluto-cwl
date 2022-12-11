#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the annotate-maf-wrapper.cwl file
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile, 
    TableReader,
    ENABLE_LARGE_TESTS,
    OFile,
    ODir
)
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
            'failed': False,
            'output_file': OFile(name = 'Sample1_hisens.ccf.maf', hash = '8cd487056bd86177d19d3dd0fe072747d31fd9b7', size = 60230, dir = output_dir),
            'stderr_txt': OFile(name = 'annotate_maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'annotate_maf_stdout.txt', dir = output_dir),
        }
        strip_related_keys = [
        ('basename', 'annotate_maf_stderr.txt', ['size', 'checksum']),
        ('basename', 'annotate_maf_stdout.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        path = os.path.join(output_dir, 'Sample1_hisens.ccf.maf')
        table_reader = TableReader(path)
        comments = table_reader.comment_lines
        fieldnames = table_reader.get_fieldnames()
        records = [ rec for rec in table_reader.read() ]
        self.assertEqual(len(records), 41)


    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_facets_annotation_wrapper(self):
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
            'failed': False,
            'output_file': OFile(name = 'Sample1_hisens.ccf.maf', hash = '7e478a8a44d27735f26e368989c672ed6ef5d52a', size = 19217199, dir = output_dir),
            'stderr_txt': OFile(name = 'annotate_maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'annotate_maf_stdout.txt', dir = output_dir),
        }
        strip_related_keys = [
        ('basename', 'annotate_maf_stderr.txt', ['size', 'checksum']),
        ('basename', 'annotate_maf_stdout.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

if __name__ == "__main__":
    unittest.main()
