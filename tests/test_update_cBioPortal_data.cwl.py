#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the udpate_cBioPortal_data.cwl
"""

import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.serializer import OFile
sys.path.pop(0)

# handle for errors arising from python3 -m unittest ...
try:
    import fixtures_cBioPortal as fxt
except ModuleNotFoundError:
    sys.path.insert(0, THIS_DIR)
    import fixtures_cBioPortal as fxt
    sys.path.pop(0)

class TestUpdate_cBioPortal_dataCWL(PlutoTestCase):
    cwl_file = CWLFile('update_cBioPortal_data.cwl')

    def setUp(self):
        super().setUp()
        self.maf_row1 = fxt.maf_row1
        self.maf_row2 = fxt.maf_row2
        self.maf_row3 = fxt.maf_row3
        self.maf_row4 = fxt.maf_row4
        self.maf_row5 = fxt.maf_row5
        self.facets_row1 = fxt.facets_row1
        self.facets_row2 = fxt.facets_row2
        self.facets_row3 = fxt.facets_row3
        self.facets_row4 = fxt.facets_row4
        self.facets_row5 = fxt.facets_row5

        self.demo_comments = fxt.demo_comments

        self.expected_row1 = fxt.expected_row1
        self.expected_row2 = fxt.expected_row2
        self.expected_row3 = fxt.expected_row3
        self.expected_row4 = fxt.expected_row4
        self.expected_row4_2 = fxt.expected_row4_2
        self.expected_row5 = fxt.expected_row5
        self.expected_row5_2 = fxt.expected_row5_2


    def test_add_clonality(self):
        """
        Test Update_cBioPortal_dataCW with tiny dataset; one sample in facets and portal maf
        """
        self.maxDiff = None
        # make sets of lines to write to tables
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = self.demo_comments)

        facets_rows = [ self.facets_row1, self.facets_row2, self.facets_row3 ]
        facets_lines = self.dicts2lines(dict_list = facets_rows, comment_list = self.demo_comments)

        input_maf = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_facets_file = self.write_table(self.tmpdir, filename = "facets.maf", lines = facets_lines)
        self.input = {
            "subcommand": "merge_mafs",
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "facets_maf":{
                  "class": "File",
                  "path": input_facets_file
                },
            "output_filename":  'output.maf',
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed': False,
            'output_file': OFile(name = 'output.maf', hash = '73273a9ac68aa9a1b6e8d8747f268ea009bf7f0d', size = 335, dir = output_dir),
            'stderr_txt': OFile(name = 'output.maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'output.maf_stdout.txt', dir = output_dir)
            }
        strip_related_keys = [
        ('basename', 'output.maf_stdout.txt', ['size', 'checksum']),
        ('basename', 'output.maf_stderr.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        mut_filepath = expected_output['output_file']['path']
        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
            self.expected_row1,
            self.expected_row2,
            self.expected_row3
        ]
        self.assertMutFileContains(mut_filepath, expected_comments, expected_mutations, identical = True)

    def test_add_clonality_two_samples(self):
        """
        Test case for adding clonality when more than one sample is in the facets and portal maf files
        """
        self.maxDiff = None
        # make sets of lines to write to tables
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4, self.maf_row5 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = self.demo_comments)

        facets_rows = [ self.facets_row1, self.facets_row2, self.facets_row3, self.facets_row4, self.facets_row5 ]
        facets_lines = self.dicts2lines(dict_list = facets_rows, comment_list = self.demo_comments)

        input_maf = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_facets_file = self.write_table(self.tmpdir, filename = "facets.maf", lines = facets_lines)
        self.input = {
            "subcommand": "merge_mafs",
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "facets_maf":{
                  "class": "File",
                  "path": input_facets_file
                },
            "output_filename":  'output.maf',
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed': False,
            'output_file': OFile(name = 'output.maf', hash = 'b266fe04825577ed0e9ca2fca9c861e8fb2f38ef', size = 450, dir = output_dir),
            'stderr_txt': OFile(name = 'output.maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'output.maf_stdout.txt', dir = output_dir),
            }
        strip_related_keys = [
        ('basename', 'output.maf_stdout.txt', ['size', 'checksum']),
        ('basename', 'output.maf_stderr.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        mut_filepath = expected_output['output_file']['path']
        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
            self.expected_row1,
            self.expected_row2,
            self.expected_row3,
            self.expected_row4,
            self.expected_row5
        ]
        self.assertMutFileContains(mut_filepath, expected_comments, expected_mutations, identical = True)


    def test_add_clonality_two_portal_one_facets(self):
        """
        Test case for adding clonality when two samples are in the portal maf and only one sample is in the facets maf
        """
        self.maxDiff = None
        # make sets of lines to write to tables
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4, self.maf_row5 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = self.demo_comments)

        facets_rows = [ self.facets_row1, self.facets_row2, self.facets_row3 ]
        facets_lines = self.dicts2lines(dict_list = facets_rows, comment_list = self.demo_comments)

        input_maf = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_facets_file = self.write_table(self.tmpdir, filename = "facets.maf", lines = facets_lines)
        self.input = {
            "subcommand": "merge_mafs",
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "facets_maf":{
                  "class": "File",
                  "path": input_facets_file
                },
            "output_filename":  'output.maf',
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed': False,
            'output_file': OFile(name = 'output.maf', hash = '79871da52b9f90b53f2bdba977416087c2c4729f', size = 450, dir = output_dir),
            'stderr_txt': OFile(name = 'output.maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'output.maf_stdout.txt', dir = output_dir),
            }
        strip_related_keys = [
        ('basename', 'output.maf_stdout.txt', ['size', 'checksum']),
        ('basename', 'output.maf_stderr.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        mut_filepath = expected_output['output_file']['path']
        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
            self.expected_row1,
            self.expected_row2,
            self.expected_row3,
            self.expected_row4_2,
            self.expected_row5_2
        ]
        self.assertMutFileContains(mut_filepath, expected_comments, expected_mutations, identical = True)


    def test_add_clonality_one_portal_two_facets(self):
        """
        Test case for adding clonality when one samples is in the portal maf and two samples are in the facets maf
        """
        self.maxDiff = None
        # make sets of lines to write to tables
        maf_rows = [ self.maf_row1, self.maf_row2, self.maf_row3 ]
        maf_lines = self.dicts2lines(dict_list = maf_rows, comment_list = self.demo_comments)

        facets_rows = [ self.facets_row1, self.facets_row2, self.facets_row3, self.facets_row4, self.facets_row5 ]
        facets_lines = self.dicts2lines(dict_list = facets_rows, comment_list = self.demo_comments)

        input_maf = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        input_facets_file = self.write_table(self.tmpdir, filename = "facets.maf", lines = facets_lines)
        self.input = {
            "subcommand": "merge_mafs",
            "input_file": {
                  "class": "File",
                  "path": input_maf
                },
            "facets_maf":{
                  "class": "File",
                  "path": input_facets_file
                },
            "output_filename":  'output.maf',
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'failed': False,
            'output_file': OFile(name = 'output.maf', hash = '73273a9ac68aa9a1b6e8d8747f268ea009bf7f0d', size = 335, dir = output_dir),
            'stderr_txt': OFile(name = 'output.maf_stderr.txt', dir = output_dir),
            'stdout_txt': OFile(name = 'output.maf_stdout.txt', dir = output_dir),
            }
        strip_related_keys = [
        ('basename', 'output.maf_stdout.txt', ['size', 'checksum']),
        ('basename', 'output.maf_stderr.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        mut_filepath = expected_output['output_file']['path']
        expected_comments = ['# comment 1', '# comment 2']
        expected_mutations = [
            self.expected_row1,
            self.expected_row2,
            self.expected_row3
        ]
        self.assertMutFileContains(mut_filepath, expected_comments, expected_mutations, identical = True)


if __name__ == "__main__":
    unittest.main()
