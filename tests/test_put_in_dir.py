#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the put_in_dir.cwl
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile,
    OFile, 
    ODir,
)
sys.path.pop(0)


class TestPutInDir(PlutoTestCase):
    cwl_file = CWLFile('put_in_dir.cwl')

    def test_put_two_files_in_dir(self):
        """
        Test that two files are put in the dir correctly
        """
        file1 = self.mkstemp(prefix = "1.")
        file2 = self.mkstemp(prefix = "2.")

        # create input data
        self.input = {
            "output_directory_name": "foo",
            "files": [
                {
                  "class": "File",
                  "path": file1
                },
                {
                  "class": "File",
                  "path": file2
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "directory": ODir(name = "foo", dir = output_dir, items = [
                OFile(name = os.path.basename(file1), size = 0, hash = 'da39a3ee5e6b4b0d3255bfef95601890afd80709'),
                OFile(name = os.path.basename(file2), size = 0, hash = 'da39a3ee5e6b4b0d3255bfef95601890afd80709'),
            ])
        }

        self.assertCWLDictEqual(output_json, expected_output)

    def test_put_one_file1_in_dir(self):
        """
        Test that one file is put in the dir correctly
        """
        file1 = self.mkstemp(prefix = "1.")
        self.input = {
            "output_directory_name": "foo",
            "files": [
                {
                  "class": "File",
                  "path": file1
                },
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "directory": ODir(name = 'foo', dir = output_dir, items = [
                OFile(name = os.path.basename(file1), size = 0, hash = 'da39a3ee5e6b4b0d3255bfef95601890afd80709')
                ])
        }

        self.assertCWLDictEqual(output_json, expected_output)




