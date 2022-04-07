#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the mergebed cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestMergeBed(PlutoTestCase):
    cwl_file = CWLFile('mergebed.cwl')

    def test_merge_bed(self):
        lines1 = [
            '1\t118166398\t118166398\n',
            '5\t35876484\t35876484\n',
            '17\t11998935\t11998935\n',
            '20\t62321135\t62321135\n'
        ]

        lines2 = [
            '5\t35876484\t35876484\n', # in common
            '7\t116418998\t116418998\n',
            '7\t151845367\t151845367\n',
            '17\t11998935\t11998935\n', # in common
        ]
        bed1 = os.path.join(self.tmpdir, "1.bed")
        bed2 = os.path.join(self.tmpdir, "2.bed")
        with open(bed1, "w") as fout:
            for line in lines1:
                fout.write(line)
        with open(bed2, "w") as fout:
            for line in lines2:
                fout.write(line)

        self.input = {
            "bed_files": [
                {"class": "File", "path": bed1},
                {"class": "File", "path": bed2}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'merged.bed'),
                'basename': 'merged.bed',
                'class': 'File',
                'checksum': 'sha1$638f9f3ef43802b8e372c3cef3848b16f2af1c66',
                'size': 149,
                'path':  os.path.join(output_dir,'merged.bed')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'merged.bed')
        with open(output_file) as fin:
            lines = [ line for line in fin ]

        expected_lines = [
            '1\t118166398\t118166398\n',
            '5\t35876484\t35876484\n',
            '17\t11998935\t11998935\n',
            '20\t62321135\t62321135\n',
            '7\t116418998\t116418998\n',
            '7\t151845367\t151845367\n',
            '17\t11998935\t11998935\n'
        ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
