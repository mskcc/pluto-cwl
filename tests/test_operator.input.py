#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis workflow Operator
"""
import os
import sys
import unittest
import json

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile 
from operators import input
sys.path.pop(0)

class TestInputGeneratePairs(PlutoTestCase):
    def setUp(self):
        # initialize the tmpdir
        super().setUp()
        self.maf1 = os.path.join(self.tmpdir, "1.maf")
        self.maf2 = os.path.join(self.tmpdir, "2.maf")
        self.pileup1 = os.path.join(self.tmpdir, "1.pileup")
        self.pileup2 = os.path.join(self.tmpdir, "2.pileup")
        self.pairs_dicts = [
            {
                'tumor_id': 'Sample1-T',
                'normal_id': 'Sample1-N',
                'pair_id': 'Sample1-T.Sample1-N',
                'pair_maf': '',
                'snp_pileup': ''
            },
            {
                'tumor_id': 'Sample2-T',
                'normal_id': 'Sample2-N',
                'pair_id': 'Sample2-T.Sample2-N',
                'pair_maf': '',
                'snp_pileup': ''
            }
        ]
        self.pair_template = {
                "pair_maf": {
                    "path": None,
                    "class": "File"
                },
                "pair_id": None,
                "tumor_id": None,
                "normal_id": None
            }

    def test_generate_pairs1(self):
        """
        Test case for generating a pairs list from a single input set
        """
        # need to make copies of self lists and dicts
        pairs_dicts = [ self.pairs_dicts[0] ]
        pairs_dicts[0] = { **pairs_dicts[0] }
        pairs_dicts[0]['pair_maf'] = self.maf1
        pairs_lines = self.dicts2lines(dict_list = pairs_dicts, comment_list = [])
        pairs_file = self.write_table(self.tmpdir, filename = "pairs.tsv", lines = pairs_lines)
        pairs = input.generate_pairs(pairs_file, self.pair_template)
        expected_pairs = [
            {
                'pair_maf': {
                    'path': self.maf1,
                    'class': 'File'
                    },
                'pair_id': 'Sample1-T.Sample1-N',
                'tumor_id': 'Sample1-T',
                'normal_id': 'Sample1-N'
            }
        ]
        self.assertEqual(pairs, expected_pairs)

    def test_generate_pairs2(self):
        """
        Test case for generating pairs list from two input pairs sets
        """
        pairs_dicts = [ *self.pairs_dicts ]
        for i, d in enumerate(pairs_dicts):
            pairs_dicts[i] = { **d }
        pairs_dicts[0]['pair_maf'] = self.maf1
        pairs_dicts[1]['pair_maf'] = self.maf2
        pairs_lines = self.dicts2lines(dict_list = pairs_dicts, comment_list = [])
        pairs_file = self.write_table(self.tmpdir, filename = "pairs.tsv", lines = pairs_lines)
        pairs = input.generate_pairs(pairs_file, self.pair_template)
        expected_pairs = [
            {
                'pair_maf': {
                    'path': self.maf1,
                    'class': 'File'
                    },
                'pair_id': 'Sample1-T.Sample1-N',
                'tumor_id': 'Sample1-T',
                'normal_id': 'Sample1-N'
            },
            {
                'pair_maf': {
                    'path': self.maf2,
                    'class': 'File'
                    },
                'pair_id': 'Sample2-T.Sample2-N',
                'tumor_id': 'Sample2-T',
                'normal_id': 'Sample2-N'
            }
        ]
        self.assertEqual(pairs, expected_pairs)

class TestInputGenerateInput(PlutoTestCase):
    def test_generate_input1(self):
        """
        Test case for generating some input data
        """
        # this will be a File input
        data_clinical_file = os.path.join(self.tmpdir, "data_clinical.txt"),

        # these paths will be read from the input file
        path1 = os.path.join(self.tmpdir, '1.txt')
        path2 = os.path.join(self.tmpdir, '2.txt')
        files_list_file = os.path.join(self.tmpdir, "some_files.txt")
        with open(files_list_file, "w") as fout:
            fout.write(path1 + '\n')
            fout.write(path2 + '\n')

        args = {
            'data_clinical_file': data_clinical_file,
            'assay_coverage': '1000',
            'is_impact': "True",
            'input_files': files_list_file
        }

        input_data = input.generate_input(
            args,
            File_keys = ['data_clinical_file'],
            bool_keys = ['is_impact'],
            array_File_keys = ['input_files']
        )
        expected_input_data = {
            'data_clinical_file': {
                'class': 'File',
                'path': data_clinical_file
            },
            'assay_coverage': '1000',
            'is_impact': True,
            'input_files': [
                {'class': 'File', 'path': path1},
                {'class': 'File', 'path': path2}
            ]
        }
        self.assertEqual(input_data, expected_input_data)

if __name__ == "__main__":
    unittest.main()
