#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase
from pluto.serializer import OFile
sys.path.pop(0)


class TestUpdateCaseList(PlutoTestCase):
    cwl_file = 'updateCaseList.cwl'

    def setUp(self):
        super().setUp()
        case_list_str = """case_list_category: all_cases_in_study
stable_id: pi_123_all
case_list_name: All Tumors
case_list_description: All tumor samples
cancer_study_identifier: pi_123
case_list_ids: Sample1\tSample2"""
        self.input_file = os.path.join(self.tmpdir, "cases.txt")
        with open(self.input_file, "w") as fout:
            fout.write(case_list_str)

    def test_update_caselist_1(self):
        """
        Test simple case list update
        """
        self.input = {
            "case_list": {"class": "File", "path": self.input_file},
            "sample_ids": ["Sample3", "Sample4"],
            "output_filename": "cases_all.txt"
        }
        output_json, output_dir = self.run_cwl()

        output_file = os.path.join(output_dir, 'cases_all.txt')

        expected_output = {
            "output_file": OFile(name = 'cases_all.txt', size = 208, hash = '59aa4b6b6695b7adfd4493390edd4038808a018f', dir = output_dir)
        }

        self.assertCWLDictEqual(output_json, expected_output)

        with open(output_file) as fin:
            text = fin.read()

        expected_text = """case_list_category: all_cases_in_study
stable_id: pi_123_all
case_list_name: All Tumors
case_list_description: All tumor samples
cancer_study_identifier: pi_123
case_list_ids: Sample1\tSample2\tSample3\tSample4
"""
        self.assertEqual(text, expected_text)


    def test_update_caselist_2(self):
        """
        Test update with no files passed
        """
        self.input = {
            "case_list": {"class": "File", "path": self.input_file},
            "sample_ids": [],
            "output_filename": "cases_all.txt"
        }
        output_json, output_dir = self.run_cwl()

        output_file = os.path.join(output_dir, 'cases_all.txt')

        expected_output = {
            "output_file": OFile(name = 'cases_all.txt', size = 192, hash = 'f1ad64f51beac01759ae690b2f787fe3978e8882', dir = output_dir)
        }

        self.assertCWLDictEqual(output_json, expected_output)

        with open(output_file) as fin:
            text = fin.read()

        expected_text = """case_list_category: all_cases_in_study
stable_id: pi_123_all
case_list_name: All Tumors
case_list_description: All tumor samples
cancer_study_identifier: pi_123
case_list_ids: Sample1\tSample2
"""
        self.assertEqual(text, expected_text)




if __name__ == "__main__":
    unittest.main()
