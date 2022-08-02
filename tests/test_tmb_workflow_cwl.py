#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis workflow cwl which uses multiple input samples
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
from pluto.serializer import OFile
sys.path.pop(0)

# handle for errors arising from python3 -m unittest ...
try:
    import fixtures_tmb as fxt
except ModuleNotFoundError:
    sys.path.insert(0, THIS_DIR)
    import fixtures_tmb as fxt
    sys.path.pop(0)


class TestTmbWorkflow(PlutoTestCase):
    cwl_file = CWLFile('tmb_workflow.cwl')

    def setUp(self):
        super().setUp() # initialize the tmpdir

        # PASS: 6
        # FAIL: 4
        self.maf_lines1 = self.dicts2lines(dict_list = fxt.rows1, comment_list = fxt.comments)

        # PASS: 5
        # FAIL: 5
        self.maf_lines2 = self.dicts2lines(dict_list = fxt.rows2, comment_list = fxt.comments)
        self.maf1 = self.write_table(self.tmpdir, filename = "input1.maf", lines = self.maf_lines1)
        self.maf2 = self.write_table(self.tmpdir, filename = "input2.maf", lines = self.maf_lines2)

    def test_tmb_workflow1(self):
        """
        Test case for running the TMB workflow on multiple files
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": self.maf1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                },
                {
                    "pair_maf": {
                        "path": self.maf2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        "pairs":[
            {
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 40, hash = "d6a57cfb5e3001697875e5b5bfae206e0f7f2310", dir = output_dir)
            },
            {
                "pair_id": "Sample2-T.Sample2-N",
                "tumor_id": "Sample2-T",
                "normal_id": "Sample2-N",
                "tmb_maf": OFile(name = "Sample2-T.Sample2-N.tmb.maf", size = 308, hash = "e6cb790e887d606e2762dd73e73fec28ecfff22b", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample2-T.Sample2-N.tmb.tsv", size = 40, hash = "851b303d4237654de54fdfdc60d41ab996f4380d", dir = output_dir)
            }
            ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['7000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)

        lines = self.read_table(expected_output["pairs"][1]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['6000.0', 'Sample2-T']
        ]
        self.assertEqual(lines, expected_lines)


    def test_tmb_workflow1_2(self):
        # one pair uses a Pooled Normal
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": self.maf1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-PooledNormal",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-PooledNormal"
                },
                {
                    "pair_maf": {
                        "path": self.maf2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "pairs":[
                {
                    "pair_id": "Sample1-T.Sample1-PooledNormal",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-PooledNormal",
                    "tmb_maf": OFile(name = "Sample1-T.Sample1-PooledNormal.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                    "tmb_tsv": OFile(name = "Sample1-T.Sample1-PooledNormal.tmb.tsv", size = 36, hash = "84ba7be1f81eba42d9845ccb18cae9d4b3d6a21f", dir = output_dir)
                },
                {
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N",
                    "tmb_maf": OFile(name = "Sample2-T.Sample2-N.tmb.maf", size = 308, hash = "e6cb790e887d606e2762dd73e73fec28ecfff22b", dir = output_dir),
                    "tmb_tsv": OFile(name = "Sample2-T.Sample2-N.tmb.tsv", size = 40, hash = "851b303d4237654de54fdfdc60d41ab996f4380d", dir = output_dir)
                }
                ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # !!!!!!!!
        # Sample1 gets value NA because its paired with a pooled normal !!
        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['NA', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)
        # !!!!!!!!

        lines = self.read_table(expected_output["pairs"][1]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['6000.0', 'Sample2-T']
        ]
        self.assertEqual(lines, expected_lines)

    def test_tmb_workflow2(self):
        """
        Test case for using a single input pair maf
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {"path": self.maf1, "class": "File"},
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [{
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 40, hash = "d6a57cfb5e3001697875e5b5bfae206e0f7f2310", dir = output_dir)
            }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['7000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)


    def test_tmb_workflow3(self):
        """
        Test case with a single real maf file
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [{
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 519440, hash = "809c3c1ac3bb750aebf22ee2f95a5ebafd41e98f", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 42, hash = "8156b9e7a0602ddd7710f002ef9385237a82c5d0", dir = output_dir)
            }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['475000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)


if __name__ == "__main__":
    unittest.main()
