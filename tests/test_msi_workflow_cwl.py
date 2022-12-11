#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the MSI analysis workflow cwl which uses multiple input samples
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
    MICROSATELLITES_LIST,
    OFile,
)
sys.path.pop(0)

class TestMsiWorkflow(PlutoTestCase):
    cwl_file = CWLFile('msi_workflow.cwl')

    # def setUp(self):
    #     # initialize the tmpdir
    #     super().setUp()

    def test_msi_workflow_demo1(self):
        """
        Test case for running the MSI workflow on single sample
        """
        normal_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")

        self.input = {
            "threads": "16",
            "microsatellites_file": {
                "class": "File",
                "path": MICROSATELLITES_LIST
            },

            "pairs": [
                {
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                }
            ],

            "normal_bam_files": [
                { "path": normal_bam,"class": "File" }
            ],

            "tumor_bam_files": [
                { "path": tumor_bam, "class": "File" }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [{
                "pair_id": "Sample1.Sample2",
                "tumor_id": "Sample1",
                "normal_id": "Sample2",
                "msi_tsv": OFile(name = 'Sample1.Sample2.msi.tsv', hash = "92576a9be4d6a36c67b26d16fdc4134b0d1b9cd9", size = 54, dir = output_dir)
                }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'Sample1.Sample2.msi.tsv')
        lines = self.read_table(output_file)

        expected_lines = [
            ['MSI_SCORE', 'MSI_STATUS', 'SAMPLE_ID'],
            ['20.90', 'Instable', 'Sample1']
            ]
        self.assertEqual(lines, expected_lines)

    # @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_msi_workflow1(self):
        """
        Test case for running the MSI workflow on multiple samples
        """
        # data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
        normal_bam2 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample35.rg.md.abra.printreads.bam")
        tumor_bam2  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample36.rg.md.abra.printreads.bam")

        self.input = {
            "microsatellites_file": {
                "class": "File",
                "path": MICROSATELLITES_LIST
            },

            "pairs": [
                {
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                },
                {
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
            ],

            "normal_bam_files": [
                { "path": normal_bam,"class": "File" },
                { "path": normal_bam2,"class": "File" }

            ],

            "tumor_bam_files": [
                { "path": tumor_bam, "class": "File" },
                { "path": tumor_bam2, "class": "File" }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        "pairs":[
            {
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "msi_tsv": OFile(name = 'Sample1-T.Sample1-N.msi.tsv', hash = "bc132f6ab9b779d7cba51e7ddfa82af724134f03", size = 56, dir = output_dir)
            },
            {
                "pair_id": "Sample2-T.Sample2-N",
                "tumor_id": "Sample2-T",
                "normal_id": "Sample2-N",
                "msi_tsv": OFile(name = 'Sample2-T.Sample2-N.msi.tsv', hash = "11fcf9459010aa5ea06e62e72155807c9723d45a", size = 56, dir = output_dir)
                }
            ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'Sample2-T.Sample2-N.msi.tsv')
        lines = self.read_table(output_file)
        expected_lines = [
        ['MSI_SCORE', 'MSI_STATUS', 'SAMPLE_ID'],
        ['40.14', 'Instable', 'Sample2-T']
        ]
        self.assertEqual(lines, expected_lines)

        output_file = os.path.join(output_dir,'Sample1-T.Sample1-N.msi.tsv')
        lines = self.read_table(output_file)
        expected_lines = [
        ['MSI_SCORE', 'MSI_STATUS', 'SAMPLE_ID'],
        ['21.97', 'Instable', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)




    # @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_msi_workflow2(self):
        """
        Test case for running the MSI workflow on single sample
        """
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")

        self.input = {
            "microsatellites_file": {
                "class": "File",
                "path": MICROSATELLITES_LIST
            },

            "pairs": [
                {
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
            ],

            "normal_bam_files": [
                { "path": normal_bam,"class": "File" }
            ],

            "tumor_bam_files": [
                { "path": tumor_bam, "class": "File" }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        "pairs":[{
            "pair_id": "Sample1-T.Sample1-N",
            "tumor_id": "Sample1-T",
            "normal_id": "Sample1-N",
            "msi_tsv": OFile(name = 'Sample1-T.Sample1-N.msi.tsv', hash = "bc132f6ab9b779d7cba51e7ddfa82af724134f03", size = 56, dir = output_dir)
            }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'Sample1-T.Sample1-N.msi.tsv')
        lines = self.read_table(output_file)

        expected_lines = [
            ['MSI_SCORE', 'MSI_STATUS', 'SAMPLE_ID'],
            ['21.97', 'Instable', 'Sample1-T']
            ]
        self.assertEqual(lines, expected_lines)


if __name__ == "__main__":
    unittest.main()
