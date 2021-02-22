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
sys.path.pop(0)

class TestMSIWorkflow(TmpDirTestCase):
    cwl_file = CWLFile('msi_workflow.cwl')
    def setUp(self):
        # initialize the tmpdir
        super().setUp()
        self.data_clinical_lines = [
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#STRING', 'STRING', 'NUMBER',],
        ['#1', '1', '1'],
        ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['Sample1-T', 'Patient1', '108'],
        ['Sample1-N', 'Patient2', '58'],
        ['Sample2-T', 'Patient3', '502'],
        ['Sample2-N', 'Patient4', '56'],
        ]

        self.data_clinical_file = write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        self.normal_bam_1= "inputs/foo_normal.rg.md.bam"
        self.tumor_bam_1 = "inputs/foo_tumor.rg.md.bam"
        self.normal_bam_2= "inputs/foo_normal_2.rg.md.bam"
        self.tumor_bam_2 = "inputs/foo_tumor.rg.md.bam"

        self.microsatellites_list = "/work/ci/resources/request_files/msisensor/b37_known_somatic_microsatellites.list"

    def test_msi_workflow1(self):
        """
        Test case for running the MSI workflow on multiple files
        """
        self.maxDiff = None
        input_json = {
            "data_clinical_file": {
                  "class": "File",
                  "path": self.data_clinical_file
                },
            "microsatellites_list": {
                  "class": "File",
                  "path": self.microsatellites_list
                },
            "pairs": [
                {
                    "normal_bam": {
                        "path": self.normal_bam_1,
                        "class": "File"
                    },
                    "tumor_bam": {
                        "path": self.tumor_bam_1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                },
                {
                    "normal_bam": {
                        "path": self.normal_bam_2,
                        "class": "File"
                    },
                    "tumor_bam": {
                        "path": self.tumor_bam_2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }

        output_json, output_dir = run_cwl(
            testcase = self,
            tmpdir = self.tmpdir,
            input_json = input_json,
            cwl_file = cwl_file,
            print_command = False,
            )

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$09243e8ed17d6c1be5307603d92de4c4f1fa6e5f',
                'size': 277,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            lines = [ l.strip().split() for l in fin ]

        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
            ['#STRING', 'STRING', 'NUMBER'],
            ['#1', '1', '1'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['Sample1-T', 'Patient1', '108', '0.000000006'],
            ['Sample1-N', 'Patient2', '58'], # tailing empty value gets stripped off
            ['Sample2-T', 'Patient3', '502', '0.000000005'],
            ['Sample2-N', 'Patient4', '56'] # tailing empty value gets stripped off
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
