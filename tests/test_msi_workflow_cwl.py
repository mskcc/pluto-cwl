#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the MSI analysis workflow cwl which uses multiple input samples
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
from pluto.settings import ENABLE_LARGE_TESTS, MICROSATELLITES_LIST
sys.path.pop(0)

class TestMsiWorkflow(PlutoTestCase):
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
        ['Sample6-T', 'Patient4', '57'],
        ['Sample7-N', 'Patient4', '58'],
        ]

    def test_msi_workflow_demo1(self):
        """
        Test case for running the MSI workflow on single sample
        """
        data_clinical_lines = [
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#STRING', 'STRING', 'NUMBER',],
        ['#1', '1', '1'],
        ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['Sample1', 'Patient1', '108'],
        ['Sample2', 'Patient2', '58'],
        ['Sample6-T', 'Patient4', '57'], # not included as input
        ['Sample7-N', 'Patient4', '58'],
        ]
        data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = data_clinical_lines)
        normal_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")

        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": data_clinical_file
                },

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
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$749f2c93dc1cc0f43bcb7e92ec5151b58f842dc3',
                'size': 357,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'data_clinical_sample.txt')
        lines = self.read_table(output_file)

        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#STRING', 'STRING', 'NUMBER', 'NUMBER','STRING'],
            ['#1', '1', '1', '0','0'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'MSI_SCORE','MSI_STATUS'],
            ['Sample1', 'Patient1',   '108',             '20.90',    'Instable'  ],
            ['Sample2', 'Patient2',   '58',              'NA',       'NA'        ],
            ['Sample6-T', 'Patient4',   '57',              'NA',       'NA'        ],
            ['Sample7-N', 'Patient4',   '58',              'NA',       'NA'        ]
            ]
        self.assertEqual(lines, expected_lines)

    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_msi_workflow1(self):
        """
        Test case for running the MSI workflow on multiple samples
        """
        data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
        normal_bam2 = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample35.rg.md.abra.printreads.bam")
        tumor_bam2  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample36.rg.md.abra.printreads.bam")

        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": data_clinical_file
                },

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
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$0a7be413f57372afd5142294467e3e91f926e9df',
                'size': 427,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'data_clinical_sample.txt')

        lines = self.read_table(output_file)

        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#STRING', 'STRING', 'NUMBER', 'NUMBER','STRING'],
            ['#1', '1', '1', '0','0'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'MSI_SCORE','MSI_STATUS'],
            ['Sample1-T', 'Patient1',   '108',             '21.97',    'Instable'  ],
            ['Sample1-N', 'Patient2',   '58',              'NA',       'NA'        ],
            ['Sample2-T', 'Patient3',   '502',             '40.14',    'Instable'  ],
            ['Sample2-N', 'Patient4',   '56',              'NA',       'NA'        ],
            ['Sample6-T', 'Patient4',   '57',              'NA',       'NA'        ],
            ['Sample7-N', 'Patient4',   '58',              'NA',       'NA'        ]
            ]
        self.assertEqual(lines, expected_lines)


    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_msi_workflow2(self):
        """
        Test case for running the MSI workflow on single sample
        """
        data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        tumor_bam  = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")

        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": data_clinical_file
                },

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
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$1de2dc1ab3d49022b451e23485ce5e9ddd7daad9',
                'size': 418,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'data_clinical_sample.txt')
        lines = self.read_table(output_file)

        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_MSI_SCORE','CMO_MSI_STATUS'],
            ['#STRING', 'STRING', 'NUMBER', 'NUMBER','STRING'],
            ['#1', '1', '1', '0','0'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'MSI_SCORE','MSI_STATUS'],
            ['Sample1-T', 'Patient1',   '108',             '21.97',    'Instable'  ],
            ['Sample1-N', 'Patient2',   '58',              'NA',       'NA'        ],
            ['Sample2-T', 'Patient3',   '502',             'NA',       'NA'        ],
            ['Sample2-N', 'Patient4',   '56',              'NA',       'NA'        ],
            ['Sample6-T', 'Patient4',   '57',              'NA',       'NA'        ],
            ['Sample7-N', 'Patient4',   '58',              'NA',       'NA'        ]
            ]
        self.assertEqual(lines, expected_lines)


if __name__ == "__main__":
    unittest.main()
