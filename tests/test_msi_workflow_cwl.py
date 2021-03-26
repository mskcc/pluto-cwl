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
from pluto.settings import DATA_SETS #, MICROSATELLITES_LIST
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

        # self.tmpdir="/work/ci/vurals/pluto-cwl/tmp"

        self.data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        self.normal_bam = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        self.tumor_bam  = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")

        self.normal_bam2 = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample35.rg.md.abra.printreads.bam")
        self.tumor_bam2  = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample36.rg.md.abra.printreads.bam")

        self.microsatellites_file = '/work/ci/resources/request_files/msisensor/microsatellites.list' #'/work/ci/vurals/pluto-cwl/b37_known_somatic_microsatellites.list'



    def test_msi_workflow1(self):
        """
        Test case for running the MSI workflow on multiple samples
        """
        self.runner_args['debug']= False #True
        self.preserve = False #True

        print("----->",self.tmpdir)

        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": self.data_clinical_file
                },

            "microsatellites_file": {
                "class": "File",
                "path": self.microsatellites_file #MICROSATELLITES_LIST
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
                { "path": self.normal_bam,"class": "File" },
                { "path": self.normal_bam2,"class": "File" }

            ],

            "tumor_bam_files": [
                { "path": self.tumor_bam, "class": "File" },
                { "path": self.tumor_bam2, "class": "File" }
            ]
        }

        output_json, output_dir = self.run_cwl()
        # print ( output_json, output_dir)

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
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']

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



    def test_msi_workflow2(self):
        """
        Test case for running the MSI workflow on single sample
        """
        self.runner_args['debug']= False #True
        self.preserve = False#True

        # print("----->",self.tmpdir)

        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": self.data_clinical_file
                },

            "microsatellites_file": {
                "class": "File",
                "path": self.microsatellites_file #MICROSATELLITES_LIST
            },

            "pairs": [
                {
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
            ],

            "normal_bam_files": [
                { "path": self.normal_bam,"class": "File" }
            ],

            "tumor_bam_files": [
                { "path": self.tumor_bam, "class": "File" }
            ]
        }

        output_json, output_dir = self.run_cwl()

        # print ( output_json, output_dir)

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
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
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
