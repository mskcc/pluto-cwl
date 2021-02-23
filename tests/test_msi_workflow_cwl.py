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

        self.data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)
        self.normal_bam = "/work/ci/vurals/msi_test/test_stuff/foo_normal.rg.md.bam"
        self.tumor_bam = "/work/ci/vurals/msi_test/test_stuff/foo_tumor.rg.md.bam"


    def test_tmb_workflow1(self):
        """
        Test case for running the TMB workflow on multiple files
        """
        self.maxDiff = None
        self.input = {
            "data_clinical_file": {
                  "class": "File",
                  "path": self.data_clinical_file
                },

            "pairs": [
                {
                    "pair_normal_bam": {
                        "path": self.normal_bam,
                        "class": "File"
                    },
                    "pair_tumor_bam": {
                        "path": self.tumor_bam,
                        "class": "File"
                    },

                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }#,
                # {
                #     "pair_maf": {
                #         "path": self.maf2,
                #         "class": "File"
                #     },
                #     "pair_id": "Sample2-T.Sample2-N",
                #     "tumor_id": "Sample2-T",
                #     "normal_id": "Sample2-N"
                # }
                ]
            }

        output_json, output_dir = self.run_cwl()

        print ( output_json, output_dir)

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$e7975b7d9f34202750c4645d02f2f98796e22c74',
                'size': 363,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        # output_file = expected_output['output_file']['path']
        #
        # lines = self.read_table(output_file)
        #
        # expected_lines = [
        #     ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
        #     ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
        #     ['#STRING', 'STRING', 'NUMBER', 'NUMBER'],
        #     ['#1', '1', '1', '1'],
        #     ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
        #     ['Sample1-T', 'Patient1', '108', '0.000000006'],
        #     ['Sample1-N', 'Patient2', '58', 'NA'],
        #     ['Sample2-T', 'Patient3', '502', '0.000000005'],
        #     ['Sample2-N', 'Patient4', '56', 'NA'],
        #     ['Sample6-T', 'Patient4', '57', 'NA'],
        #     ['Sample7-N', 'Patient4', '58', 'NA']
        #     ]
        # self.assertEqual(lines, expected_lines)






    #
    # def test_tmb_workflow2(self):
    #     """
    #     Test case for using a single input pair maf
    #     """
    #     self.input = {
    #         "data_clinical_file": {
    #               "class": "File",
    #               "path": self.data_clinical_file
    #             },
    #         "assay_coverage":  '1000',
    #         "pairs": [
    #             {
    #                 "pair_maf": {
    #                     "path": self.maf1,
    #                     "class": "File"
    #                 },
    #                 "pair_id": "Sample1-T.Sample1-N",
    #                 "tumor_id": "Sample1-T",
    #                 "normal_id": "Sample1-N"
    #             }
    #             ]
    #         }
    #
    #     output_json, output_dir = self.run_cwl()
    #
    #     expected_output = {
    #         'output_file': {
    #             'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
    #             'basename': 'data_clinical_sample.txt',
    #             'class': 'File',
    #             'checksum': 'sha1$51dfe08e1f1da012880a6d55bd516400fe36cd5a',
    #             'size': 354,
    #             'path':  os.path.join(output_dir,'data_clinical_sample.txt')
    #             }
    #         }
    #     self.assertDictEqual(output_json, expected_output)
    #
    #     output_file = expected_output['output_file']['path']
    #
    #     lines = self.read_table(output_file)
    #
    #     expected_lines = [
    #         ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['#STRING', 'STRING', 'NUMBER', 'NUMBER'],
    #         ['#1', '1', '1', '1'],
    #         ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['Sample1-T', 'Patient1', '108', '0.000000006'],
    #         ['Sample1-N', 'Patient2', '58', 'NA'],
    #         ['Sample2-T', 'Patient3', '502', 'NA'],
    #         ['Sample2-N', 'Patient4', '56', 'NA'],
    #         ['Sample6-T', 'Patient4', '57', 'NA'],
    #         ['Sample7-N', 'Patient4', '58', 'NA']
    #         ]
    #     self.assertEqual(lines, expected_lines)
    #
    # def test_tmb_workflow3(self):
    #     """
    #     Test case with a single real maf file
    #     """
    #     self.input = {
    #         "data_clinical_file": {
    #               "class": "File",
    #               "path": self.data_clinical_file
    #             },
    #         "assay_coverage":  '1000',
    #         "pairs": [
    #             {
    #                 "pair_maf": {
    #                     "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
    #                     "class": "File"
    #                 },
    #                 "pair_id": "Sample1-T.Sample1-N",
    #                 "tumor_id": "Sample1-T",
    #                 "normal_id": "Sample1-N"
    #             }
    #             ]
    #         }
    #
    #     output_json, output_dir = self.run_cwl()
    #
    #     expected_output = {
    #         'output_file': {
    #             'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
    #             'basename': 'data_clinical_sample.txt',
    #             'class': 'File',
    #             'checksum': 'sha1$2ba1761afc6dba46d9b7d699ff29c52e3e04d5b5',
    #             'size': 354,
    #             'path':  os.path.join(output_dir,'data_clinical_sample.txt')
    #             }
    #         }
    #     self.assertDictEqual(output_json, expected_output)
    #
    #     output_file = expected_output['output_file']['path']
    #
    #     lines = self.read_table(output_file)
    #
    #     expected_lines = [
    #         ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['#STRING', 'STRING', 'NUMBER', 'NUMBER'],
    #         ['#1', '1', '1', '1'],
    #         ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
    #         ['Sample1-T', 'Patient1', '108', '0.000000013'],
    #         ['Sample1-N', 'Patient2', '58', 'NA'],
    #         ['Sample2-T', 'Patient3', '502', 'NA'],
    #         ['Sample2-N', 'Patient4', '56', 'NA'],
    #         ['Sample6-T', 'Patient4', '57', 'NA'],
    #         ['Sample7-N', 'Patient4', '58', 'NA']
    #         ]
    #     self.assertEqual(lines, expected_lines)
    #
    # def test_tmb_workflow4(self):
    #     """
    #     Test case with a real maf and real data clinical file
    #     """
    #     self.maxDiff = None
    #     data_clinical_file = os.path.join(self.DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
    #     self.input = {
    #         "data_clinical_file": {
    #               "class": "File",
    #               "path": data_clinical_file
    #             },
    #         "assay_coverage":  '1000',
    #         "pairs": [
    #             {
    #                 "pair_maf": {
    #                     "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
    #                     "class": "File"
    #                 },
    #                 "pair_id": "Sample1.Sample2",
    #                 "tumor_id": "Sample1",
    #                 "normal_id": "Sample2"
    #             }
    #             ]
    #         }
    #
    #     output_json, output_dir = self.run_cwl()
    #
    #     expected_output = {
    #         'output_file': {
    #             'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
    #             'basename': 'data_clinical_sample.txt',
    #             'class': 'File',
    #             'checksum': 'sha1$347a4f54e4490cac5ef1f67480957f0ea9337fbb',
    #             'size': 6487,
    #             'path':  os.path.join(output_dir,'data_clinical_sample.txt')
    #             }
    #         }
    #     self.assertDictEqual(output_json, expected_output)
    #
    #     output_file = expected_output['output_file']['path']
    #     table_reader = TableReader(output_file)
    #     comments = table_reader.comment_lines
    #     fieldnames = table_reader.get_fieldnames()
    #     records = [ rec for rec in table_reader.read() ]
    #
    #     expected_comments = [
    #     '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tSEX\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tCMO_TMB_SCORE\n',
    #     '#SAMPLE_ID\tIGO_ID\tPATIENT_ID\tCOLLAB_ID\tSAMPLE_TYPE\tSAMPLE_CLASS\tGENE_PANEL\tONCOTREE_CODE\tSPECIMEN_PRESERVATION_TYPE\tSEX\tTISSUE_SITE\tREQUEST_ID\tPROJECT_ID\tPIPELINE\tPIPELINE_VERSION\tCMO_TMB_SCORE\n',
    #     '#STRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tSTRING\tNUMBER\n',
    #     '#1\t1\t1\t0\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n'
    #     ]
    #     self.assertEqual(comments, expected_comments)
    #
    #     tmbs = {}
    #     for record in records:
    #         tmbs[record['SAMPLE_ID']] = record['CMO_TMB_SCORE']
    #
    #     expected_tmbs = {
    #     'Sample46': 'NA', 'Sample44': 'NA', 'Sample80': 'NA', 'Sample20': 'NA', 'Sample38': 'NA', 'Sample26': 'NA', 'Sample94': 'NA', 'Sample48': 'NA', 'Sample68': 'NA', 'Sample90': 'NA', 'Sample18': 'NA', 'Sample54': 'NA', 'Sample52': 'NA', 'Sample86': 'NA', 'Sample30': 'NA', 'Sample78': 'NA', 'Sample84': 'NA', 'Sample82': 'NA', 'Sample6': 'NA', 'Sample96': 'NA', 'Sample72': 'NA', 'Sample56': 'NA', 'Sample64': 'NA', 'Sample58': 'NA', 'Sample92': 'NA', 'Sample62': 'NA', 'Sample8': 'NA', 'Sample24': 'NA', 'Sample12': 'NA', 'Sample16': 'NA', 'Sample88': 'NA', 'Sample22': 'NA', 'Sample42': 'NA', 'Sample76': 'NA', 'Sample28': 'NA', 'Sample74': 'NA', 'Sample50': 'NA', 'Sample60': 'NA', 'Sample10': 'NA', 'Sample36': 'NA', 'Sample34': 'NA', 'Sample40': 'NA', 'Sample66': 'NA', 'Sample14': 'NA', 'Sample32': 'NA', 'Sample70': 'NA', 'Sample4': 'NA', 'Sample1': '0.000000013'
    #     }
    #     self.assertEqual(tmbs, expected_tmbs)


if __name__ == "__main__":
    unittest.main()
