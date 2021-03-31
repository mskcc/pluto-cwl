#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the MSI analysis cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.settings import ENABLE_LARGE_TESTS
sys.path.pop(0)

class TestMSI(PlutoTestCase):
    cwl_file = CWLFile('msi.cwl')

    def test_msi_demo1(self):
        """
        Test case for running MSI on small demo dataset
        """
        self.maxDiff = None
        tumor_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")
        normal_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        microsatellites_file = self.DATA_SETS['demo']['microsatellites_file']

        self.input = {
              "d":    {"class": "File", "path": microsatellites_file },
              "n":    {"class": "File", "path": normal_bam },
              "t":    {"class": "File", "path": tumor_bam },
              "o": "Sample1.Sample2"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, "Sample1.Sample2"),
                'basename': "Sample1.Sample2",
                'class': 'File',
                'checksum': 'sha1$d307badbf88a82ce0c368d5d41d816d0e053da4c',
                'size': 59,
                'path': os.path.join(output_dir, "Sample1.Sample2")
                },
            'dis_file': {
                'basename': 'Sample1.Sample2_dis',
                'checksum': 'sha1$0ba8807a47e0cfe14533e4bcaf80219687b59857',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir, 'Sample1.Sample2_dis'),
                'path': os.path.join(output_dir, 'Sample1.Sample2_dis'),
                'size': 14687751
                },
           'somatic_file': {
               'basename': 'Sample1.Sample2_somatic',
                'checksum': 'sha1$e9391ea6b9063aafe8f9718932abe2429fc253da',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir, 'Sample1.Sample2_somatic'),
                'path': os.path.join(output_dir, 'Sample1.Sample2_somatic'),
                'size': 397
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        lines = self.read_table(output_file)

        expected_lines = [
            ['Total_Number_of_Sites', 'Number_of_Somatic_Sites', '%'],
            ['25', '7', '28.00']
            ]
        self.assertEqual(lines, expected_lines)

    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_msi1(self):
        """
        Test case for the MSI analysis workflow
        """
        self.maxDiff = None
        tumor_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
        normal_bam = os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        microsatellites_file = self.DATA_SETS['Proj_08390_G']['microsatellites_file']

        self.input = {
              "d":    {"class": "File", "path": microsatellites_file },
              "n":    {"class": "File", "path": normal_bam },
              "t":    {"class": "File", "path": tumor_bam },
              "o": "Sample24.Sample23"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, "Sample24.Sample23"),
                'basename': "Sample24.Sample23",
                'class': 'File',
                'checksum': 'sha1$8d5856202b859fc0da427bf9069eb870f2adcd55',
                'size': 62,
                'path': os.path.join(output_dir, "Sample24.Sample23")
                },
            'dis_file': {
                'basename': 'Sample24.Sample23_dis',
                'checksum': 'sha1$4129165a2034041595052a4e94c28c58c8d44a00',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_dis'),
                'path': os.path.join(output_dir, 'Sample24.Sample23_dis'),
                'size': 772099729
                },
           'somatic_file': {
               'basename': 'Sample24.Sample23_somatic',
                'checksum': 'sha1$ad5556df9f02aee0e43d74a4b7855d5d8cfc4def',
                'class': 'File',
                'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_somatic'),
                'path': os.path.join(output_dir, 'Sample24.Sample23_somatic'),
                'size': 8432
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        lines = self.read_table(output_file)

        expected_lines = [
            ['Total_Number_of_Sites', 'Number_of_Somatic_Sites', '%'],
            ['628', '138', '21.97']
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
