#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the portal-workflow.cwl
"""
import os
import sys
import unittest



from pluto import (
    CWLFile, 
    PlutoTestCase,
    DATA_SETS, 
    KNOWN_FUSIONS_FILE, 
    ENABLE_LARGE_TESTS,
    OFile, 
    ODir
)



class TestFusionToSV(PlutoTestCase):
    cwl_file = CWLFile('fusion_to_sv.cwl')

    def test_fusion_to_sv(self):
        """
        Test fusion to sv conversion
        """
        fusion_file = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt")

        self.input = {
            "fusion_file": {
                  "class": "File",
                  "path": fusion_file
                },
            "output_filename": "data_SV.txt"

        }        

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(
                name='data_SV.txt', size=1103, hash='02fda70b7838931321544f6797de4782adaf1a46', dir=output_dir)

            
        }

        self.maxDiff = None
        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

       



