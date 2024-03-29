#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the fusion_filter.cwl
"""
import os
import sys
from datasets import (
    DATA_SETS,
    KNOWN_FUSIONS_FILE,
)
from pluto import (
    CWLFile,
    PlutoTestCase,
    CWL_ARGS,
    OFile
)


class TestFusionFilter(PlutoTestCase):
    cwl_file = CWLFile('fusion_filter.cwl')

    def test_fusion_filter1(self):
        """
        """
        fusion_file = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt")

        self.input = {
            "fusions_file": {
                  "class": "File",
                  "path": fusion_file
                },
            "output_filename": "data_fusions.txt",
            "known_fusions_file": {
                "class": "File",
                "path": KNOWN_FUSIONS_FILE
            }
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "output_file": OFile(name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)
