#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the index_sample_bam cwl
"""
import os
import sys
from datasets import (
    DATA_SETS,
)
from pluto import (
    PlutoTestCase,
    CWLFile,
    OFile
)


class TestIndexSampleBam(PlutoTestCase):
    cwl_file = CWLFile('index_sample_bam.cwl')

    def test_index1(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
        sample1_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample1.UnitTest01.bam')

        self.input = {
        "sample": { # "types.yml#FilloutIndexSample"
            "sample_id": "Sample1",
            "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
            "sample_type": "research",
            "prefilter": True,
            "maf_file": { "class": "File", "path": sample1_maf },
            "bam_file": { "class": "File", "path": sample1_bam }
            }
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'sample': {
                'bam_file': OFile(name = "Sample1.UnitTest01.bam", hash = "c77293b997de1f12e6ba155cde37eb06210af711", size = 4454012599, dir = output_dir,
                    secondaryFiles = [OFile(name = "Sample1.UnitTest01.bai", hash = "fc2d6cc4ebeac0fca2b6c5f19ab2a9930cb8f349", size = 4182280, dir = output_dir)]),
                'maf_file': OFile(name = "Sample1.FillOutUnitTest01.muts.maf", hash = "c49c5ea6e0f755d5229ff761e0cb10bc7df72a25", size = 101253, dir = output_dir),
                'normal_id': 'FROZENPOOLEDNORMAL_IMPACT505_V2',
                'prefilter': True,
                'sample_id': 'Sample1',
                'sample_type': 'research'
            },
        }
        self.assertCWLDictEqual(output_json, expected_output)
