#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile,
    DATA_SETS,
    OFile
)
sys.path.pop(0)

class TestFilloutClinicalFilter(PlutoTestCase):
    cwl_file = CWLFile('fillout_clinical_filter.cwl')

    def test_fillout_clinical_filter1(self):
        """
        """
        # NOTE: maf_file and bam_file do not actually get used for anything but are required by the type schema
        sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
        sample1_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample1.UnitTest01.bam')
        sample_ids = ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"]
        sample = {
                "sample_id": sample_ids[0],
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research", # NOTE: does not actually get used here, CWL uses clinical_sample_ids instead !!
                "maf_file": { "class": "File", "path": sample1_maf },
                "bam_file": { "class": "File", "path": sample1_bam }
            }

        self.input = {
            "sample": sample,
            "clinical_sample_ids": ["Sample3", "Sample4"],
            "fillout_vcf": {"class": "File", "path": os.path.join(self.DATA_SETS['Fillout01']['VCF_DIR'], "fillout.merged.sources.vcf")}
        }

        output_json, output_dir = self.run_cwl()
        expected_output = {
            "filtered_vcf": OFile(name = "Sample1.filtered.vcf", size = 11307, hash = "4787bef66c23be6fccdd75ebb03028960dbf5aca", dir = output_dir),
            "sample": {
                "bam_file": OFile(name = "Sample1.UnitTest01.bam", size = 4454012599, hash = "c77293b997de1f12e6ba155cde37eb06210af711", dir = output_dir,
                    secondaryFiles = [OFile(name = "Sample1.UnitTest01.bam.bai", size = 4182280, hash = "fc2d6cc4ebeac0fca2b6c5f19ab2a9930cb8f349", dir = output_dir)]),
                "filtered_maf": None,
                "filtered_vcf": OFile(name = "Sample1.filtered.vcf", hash = "c4866e67fc8f399212fc2bd9c774dd6118bf7b53", size = 11307, dir = output_dir),
                "maf_file": OFile(name = "Sample1.FillOutUnitTest01.muts.maf", size = 101253, hash = "c49c5ea6e0f755d5229ff761e0cb10bc7df72a25", dir = output_dir),
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_id": "Sample1",
                "sample_type": "research",
                "unfiltered_maf": None,
                "unfiltered_vcf": OFile(name = "Sample1.vcf", size = 20365, hash = "e1f905c248f253cefbee559781def5ad3979b0f2", dir = output_dir),
            },
            "unfiltered_vcf": OFile(name="Sample1.vcf", size = 20365, hash = "e1f905c248f253cefbee559781def5ad3979b0f2", dir = output_dir),
        }
        # fields inside the vcf are not static due to timestamps, etc..
        strip_related_keys = [
        ('basename', 'Sample1.filtered.vcf', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)
        self.assertNumMutationsHash(expected_output["filtered_vcf"]["path"], 45, "b60a9eb9127ee95f16ce2c566ae1d0df")

if __name__ == "__main__":
    unittest.main()
