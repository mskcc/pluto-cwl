#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

$ CWL_ENGINE=Toil PRINT_TESTNAME=T python3 tests/test_fillout_singleton_processing.py
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

class TestSamplesFilloutIndexBatch(PlutoTestCase):
    cwl_file = CWLFile('fillout_singleton_processing.cwl')

    # these will get assigned during setup
    tmpdir = None
    output_json = None
    output_dir = None
    expected_output = None
    output_file = None
    filtered_output_path = None
    portal_output_path = None
    uncalled_output_path = None

    @classmethod
    def setUpClass(cls):
        """
        Need to make sure the pipeline only runs one time, then make the outputs available to the
        individual test cases

        NOTE: class gets imported once, but each test case gets executed in a separate instance of the class
        so we cannot assign values to instance "self", must be assigned to "cls"

        NOTE: dont think we can use super() here

        TODO: consider making this a standard part of PlutoTestCase somehow
        """
        # need to make an instance of the test case class in order to run it
        tc = cls()
        tc.setUp()
        output_json, output_dir = tc.setUpRun()

        # store the outputs on the class itself
        cls.tc = tc
        cls.tmpdir = tc.tmpdir
        cls.output_json = output_json
        cls.output_dir = output_dir
        cls.expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }
        cls.output_file = os.path.join(output_dir,'output.maf')
        cls.filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        cls.portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        cls.uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

    @classmethod
    def tearDownClass(cls):
        cls.tc.tearDown() # cls.rmtree(cls.tmpdir)

    def setUpRun(self):
        """
        Need to run setup methods then run the pipeline then return the results
        """
        sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')

        # NOTE: .bam files do not get used but are required for the type schema; "types.yml#FilloutIndexSample[]"
        sample3_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample3.UnitTest01.bam')
        sample4_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample4.UnitTest01.bam')
        sample5_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample5.UnitTest01.bam')

        samples = [
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample3_maf },
                "bam_file": { "class": "File", "path": sample3_bam }
            },
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample4_maf },
                "bam_file": { "class": "File", "path": sample4_bam }
            },
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample5_maf },
                "bam_file": { "class": "File", "path": sample5_bam }
            },
        ]
        self.input = {
            "samples": samples, # "types.yml#FilloutIndexSample[]"
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()
        return(output_json, output_dir)


    def test_CWLDict(self):
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(self.output_json, self.expected_output, related_keys = strip_related_keys)

    def test_muts_output(self):
        self.assertNumMutationsHash(self.output_file, 186, '39d85855a1b11f0f2b03d5c699811868')

    def test_muts_filtered_output(self):
        self.assertNumMutationsHash(self.filtered_output_path, 186, 'a4302c1291fd5399b8f5cf8487de93c0')

    def test_portal_muts(self):
        self.assertNumMutationsHash(self.portal_output_path, 77, 'bcbe84589dd21e02c0e8c4eb64e45925')

    def test_uncalled_muts(self):
        self.assertNumMutationsHash(self.uncalled_output_path, 109, '788fe2049c3cb23d0232f046fcb98ba3')





