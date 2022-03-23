#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_workflow cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader, md5_obj
from pluto.settings import ENABLE_LARGE_TESTS, DATA_SETS
sys.path.pop(0)

from fixtures_fillout import rows

class TestSamplesFillout(PlutoTestCase):
    cwl_file = CWLFile('samples_fillout_workflow.cwl')

    def setUp(self):
        super().setUp()

        # Sample24
        lines1 = self.dicts2lines([ rows.r1, rows.r2 ], comment_list = rows.comments)
        self.maf1 = self.write_table(tmpdir = self.tmpdir, filename = "1.maf", lines = lines1)

        # Sample23
        lines2 = self.dicts2lines([ rows.r3, rows.r4 ], comment_list = rows.comments)
        self.maf2 = self.write_table(tmpdir = self.tmpdir, filename = "2.maf", lines = lines2)

    def test_small_1(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True
        # self.preserve = True
        # print(self.tmpdir)

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample24",
                    "normal_id": "Sample24-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": self.maf1 },
                    "bam_file": { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") }
                },
                {
                    "sample_id": "Sample23",
                    "normal_id": "Sample23-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": self.maf2 },
                    "bam_file": { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')
        output_portal_path = os.path.join(output_dir,'fillout.portal.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$70521e72bcf0b5c6ef8660178a8e54fa113d3efa',
                'size': 8202,
                'path':  output_path
                },
            'filtered_file': {
                'location': 'file://' + output_filtered_path,
                'basename': 'output.filtered.maf',
                'class': 'File',
                'checksum': 'sha1$70521e72bcf0b5c6ef8660178a8e54fa113d3efa',
                'size': 8202,
                'path':  output_filtered_path
                },
            'portal_file': {
                'basename': 'fillout.portal.maf',
                'checksum': 'sha1$399b7a685d755277c1842d3037d5102af32fc60f',
                'class': 'File',
                'location': 'file://' + output_portal_path,
                'path': output_portal_path,
                'size': 1380
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        reader = TableReader(output_path)
        comments = reader.comment_lines
        fieldnames = reader.get_fieldnames()
        records = [ rec for rec in reader.read() ]

        self.assertTrue(len(records) == 4)

        # subset the records to make sure the expected entries are present
        keep_cols = [
            'Hugo_Symbol',
            'Chromosome',
            'Start_Position',
            'End_Position',
            'Tumor_Sample_Barcode',
            't_depth',
            't_ref_count',
            't_alt_count',
            't_FL_VF',
            'is_fillout',
            'SRC'
            ]
        r = []
        for rec in records:
            r.append( dict((k, rec[k]) for k in keep_cols ) )

        expected_records = [
        {'Hugo_Symbol': 'KMT2C', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Tumor_Sample_Barcode': 'Sample24', 't_depth': '72', 't_ref_count': '68', 't_alt_count': '4', 't_FL_VF': '0.0555556', 'is_fillout': 'True', 'SRC': 'Sample23,'},
        {'Hugo_Symbol': 'RTEL1', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Tumor_Sample_Barcode': 'Sample24', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 't_FL_VF': '0', 'is_fillout': 'False', 'SRC': 'Sample24,'},
        {'Hugo_Symbol': 'KMT2C', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Tumor_Sample_Barcode': 'Sample23', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 't_FL_VF': '0', 'is_fillout': 'False', 'SRC': 'Sample23,'},
        {'Hugo_Symbol': 'RTEL1', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Tumor_Sample_Barcode': 'Sample23', 't_depth': '184', 't_ref_count': '184', 't_alt_count': '0', 't_FL_VF': '0', 'is_fillout': 'True', 'SRC': 'Sample24,'}
        ]
        self.assertEqual(r, expected_records)

        reader = TableReader(output_portal_path)
        comments = reader.comment_lines
        fieldnames = reader.get_fieldnames()
        records = [ rec for rec in reader.read() ]

        self.assertTrue(len(records) == 4)


    def test_small_1_clinical(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True
        # self.preserve = True
        # print(self.tmpdir)

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample24",
                    "normal_id": "Sample24-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": self.maf1 },
                    "bam_file": { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") }
                },
                {
                    "sample_id": "Sample23",
                    "normal_id": "Sample23-N",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": self.maf2 },
                    "bam_file": { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$70521e72bcf0b5c6ef8660178a8e54fa113d3efa',
                'size': 8202,
                'path':  output_path
                },
            'filtered_file': {
                'location': 'file://' + output_filtered_path,
                'basename': 'output.filtered.maf',
                'class': 'File',
                'checksum': 'sha1$cc604784c8b6bd0417f51fd32a912943b2375feb',
                'size': 8208,
                'path':  output_filtered_path
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        reader = TableReader(output_path)
        comments = reader.comment_lines
        fieldnames = reader.get_fieldnames()
        records = [ rec for rec in reader.read() ]

        self.assertTrue(len(records) == 4)

        # subset the records to make sure the expected entries are present
        keep_cols = [
            'Hugo_Symbol',
            'Chromosome',
            'Start_Position',
            'End_Position',
            'Tumor_Sample_Barcode',
            't_depth',
            't_ref_count',
            't_alt_count',
            't_FL_VF',
            'is_fillout',
            'SRC'
            ]
        r = []
        for rec in records:
            r.append( dict((k, rec[k]) for k in keep_cols ) )

        expected_records = [
        {'Hugo_Symbol': 'KMT2C', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Tumor_Sample_Barcode': 'Sample24', 't_depth': '72', 't_ref_count': '68', 't_alt_count': '4', 't_FL_VF': '0.0555556', 'is_fillout': 'True', 'SRC': 'Sample23,'},
        {'Hugo_Symbol': 'RTEL1', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Tumor_Sample_Barcode': 'Sample24', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 't_FL_VF': '0', 'is_fillout': 'False', 'SRC': 'Sample24,'},
        {'Hugo_Symbol': 'KMT2C', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Tumor_Sample_Barcode': 'Sample23', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 't_FL_VF': '0', 'is_fillout': 'False', 'SRC': 'Sample23,'},
        {'Hugo_Symbol': 'RTEL1', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Tumor_Sample_Barcode': 'Sample23', 't_depth': '184', 't_ref_count': '184', 't_alt_count': '0', 't_FL_VF': '0', 'is_fillout': 'True', 'SRC': 'Sample24,'}
        ]
        self.assertEqual(r, expected_records)


    def test_small_2(self):
        """
        Test case with small variant set that should have filters applied inside the pipeline

        run the fillout like this:

        sample_id       sample_type
        Sample3 research
        Sample5 research
        Sample2 research
        Sample1 research
        Sample4 research
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True
        # self.preserve = True
        # print(self.tmpdir)

        sample1_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample1.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample2_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample2.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample3_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample3.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample4.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample5.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')

        sample1_bam = os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample1.rg.md.abra.printreads.bam')
        sample2_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample2.rg.md.abra.printreads.bam')
        sample3_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample3.rg.md.abra.printreads.bam')
        sample4_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample4.rg.md.abra.printreads.bam')
        sample5_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample4.rg.md.abra.printreads.bam')

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample1_maf },
                    "bam_file": { "class": "File", "path": sample1_bam }
                },
                {
                    "sample_id": "Sample2",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample2_maf },
                    "bam_file": { "class": "File", "path": sample2_bam }
                },
                {
                    "sample_id": "Sample3",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample3_maf },
                    "bam_file": { "class": "File", "path": sample3_bam }
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample4_maf },
                    "bam_file": { "class": "File", "path": sample4_bam }
                },
                {
                    "sample_id": "Sample5",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample5_maf },
                    "bam_file": { "class": "File", "path": sample5_bam }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$2513c14c720e9e1ba02bb4a61fe0f31a80f60d12',
                # 'size': 114008492,
                'path':  output_path
                },
            'filtered_file': {
                'basename': 'output.filtered.maf',
                # 'checksum': 'sha1$90ec78a051db9186d7645cc76f125e6f20ccd077',
                'class': 'File',
                'location': 'file://' + output_filtered_path,
                'path': output_filtered_path
                # 'size': 115187724
                }
            }
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['filtered_file'].pop('checksum')
        output_json['filtered_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)
        # all_effects field is variable and changes bytes and checksum
        # need to check number of variant outputs instead

        comments, mutations = self.load_mutations(output_path, strip = True)
        self.assertEqual(len(mutations), 126975)
        hash = md5_obj(mutations)
        expected_hash = '80a0695ce2a2ee8a784b6092e36e0dd4'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(output_filtered_path, strip = True)
        self.assertEqual(len(mutations), 126975)
        hash = md5_obj(mutations)
        expected_hash = '80a0695ce2a2ee8a784b6092e36e0dd4'
        self.assertEqual(hash, expected_hash)


    def test_small_2_clinical(self):
        """
        Test case with small variant set that should have filters applied inside the pipeline

        run the fillout like this:

        sample_id       sample_type
        Sample3 clinical
        Sample5 research
        Sample2 research
        Sample1 research
        Sample4 clinical

        This test takes about 41min to complete
        """
        self.maxDiff = None
        self.runner_args['use_cache'] = False # do not use cache because it breaks for some reason
        self.runner_args['debug'] = True
        self.runner_args['js_console'] = True
        # self.preserve = True
        # print(self.tmpdir)

        sample1_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample1.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample2_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample2.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample3_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample3.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample4_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample4.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')
        sample5_maf = os.path.join(DATA_SETS['07618_AG']['MAF_DIR'], 'Sample5.FROZENPOOLEDNORMAL_IMPACT505_V2.muts.maf')

        sample1_bam = os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample1.rg.md.abra.printreads.bam')
        sample2_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample2.rg.md.abra.printreads.bam')
        sample3_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample3.rg.md.abra.printreads.bam')
        sample4_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample4.rg.md.abra.printreads.bam')
        sample5_bam =os.path.join(DATA_SETS['07618_AG']['BAM_DIR'], 'Sample4.rg.md.abra.printreads.bam')

        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample1_maf },
                    "bam_file": { "class": "File", "path": sample1_bam }
                },
                {
                    "sample_id": "Sample2",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample2_maf },
                    "bam_file": { "class": "File", "path": sample2_bam }
                },
                {
                    "sample_id": "Sample3",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": sample3_maf },
                    "bam_file": { "class": "File", "path": sample3_bam }
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": sample4_maf },
                    "bam_file": { "class": "File", "path": sample4_bam }
                },
                {
                    "sample_id": "Sample5",
                    "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": sample5_maf },
                    "bam_file": { "class": "File", "path": sample5_bam }
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$2513c14c720e9e1ba02bb4a61fe0f31a80f60d12',
                # 'size': 114008492,
                'path':  output_path
                },
            'filtered_file': {
                'basename': 'output.filtered.maf',
                # 'checksum': 'sha1$7800c1244d1b60b82e86f2fd3db87e1aff93afbc',
                'class': 'File',
                'location': 'file://' + filtered_output_path,
                'path': filtered_output_path
                # 'size': 110956123
                }
            }
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['filtered_file'].pop('checksum')
        output_json['filtered_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)
        # all_effects field is variable and changes bytes and checksum
        # need to check number of variant outputs instead

        comments, mutations = self.load_mutations(output_path, strip = True)
        self.assertEqual(len(mutations), 126975)
        hash = md5_obj(mutations)
        expected_hash = '80a0695ce2a2ee8a784b6092e36e0dd4'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(filtered_output_path, strip = True)
        self.assertEqual(len(mutations), 121699)
        hash = md5_obj(mutations)
        expected_hash = '4b55be411af84915ab07e01bb04d0619'
        self.assertEqual(hash, expected_hash)






    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_large_1(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        This test uses full samples
        """
        self.maxDiff = None
        maf1 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        maf24 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample24.Sample23.muts.maf")
        bam1 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample1.bam")
        bam24 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample24.bam")
        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "Sample1-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": maf1 },
                    "bam_file": { "class": "File", "path": bam1 },
                },
                {
                    "sample_id": "Sample24",
                    "normal_id": "Sample24-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": maf24 },
                    "bam_file": { "class": "File", "path": bam24 },
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_1']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')
        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$be8534bcaf326de029790a832ab5b44a17a03d22',
                # 'size': 40194610,
                'path':  output_path
                },
            'filtered_file': {
                'basename': 'output.filtered.maf',
                # 'checksum': 'sha1$f2dafd621df351af24820791f31bfc01d8dcd6ca',
                'class': 'File',
                'location': 'file://' + output_filtered_path,
                # 'size': 41164887
            }
            }
        # NOTE: for some reason, this file keeps coming out with different annotations for 'splice_acceptor_variant' or `splice_donor_variant`
        # this keeps changing the byte size and checksum so need to remove those here for now
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['filtered_file'].pop('checksum')
        output_json['filtered_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        # Need to remove some extra fields because they are inconsistent on the output maf file
        comments, mutations = self.load_mutations(output_path, strip = True)
        self.assertEqual(len(mutations), 38920)
        hash = md5_obj(mutations)
        expected_hash = '0cb1de1922f023a7157f5db273c9fe00'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(output_filtered_path, strip = True)
        self.assertEqual(len(mutations), 38920)
        hash = md5_obj(mutations)
        expected_hash = '0cb1de1922f023a7157f5db273c9fe00'
        self.assertEqual(hash, expected_hash)



    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_large_2(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        This test uses full samples
        """
        self.maxDiff = None
        maf1 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        maf4 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.muts.maf")
        bam1 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample1.bam")
        bam4 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample4.bam")
        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "Sample1-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": maf1 },
                    "bam_file": { "class": "File", "path": bam1 },
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "Sample4-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": maf4 },
                    "bam_file": { "class": "File", "path": bam4 },
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_1']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$2f60f58389ec65af87612c7532ad28b882fb84ba',
                # 'size': 26238820,
                'path':  output_path
                },
                'filtered_file': {
                    'basename': 'output.filtered.maf',
                    # 'checksum': 'sha1$f2dafd621df351af24820791f31bfc01d8dcd6ca',
                    'class': 'File',
                    'location': 'file://' + output_filtered_path,
                    # 'size': 41164887
                }
            }
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['filtered_file'].pop('checksum')
        output_json['filtered_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(output_path, strip = True)
        self.assertEqual(len(mutations), 26404)
        hash = md5_obj(mutations)
        expected_hash = 'eecc16c5910c9031831cbf4c44848796'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(output_filtered_path, strip = True)
        self.assertEqual(len(mutations), 26404)
        hash = md5_obj(mutations)
        expected_hash = 'eecc16c5910c9031831cbf4c44848796'
        self.assertEqual(hash, expected_hash)

    # @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_large_2_clinical(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        This test uses full samples

        takes about 9 minutes to complete
        """
        self.maxDiff = None
        maf1 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        maf4 = os.path.join(self.DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.muts.maf")
        bam1 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample1.bam")
        bam4 = os.path.join(self.DATA_SETS['Proj_1']['BAM_DIR'], "Sample4.bam")
        self.input = {
            "samples": [
                {
                    "sample_id": "Sample1",
                    "normal_id": "Sample1-N",
                    "sample_type": "research",
                    "maf_file": { "class": "File", "path": maf1 },
                    "bam_file": { "class": "File", "path": bam1 },
                },
                {
                    "sample_id": "Sample4",
                    "normal_id": "Sample4-N",
                    "sample_type": "clinical",
                    "maf_file": { "class": "File", "path": maf4 },
                    "bam_file": { "class": "File", "path": bam4 },
                },
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_1']['REF_FASTA']}
        }

        output_json, output_dir = self.run_cwl()
        output_path = os.path.join(output_dir,'output.maf')
        output_filtered_path = os.path.join(output_dir,'output.filtered.maf')

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$2f60f58389ec65af87612c7532ad28b882fb84ba',
                # 'size': 26238820,
                'path':  output_path
                },
                'filtered_file': {
                    'basename': 'output.filtered.maf',
                    # 'checksum': 'sha1$f2dafd621df351af24820791f31bfc01d8dcd6ca',
                    'class': 'File',
                    'location': 'file://' + output_filtered_path,
                    # 'size': 41164887
                }
            }
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['filtered_file'].pop('checksum')
        output_json['filtered_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(output_path, strip = True)
        self.assertEqual(len(mutations), 26404)
        hash = md5_obj(mutations)
        expected_hash = 'eecc16c5910c9031831cbf4c44848796'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(output_filtered_path, strip = True)
        self.assertEqual(len(mutations), 26403) # 25946
        hash = md5_obj(mutations)
        expected_hash = 'ad7dee2dc0c96d04caaf13bdfcc2d45b'
        self.assertEqual(hash, expected_hash)



if __name__ == "__main__":
    unittest.main()
