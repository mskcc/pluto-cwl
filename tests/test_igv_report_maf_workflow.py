#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the igv-report_maf_workflow cwl

NOTE: for some reason, some of the variants used in the setUp method here are getting skipped along the workflow when getting converted to vcf
TODO: find a better way to convert maf to vcf, or just use maf with newer version of igv-reports
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
sys.path.pop(0)

class TestIGVReportMaf(PlutoTestCase):
    cwl_file = CWLFile('igv-report_maf_workflow.cwl')

    def setUp(self):
        super().setUp()
        self.comments = [
        ['# comment 1'],
        ['# comment 2']
        ]
        self.maf_row1 = OrderedDict([
        ('Hugo_Symbol', 'RTEL1'),
        ('Entrez_Gene_Id', '51750'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '20'),
        ('Start_Position', '62321135'),
        ('End_Position', '62321135'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row2 = OrderedDict([ # NOTE: This one is getting skipped in output
        ('Hugo_Symbol', 'FAM46C'),
        ('Entrez_Gene_Id', '54855'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '1'),
        ('Start_Position', '118166398'),
        ('End_Position', '118166398'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row3 = OrderedDict([ # NOTE: This one is getting skipped in output
        ('Hugo_Symbol', 'IL7R'),
        ('Entrez_Gene_Id', '3575'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '5'),
        ('Start_Position', '35876484'),
        ('End_Position', '35876484'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])
        self.maf_row4 = OrderedDict([
        ('Hugo_Symbol', 'KMT2C'),
        ('Entrez_Gene_Id', '58508'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '7'),
        ('Start_Position', '151845367'),
        ('End_Position', '151845367'),
        ('Variant_Classification', 'Silent'),
        ('Reference_Allele', 'G'),
        ('Tumor_Seq_Allele1', 'G'),
        ('Tumor_Seq_Allele2', 'A'),
        ('n_alt_count', '1'),
        ('Matched_Norm_Sample_Barcode', '.'),
        ('t_alt_count', '142'),
        ('t_ref_count', '511'),
        ('n_ref_count', '212'),
        ('Tumor_Sample_Barcode', '.')
        ])

        rows1 = [ self.maf_row1, self.maf_row2 ]
        lines1 = self.dicts2lines(rows1, comment_list = self.comments)
        self.maf1 = self.write_table(tmpdir = self.tmpdir, filename = "1.maf", lines = lines1)

        rows2 = [ self.maf_row3, self.maf_row4 ]
        lines2 = self.dicts2lines(rows2, comment_list = self.comments)
        self.maf2 = self.write_table(tmpdir = self.tmpdir, filename = "2.maf", lines = lines2)

        rows3 = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4 ]
        lines3 = self.dicts2lines(rows3, comment_list = self.comments)
        self.maf3 = self.write_table(tmpdir = self.tmpdir, filename = "3.maf", lines = lines3)

    def test_igv_report_demo1(self):
        """
        Test case with small bam files to simulate how files would be passed in from the main workflow
        NOTE: the maf file used should contain variants for which the bam has reads otherwise no track may appear
        """
        maf_file = os.path.join(self.DATA_SETS['demo']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        normal_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample2.bam")
        tumor_bam = os.path.join(self.DATA_SETS['demo']['BAM_DIR'], "Sample1.bam")
        report_filename = "Sample1.Sample2.igv.html"

        self.input = {
        # use the same single maf file for sites and for tracks to simulate how a tumor/normal pair would be handled
            "sites_maf": {"class": "File", "path": maf_file},
            "maf_files": [
                {"class": "File", "path": maf_file}
            ],
            "bam_files":[
                { "class": "File", "path": tumor_bam },
                { "class": "File", "path": normal_bam }
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
            "report_filename": report_filename
        }

        output_json, output_dir = self.run_cwl()

        # # do not include size and checksum since they are not consistent with .gz
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, report_filename),
                'basename': report_filename,
                'class': 'File',
                # 'checksum': 'sha1$bd56f50dc16e873707dc37b2a0ae96d36258ec1f',
                # 'size': 2596168,
                'path': os.path.join(output_dir, report_filename) }}

        self.assertEqual(output_json, expected_output)

    def test_igv_report(self):
        """
        Make an IGV report from a maf, a list of mafs, and bams
        """
        self.maxDiff = None
        self.input = {
            "sites_maf": {"class": "File", "path": self.maf3},
            "maf_files": [
                {"class": "File", "path": self.maf1},
                {"class": "File", "path": self.maf2}
            ],
            "bam_files":[
                { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") },
                { "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") }
            ],
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()

        # # do not include size and checksum since they are not consistent with .gz
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'igv.html'),
                'basename': 'igv.html',
                'class': 'File',
                # 'checksum': 'sha1$dcda65da7665683dcf97ed9c3989fd75c7a839c8',
                # 'size': 14885,
                'path': os.path.join(output_dir,'igv.html') }}

        self.assertEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
