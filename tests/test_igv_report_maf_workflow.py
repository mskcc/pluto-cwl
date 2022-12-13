#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the igv-report_maf_workflow cwl
"""
import os
import sys
from collections import OrderedDict
from datasets import (
    DATA_SETS,
)
from pluto import (
    PlutoTestCase,
    CWLFile,
    TableReader
)


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
        self.maf_row2 = OrderedDict([
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
        self.maf_row3 = OrderedDict([
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

    def test_igv_report(self):
        """
        Make an IGV report from a maf, a list of mafs, and bams
        """
        self.skipTest("Fix unicode error")
        self.maxDiff = None
        self.input = {
            "sites_maf": {"class": "File", "path": self.maf3},
            "maf_files": [
                {"class": "File", "path": self.maf1},
                {"class": "File", "path": self.maf2}
            ],
            "bam_files":[
                { "class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") },
                { "class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") }
            ],
            "ref_fasta": {"class": "File", "path": DATA_SETS['Proj_08390_G']['REF_FASTA']},
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

        self.assertCWLDictEqual(output_json, expected_output)
