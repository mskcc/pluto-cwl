#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the consensus_maf cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestConcatMafs(PlutoTestCase):
    cwl_file = CWLFile('consensus_maf.cwl')

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
        ('Variant_Classification', 'Silent'), # all the columns after this do not matter in the output and will get blank'd
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
        self.maf_row5 = OrderedDict([
        ('Hugo_Symbol', 'MET'),
        ('Entrez_Gene_Id', '4233'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '7'),
        ('Start_Position', '116418998'),
        ('End_Position', '116418998'),
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
        self.maf_row6 = OrderedDict([
        ('Hugo_Symbol', 'MAP2K4'),
        ('Entrez_Gene_Id', '6416'),
        ('Center', 'mskcc.org'),
        ('NCBI_Build', 'GRCh37'),
        ('Chromosome', '17'),
        ('Start_Position', '11998935'),
        ('End_Position', '11998935'),
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

    def test_consensus_maf1(self):
        """
        Test case for making a consensus maf from two maf files
        """
        self.maxDiff = None

        rows1 = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4 ]
        lines1 = self.dicts2lines(rows1, comment_list = self.comments)
        maf1 = self.write_table(tmpdir = self.tmpdir, filename = "1.maf", lines = lines1)

        rows2 = [ self.maf_row3, self.maf_row4, self.maf_row5, self.maf_row6 ]
        lines2 = self.dicts2lines(rows2, comment_list = self.comments)
        maf2 = self.write_table(tmpdir = self.tmpdir, filename = "2.maf", lines = lines2)

        self.input = {
            'maf_files': [
                {'class': 'File', 'path': maf1},
                {'class': 'File', 'path': maf2},
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'dedup.maf'),
                'basename': 'dedup.maf',
                'class': 'File',
                'checksum': 'sha1$2050566c000881a879498b55e3eabd42c9a29f60',
                'size': 652,
                'path':  os.path.join(output_dir,'dedup.maf')
                }
            }

        self.assertEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']

        with open(output_file) as f:
            lines = [ line for line in f ]

        expected_lines = [
            '# comment 1\n',
            '# comment 2\n',
            'Hugo_Symbol\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\n',
            'FAM46C\tmskcc.org\tGRCh37\t1\t118166398\t118166398\tSilent\tG\tG\tA\t\t\t\t\t\t\n',
            'IL7R\tmskcc.org\tGRCh37\t5\t35876484\t35876484\tSilent\tG\tG\tA\t\t\t\t\t\t\n',
            'MET\tmskcc.org\tGRCh37\t7\t116418998\t116418998\tSilent\tG\tG\tA\t\t\t\t\t\t\n',
            'KMT2C\tmskcc.org\tGRCh37\t7\t151845367\t151845367\tSilent\tG\tG\tA\t\t\t\t\t\t\n',
            'MAP2K4\tmskcc.org\tGRCh37\t17\t11998935\t11998935\tSilent\tG\tG\tA\t\t\t\t\t\t\n',
            'RTEL1\tmskcc.org\tGRCh37\t20\t62321135\t62321135\tSilent\tG\tG\tA\t\t\t\t\t\t\n'
        ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
