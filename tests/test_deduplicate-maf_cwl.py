#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the deduplicate-maf cwl
"""
import os
import sys
from collections import OrderedDict
from pluto import (
    PlutoTestCase,
    CWLFile,
)

class TestDedupMaf(PlutoTestCase):
    cwl_file = CWLFile('deduplicate-maf.cwl')

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
        ('Variant_Classification', 'Silent'), # all the columns after this do not matter
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

    def test_dedup_maf1(self):
        """
        Test case for deduplicating and sorting a maf file
        """
        self.maxDiff = None
        rows1 = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row3, self.maf_row2, self.maf_row4 ]
        lines1 = self.dicts2lines(rows1, comment_list = self.comments)
        maf1 = self.write_table(tmpdir = self.tmpdir, filename = "1.maf", lines = lines1)

        self.input = {
            'input_file': {'class': 'File', 'path': maf1}
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'dedup.maf'),
                'basename': 'dedup.maf',
                'class': 'File',
                'checksum': 'sha1$7a9d60950b50f663c9fcd17e219db26412a3d27e',
                'size': 605,
                'path':  os.path.join(output_dir,'dedup.maf')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'dedup.maf')

        with open(output_file) as f:
            lines = [ line for line in f ]
        expected_lines = [
            '# comment 1\n',
            '# comment 2\n',
            'Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tn_alt_count\tMatched_Norm_Sample_Barcode\tt_alt_count\tt_ref_count\tn_ref_count\tTumor_Sample_Barcode\n',
            'FAM46C\t54855\tmskcc.org\tGRCh37\t1\t118166398\t118166398\tSilent\tG\tG\tA\t1\t.\t142\t511\t212\t.\n',
            'IL7R\t3575\tmskcc.org\tGRCh37\t5\t35876484\t35876484\tSilent\tG\tG\tA\t1\t.\t142\t511\t212\t.\n',
            'KMT2C\t58508\tmskcc.org\tGRCh37\t7\t151845367\t151845367\tSilent\tG\tG\tA\t1\t.\t142\t511\t212\t.\n',
            'RTEL1\t51750\tmskcc.org\tGRCh37\t20\t62321135\t62321135\tSilent\tG\tG\tA\t1\t.\t142\t511\t212\t.\n'
        ]
        self.assertEqual(lines, expected_lines)
