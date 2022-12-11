#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the fillout_workflow cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile
)
sys.path.pop(0)

class TestFilloutWorkflow(PlutoTestCase):
    cwl_file = CWLFile('fillout_workflow.cwl')

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

        rows = [ self.maf_row1, self.maf_row2, self.maf_row3, self.maf_row4 ]
        lines = self.dicts2lines(rows, comment_list = self.comments)
        self.maf = self.write_table(tmpdir = self.tmpdir, filename = "consensus.maf", lines = lines)

    def test_run_fillout_workflow(self):
        """
        Test case for running the cohort fillout for multiple samples against a single maf
        """
        # self.preserve = True
        # self.maxDiff = None
        self.input = {
            "maf_file": {"class": "File", "path": self.maf},
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
            # NOTE: Must be in the same order as the sample bam files!
            "sample_ids":["Sample24", "Sample23"],
            # NOTE: Each must have a .bai file as well
            "bam_files": [
                {
                    "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
                },
                {
                    "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.maf'),
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$cad1317ab7c2940f11d91ce72cdb8a708c33108e',
                'size': 1871,
                'path':  os.path.join(output_dir,'output.maf')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'output.maf')

        lines = self.read_table(output_file)

        expected_lines = [
        ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Sequencer', 't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count', 'Caller', 't_total_count', 't_variant_frequency', 't_total_count_forward', 't_ref_count_forward', 't_alt_count_forward'],
        ['FAM46C', '.', 'mskcc.org', 'GRCh37', '1', '118166398', '118166398', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample24', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '0', '41', '.', '.', '.', '41', '1', '22', '0', '22'],
        ['IL7R', '.', 'mskcc.org', 'GRCh37', '5', '35876484', '35876484', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample24', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '0', '0', '.', '.', '.', '52', '0', '26', '0', '0'],
        ['KMT2C', '.', 'mskcc.org', 'GRCh37', '7', '151845367', '151845367', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample24', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '68', '4', '.', '.', '.', '72', '0.0555556', '34', '32', '2'],
        ['RTEL1', '.', 'mskcc.org', 'GRCh37', '20', '62321135', '62321135', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample24', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '129', '0', '.', '.', '.', '129', '0', '66', '66', '0'],
        ['FAM46C', '.', 'mskcc.org', 'GRCh37', '1', '118166398', '118166398', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample23', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '0', '40', '.', '.', '.', '40', '1', '20', '0', '20'],
        ['IL7R', '.', 'mskcc.org', 'GRCh37', '5', '35876484', '35876484', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample23', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '0', '0', '.', '.', '.', '120', '0', '58', '0', '0'],
        ['KMT2C', '.', 'mskcc.org', 'GRCh37', '7', '151845367', '151845367', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample23', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '91', '0', '.', '.', '.', '91', '0', '43', '43', '0'],
        ['RTEL1', '.', 'mskcc.org', 'GRCh37', '20', '62321135', '62321135', '+', 'Silent', 'SNP', 'G', 'A', '.', '.', '.', 'Sample23', 'Normal', '.', '.', '.', '.', '.', '.', '.', '.', 'UNPAIRED', '.', '.', '.', '.', '.', '.', '184', '0', '.', '.', '.', '184', '0', '89', '89', '0']
        ]

        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
