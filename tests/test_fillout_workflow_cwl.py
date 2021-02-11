#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the fillout_workflow cwl
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

class TestGetBaseCounts(PlutoTestCase):
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
        self.maxDiff = None
        self.input = {
            "maf_file": {"class": "File", "path": self.maf},
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
            "bams": [
                {
                "bam_file": {
                    "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
                    },
                "sample_id": "Sample24"
                },
                {
                "bam_file": {
                    "class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
                    },
                "sample_id": "Sample23"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'fillout.maf'),
                'basename': 'fillout.maf',
                'class': 'File',
                'checksum': 'sha1$60d27c22093c202607cacc91834d630fc8239352',
                'size': 1201,
                'path':  os.path.join(output_dir,'fillout.maf')
                }
            }
        self.assertEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']

        with open(output_file) as f:
            lines = [ line for line in f ]

        print(lines)

        expected_lines = [
            'Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\tdbSNP_Val_Status\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tt_ref_count\tt_alt_count\tn_ref_count\tn_alt_count\tCaller\tt_total_count\tt_variant_frequency\tt_total_count_forward\tt_ref_count_forward\tt_alt_count_forward\n',
            'FAM46C\t\tmskcc.org\tGRCh37\t1\t118166398\t118166398\t+\tSilent\tSNP\tG\tA\t\t\t\tSample24\tNormal\t\t\t\t\t\t\t\t\tUNPAIRED\t\t\t\t\t\t\t0\t41\t\t\t\t41\t1\t22\t0\t22\n',
            'IL7R\t\tmskcc.org\tGRCh37\t5\t35876484\t35876484\t+\tSilent\tSNP\tG\tA\t\t\t\tSample24\tNormal\t\t\t\t\t\t\t\t\tUNPAIRED\t\t\t\t\t\t\t0\t0\t\t\t\t52\t0\t26\t0\t0\n',
            'KMT2C\t\tmskcc.org\tGRCh37\t7\t151845367\t151845367\t+\tSilent\tSNP\tG\tA\t\t\t\tSample24\tNormal\t\t\t\t\t\t\t\t\tUNPAIRED\t\t\t\t\t\t\t68\t4\t\t\t\t72\t0.0555556\t34\t32\t2\n',
            'RTEL1\t\tmskcc.org\tGRCh37\t20\t62321135\t62321135\t+\tSilent\tSNP\tG\tA\t\t\t\tSample24\tNormal\t\t\t\t\t\t\t\t\tUNPAIRED\t\t\t\t\t\t\t129\t0\t\t\t\t129\t0\t66\t66\t0\n'
        ]

        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
