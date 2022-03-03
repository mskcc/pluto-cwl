#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common location for some Python object fixtures that are used in various test cases for fillout workflows

A lot of these will be individual maf variant records to use for test cases
"""
import os
import sys
from collections import OrderedDict

comments = [
['# comment 1'],
['# comment 2']
]

# Sample24 - research sample
row1 = OrderedDict([
('Hugo_Symbol', 'RTEL1'),
('Entrez_Gene_Id', '51750'),
('Center', 'mskcc.org'),
('NCBI_Build', 'GRCh37'),
('Chromosome', '20'),
('Start_Position', '62321135'),
('End_Position', '62321135'),
('Matched_Norm_Sample_Barcode', 'Sample24-N'),
('Tumor_Sample_Barcode', 'Sample24'),
('Variant_Classification', 'Silent'), # all the columns after this do not matter # NOTE: is that true?
('Reference_Allele', 'G'),
('Tumor_Seq_Allele1', 'G'),
('Tumor_Seq_Allele2', 'A'),
('n_alt_count', '1'),
('t_alt_count', '142'),
('t_ref_count', '511'),
('n_ref_count', '212')
])
row2 = OrderedDict([
('Hugo_Symbol', 'FAM46C'),
('Entrez_Gene_Id', '54855'),
('Center', 'mskcc.org'),
('NCBI_Build', 'GRCh37'),
('Chromosome', '1'),
('Start_Position', '118166398'),
('End_Position', '118166398'),
('Matched_Norm_Sample_Barcode', 'Sample24-N'),
('Tumor_Sample_Barcode', 'Sample24'),
('Variant_Classification', 'Silent'),
('Reference_Allele', 'G'),
('Tumor_Seq_Allele1', 'G'),
('Tumor_Seq_Allele2', 'A'),
('n_alt_count', '1'),
('t_alt_count', '142'),
('t_ref_count', '511'),
('n_ref_count', '212')
])


# Sample23 - research sample
row3 = OrderedDict([
('Hugo_Symbol', 'IL7R'),
('Entrez_Gene_Id', '3575'),
('Center', 'mskcc.org'),
('NCBI_Build', 'GRCh37'),
('Chromosome', '5'),
('Start_Position', '35876484'),
('End_Position', '35876484'),
('Matched_Norm_Sample_Barcode', 'Sample23-N'),
('Tumor_Sample_Barcode', 'Sample23'),
('Variant_Classification', 'Silent'),
('Reference_Allele', 'G'),
('Tumor_Seq_Allele1', 'G'),
('Tumor_Seq_Allele2', 'A'),
('n_alt_count', '1'),
('t_alt_count', '142'),
('t_ref_count', '511'),
('n_ref_count', '212')
])
row4 = OrderedDict([
('Hugo_Symbol', 'KMT2C'),
('Entrez_Gene_Id', '58508'),
('Center', 'mskcc.org'),
('NCBI_Build', 'GRCh37'),
('Chromosome', '7'),
('Start_Position', '151845367'),
('End_Position', '151845367'),
('Matched_Norm_Sample_Barcode', 'Sample23-N'),
('Tumor_Sample_Barcode', 'Sample23'),
('Variant_Classification', 'Silent'),
('Reference_Allele', 'G'),
('Tumor_Seq_Allele1', 'G'),
('Tumor_Seq_Allele2', 'A'),
('n_alt_count', '1'),
('t_alt_count', '142'),
('t_ref_count', '511'),
('n_ref_count', '212')
])

class Rows(object):
    def __init__(self):
        self.comments = comments
        self.r1 = row1
        self.r2 = row2
        self.r3 = row3
        self.r4 = row4

rows = Rows()
