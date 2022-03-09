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

    def test_run_fillout_workflow_small_1(self):
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

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                'checksum': 'sha1$7932ae9938a5686f6328f143a6c82308877cb822',
                'size': 8008,
                'path':  output_path
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        reader = TableReader(output_path)
        comments = reader.comment_lines
        fieldnames = reader.get_fieldnames()
        records = [ rec for rec in reader.read() ]

        self.assertTrue(len(records) == 4)

        expected_records = [
        {'Hugo_Symbol': 'KMT2C', 'Entrez_Gene_Id': '58508', 'Center': '.', 'NCBI_Build': 'GRCh37', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Strand': '+', 'Variant_Classification': 'Missense_Mutation', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'G', 'Tumor_Seq_Allele2': 'A', 'dbSNP_RS': '', 'dbSNP_Val_Status': '', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Sample24-N', 'Match_Norm_Seq_Allele1': 'G', 'Match_Norm_Seq_Allele2': 'G', 'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '', 'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '', 'Verification_Status': '', 'Validation_Status': '', 'Mutation_Status': '', 'Sequencing_Phase': '', 'Sequence_Source': '', 'Validation_Method': '', 'Score': '', 'BAM_File': '', 'Sequencer': '', 'Tumor_Sample_UUID': '', 'Matched_Norm_Sample_UUID': '', 'HGVSc': 'c.13645C>T', 'HGVSp': 'p.Arg4549Cys', 'HGVSp_Short': 'p.R4549C', 'Transcript_ID': 'ENST00000262189', 'Exon_Number': '52/59', 't_depth': '72', 't_ref_count': '68', 't_alt_count': '4', 'n_depth': '', 'n_ref_count': '', 'n_alt_count': '', 'all_effects': 'KMT2C,missense_variant,p.Arg4606Cys,ENST00000355193,;KMT2C,missense_variant,p.Arg4549Cys,ENST00000262189,NM_170606.2;KMT2C,missense_variant,p.Arg2110Cys,ENST00000360104,;KMT2C,missense_variant,p.Arg1166Cys,ENST00000424877,;KMT2C,downstream_gene_variant,,ENST00000418061,;KMT2C,downstream_gene_variant,,ENST00000485241,;KMT2C,3_prime_UTR_variant,,ENST00000558084,;KMT2C,non_coding_transcript_exon_variant,,ENST00000473186,;', 'Allele': 'A', 'Gene': 'ENSG00000055609', 'Feature': 'ENST00000262189', 'Feature_type': 'Transcript', 'Consequence': 'missense_variant', 'cDNA_position': '13864/16862', 'CDS_position': '13645/14736', 'Protein_position': '4549/4911', 'Amino_acids': 'R/C', 'Codons': 'Cgc/Tgc', 'Existing_variation': 'COSM245709,COSM245710', 'ALLELE_NUM': '1', 'DISTANCE': '', 'STRAND_VEP': '-1', 'SYMBOL': 'KMT2C', 'SYMBOL_SOURCE': 'HGNC', 'HGNC_ID': '13726', 'BIOTYPE': 'protein_coding', 'CANONICAL': 'YES', 'CCDS': 'CCDS5931.1', 'ENSP': 'ENSP00000262189', 'SWISSPROT': 'Q8NEZ4', 'TREMBL': 'Q6N019,Q75MN6,H0YMU7', 'UNIPARC': 'UPI0000141B9F', 'RefSeq': 'NM_170606.2', 'SIFT': '', 'PolyPhen': 'probably_damaging(0.999)', 'EXON': '52/59', 'INTRON': '', 'DOMAINS': 'PROSITE_profiles:PS51542,hmmpanther:PTHR22884,hmmpanther:PTHR22884:SF305', 'AF': '', 'AFR_AF': '', 'AMR_AF': '', 'ASN_AF': '', 'EAS_AF': '', 'EUR_AF': '', 'SAS_AF': '', 'AA_AF': '', 'EA_AF': '', 'CLIN_SIG': '', 'SOMATIC': '1,1', 'PUBMED': '', 'MOTIF_NAME': '', 'MOTIF_POS': '', 'HIGH_INF_POS': '', 'MOTIF_SCORE_CHANGE': '', 'IMPACT': 'MODERATE', 'PICK': '1', 'VARIANT_CLASS': 'SNV', 'TSL': '', 'HGVS_OFFSET': '', 'PHENO': '1,1', 'MINIMISED': '', 'ExAC_AF': '', 'ExAC_AF_AFR': '', 'ExAC_AF_AMR': '', 'ExAC_AF_EAS': '', 'ExAC_AF_FIN': '', 'ExAC_AF_NFE': '', 'ExAC_AF_OTH': '', 'ExAC_AF_SAS': '', 'GENE_PHENO': '', 'FILTER': '.', 'flanking_bps': 'CGA', 'vcf_id': '.', 'vcf_qual': '.', 'ExAC_AF_Adj': '', 'ExAC_AC_AN_Adj': '', 'ExAC_AC_AN': '', 'ExAC_AC_AN_AFR': '', 'ExAC_AC_AN_AMR': '', 'ExAC_AC_AN_EAS': '', 'ExAC_AC_AN_FIN': '', 'ExAC_AC_AN_NFE': '', 'ExAC_AC_AN_OTH': '', 'ExAC_AC_AN_SAS': '', 'ExAC_FILTER': '', 'gnomAD_AF': '', 'gnomAD_AFR_AF': '', 'gnomAD_AMR_AF': '', 'gnomAD_ASJ_AF': '', 'gnomAD_EAS_AF': '', 'gnomAD_FIN_AF': '', 'gnomAD_NFE_AF': '', 'gnomAD_OTH_AF': '', 'gnomAD_SAS_AF': '', 'vcf_pos': '151845367', 'AC': '1', 'AN': '2', 'SRC': 'Sample23,', 't_GT': './.', 'n_GT': '', 't_FL_AD': '4', 'n_FL_AD': '', 't_FL_ADN': '2', 'n_FL_ADN': '', 't_FL_ADP': '2', 'n_FL_ADP': '', 't_FL_DP': '72', 'n_FL_DP': '', 't_FL_DPN': '38', 'n_FL_DPN': '', 't_FL_DPP': '34', 'n_FL_DPP': '', 't_FL_RD': '68', 'n_FL_RD': '', 't_FL_RDN': '36', 'n_FL_RDN': '', 't_FL_RDP': '32', 'n_FL_RDP': '', 't_FL_VF': '0.0555556', 'n_FL_VF': '', 't_AD': '', 'n_AD': '', 't_DP': '.', 'n_DP': '', 't_depth_sample': '', 't_ref_count_sample': '', 't_alt_count_sample': '', 'is_fillout': 'True'},

        {'Hugo_Symbol': 'RTEL1', 'Entrez_Gene_Id': '51750', 'Center': '.', 'NCBI_Build': 'GRCh37', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'G', 'Tumor_Seq_Allele2': 'A', 'dbSNP_RS': 'rs746824222', 'dbSNP_Val_Status': '', 'Tumor_Sample_Barcode': 'Sample24', 'Matched_Norm_Sample_Barcode': 'Sample24-N', 'Match_Norm_Seq_Allele1': 'G', 'Match_Norm_Seq_Allele2': 'G', 'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '', 'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '', 'Verification_Status': '', 'Validation_Status': '', 'Mutation_Status': '', 'Sequencing_Phase': '', 'Sequence_Source': '', 'Validation_Method': '', 'Score': '', 'BAM_File': '', 'Sequencer': '', 'Tumor_Sample_UUID': '', 'Matched_Norm_Sample_UUID': '', 'HGVSc': 'c.2130G>A', 'HGVSp': 'p.=', 'HGVSp_Short': 'p.Q710=', 'Transcript_ID': 'ENST00000508582', 'Exon_Number': '24/35', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 'n_depth': '', 'n_ref_count': '', 'n_alt_count': '', 'all_effects': 'RTEL1,synonymous_variant,p.=,ENST00000318100,;RTEL1,synonymous_variant,p.=,ENST00000370018,NM_032957.4,NM_016434.3;RTEL1,synonymous_variant,p.=,ENST00000360203,NM_001283009.1;RTEL1,synonymous_variant,p.=,ENST00000508582,;RTEL1,synonymous_variant,p.=,ENST00000425905,;RTEL1,upstream_gene_variant,,ENST00000370003,;RTEL1-TNFRSF6B,synonymous_variant,p.=,ENST00000482936,;RTEL1-TNFRSF6B,synonymous_variant,p.=,ENST00000492259,;RTEL1-TNFRSF6B,non_coding_transcript_exon_variant,,ENST00000480273,;RTEL1-TNFRSF6B,non_coding_transcript_exon_variant,,ENST00000496281,;RTEL1,upstream_gene_variant,,ENST00000496816,;', 'Allele': 'A', 'Gene': 'ENSG00000258366', 'Feature': 'ENST00000508582', 'Feature_type': 'Transcript', 'Consequence': 'synonymous_variant', 'cDNA_position': '2476/4273', 'CDS_position': '2130/3732', 'Protein_position': '710/1243', 'Amino_acids': 'Q', 'Codons': 'caG/caA', 'Existing_variation': 'rs746824222', 'ALLELE_NUM': '1', 'DISTANCE': '', 'STRAND_VEP': '1', 'SYMBOL': 'RTEL1', 'SYMBOL_SOURCE': 'HGNC', 'HGNC_ID': '15888', 'BIOTYPE': 'protein_coding', 'CANONICAL': 'YES', 'CCDS': 'CCDS13530.3', 'ENSP': 'ENSP00000424307', 'SWISSPROT': 'Q9NZ71', 'TREMBL': '', 'UNIPARC': 'UPI00019B2219', 'RefSeq': '', 'SIFT': '', 'PolyPhen': '', 'EXON': '24/35', 'INTRON': '', 'DOMAINS': 'Superfamily_domains:SSF52540,SMART_domains:SM00491,Pfam_domain:PF13307,TIGRFAM_domain:TIGR00604,hmmpanther:PTHR11472:SF4,hmmpanther:PTHR11472', 'AF': '', 'AFR_AF': '', 'AMR_AF': '', 'ASN_AF': '', 'EAS_AF': '', 'EUR_AF': '', 'SAS_AF': '', 'AA_AF': '', 'EA_AF': '', 'CLIN_SIG': '', 'SOMATIC': '', 'PUBMED': '', 'MOTIF_NAME': '', 'MOTIF_POS': '', 'HIGH_INF_POS': '', 'MOTIF_SCORE_CHANGE': '', 'IMPACT': 'LOW', 'PICK': '1', 'VARIANT_CLASS': 'SNV', 'TSL': '', 'HGVS_OFFSET': '', 'PHENO': '', 'MINIMISED': '', 'ExAC_AF': '', 'ExAC_AF_AFR': '', 'ExAC_AF_AMR': '', 'ExAC_AF_EAS': '', 'ExAC_AF_FIN': '', 'ExAC_AF_NFE': '', 'ExAC_AF_OTH': '', 'ExAC_AF_SAS': '', 'GENE_PHENO': '1', 'FILTER': '.', 'flanking_bps': 'AGG', 'vcf_id': '.', 'vcf_qual': '.', 'ExAC_AF_Adj': '', 'ExAC_AC_AN_Adj': '', 'ExAC_AC_AN': '', 'ExAC_AC_AN_AFR': '', 'ExAC_AC_AN_AMR': '', 'ExAC_AC_AN_EAS': '', 'ExAC_AC_AN_FIN': '', 'ExAC_AC_AN_NFE': '', 'ExAC_AC_AN_OTH': '', 'ExAC_AC_AN_SAS': '', 'ExAC_FILTER': '', 'gnomAD_AF': '', 'gnomAD_AFR_AF': '', 'gnomAD_AMR_AF': '', 'gnomAD_ASJ_AF': '', 'gnomAD_EAS_AF': '', 'gnomAD_FIN_AF': '', 'gnomAD_NFE_AF': '', 'gnomAD_OTH_AF': '', 'gnomAD_SAS_AF': '', 'vcf_pos': '62321135', 'AC': '1', 'AN': '2', 'SRC': 'Sample24,', 't_GT': '0/1', 'n_GT': '', 't_FL_AD': '0', 'n_FL_AD': '', 't_FL_ADN': '0', 'n_FL_ADN': '', 't_FL_ADP': '0', 'n_FL_ADP': '', 't_FL_DP': '129', 'n_FL_DP': '', 't_FL_DPN': '63', 'n_FL_DPN': '', 't_FL_DPP': '66', 'n_FL_DPP': '', 't_FL_RD': '129', 'n_FL_RD': '', 't_FL_RDN': '63', 'n_FL_RDN': '', 't_FL_RDP': '66', 'n_FL_RDP': '', 't_FL_VF': '0', 'n_FL_VF': '', 't_AD': '511,142', 'n_AD': '', 't_DP': '653', 'n_DP': '', 't_depth_sample': '653', 't_ref_count_sample': '511', 't_alt_count_sample': '142', 'is_fillout': 'False'},

        {'Hugo_Symbol': 'KMT2C', 'Entrez_Gene_Id': '58508', 'Center': '.', 'NCBI_Build': 'GRCh37', 'Chromosome': '7', 'Start_Position': '151845367', 'End_Position': '151845367', 'Strand': '+', 'Variant_Classification': 'Missense_Mutation', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'G', 'Tumor_Seq_Allele2': 'A', 'dbSNP_RS': '', 'dbSNP_Val_Status': '', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Sample23-N', 'Match_Norm_Seq_Allele1': 'G', 'Match_Norm_Seq_Allele2': 'G', 'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '', 'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '', 'Verification_Status': '', 'Validation_Status': '', 'Mutation_Status': '', 'Sequencing_Phase': '', 'Sequence_Source': '', 'Validation_Method': '', 'Score': '', 'BAM_File': '', 'Sequencer': '', 'Tumor_Sample_UUID': '', 'Matched_Norm_Sample_UUID': '', 'HGVSc': 'c.13645C>T', 'HGVSp': 'p.Arg4549Cys', 'HGVSp_Short': 'p.R4549C', 'Transcript_ID': 'ENST00000262189', 'Exon_Number': '52/59', 't_depth': '653', 't_ref_count': '511', 't_alt_count': '142', 'n_depth': '', 'n_ref_count': '', 'n_alt_count': '', 'all_effects': 'KMT2C,missense_variant,p.Arg4606Cys,ENST00000355193,;KMT2C,missense_variant,p.Arg4549Cys,ENST00000262189,NM_170606.2;KMT2C,missense_variant,p.Arg2110Cys,ENST00000360104,;KMT2C,missense_variant,p.Arg1166Cys,ENST00000424877,;KMT2C,downstream_gene_variant,,ENST00000418061,;KMT2C,downstream_gene_variant,,ENST00000485241,;KMT2C,3_prime_UTR_variant,,ENST00000558084,;KMT2C,non_coding_transcript_exon_variant,,ENST00000473186,;', 'Allele': 'A', 'Gene': 'ENSG00000055609', 'Feature': 'ENST00000262189', 'Feature_type': 'Transcript', 'Consequence': 'missense_variant', 'cDNA_position': '13864/16862', 'CDS_position': '13645/14736', 'Protein_position': '4549/4911', 'Amino_acids': 'R/C', 'Codons': 'Cgc/Tgc', 'Existing_variation': 'COSM245709,COSM245710', 'ALLELE_NUM': '1', 'DISTANCE': '', 'STRAND_VEP': '-1', 'SYMBOL': 'KMT2C', 'SYMBOL_SOURCE': 'HGNC', 'HGNC_ID': '13726', 'BIOTYPE': 'protein_coding', 'CANONICAL': 'YES', 'CCDS': 'CCDS5931.1', 'ENSP': 'ENSP00000262189', 'SWISSPROT': 'Q8NEZ4', 'TREMBL': 'Q6N019,Q75MN6,H0YMU7', 'UNIPARC': 'UPI0000141B9F', 'RefSeq': 'NM_170606.2', 'SIFT': '', 'PolyPhen': 'probably_damaging(0.999)', 'EXON': '52/59', 'INTRON': '', 'DOMAINS': 'PROSITE_profiles:PS51542,hmmpanther:PTHR22884,hmmpanther:PTHR22884:SF305', 'AF': '', 'AFR_AF': '', 'AMR_AF': '', 'ASN_AF': '', 'EAS_AF': '', 'EUR_AF': '', 'SAS_AF': '', 'AA_AF': '', 'EA_AF': '', 'CLIN_SIG': '', 'SOMATIC': '1,1', 'PUBMED': '', 'MOTIF_NAME': '', 'MOTIF_POS': '', 'HIGH_INF_POS': '', 'MOTIF_SCORE_CHANGE': '', 'IMPACT': 'MODERATE', 'PICK': '1', 'VARIANT_CLASS': 'SNV', 'TSL': '', 'HGVS_OFFSET': '', 'PHENO': '1,1', 'MINIMISED': '', 'ExAC_AF': '', 'ExAC_AF_AFR': '', 'ExAC_AF_AMR': '', 'ExAC_AF_EAS': '', 'ExAC_AF_FIN': '', 'ExAC_AF_NFE': '', 'ExAC_AF_OTH': '', 'ExAC_AF_SAS': '', 'GENE_PHENO': '', 'FILTER': '.', 'flanking_bps': 'CGA', 'vcf_id': '.', 'vcf_qual': '.', 'ExAC_AF_Adj': '', 'ExAC_AC_AN_Adj': '', 'ExAC_AC_AN': '', 'ExAC_AC_AN_AFR': '', 'ExAC_AC_AN_AMR': '', 'ExAC_AC_AN_EAS': '', 'ExAC_AC_AN_FIN': '', 'ExAC_AC_AN_NFE': '', 'ExAC_AC_AN_OTH': '', 'ExAC_AC_AN_SAS': '', 'ExAC_FILTER': '', 'gnomAD_AF': '', 'gnomAD_AFR_AF': '', 'gnomAD_AMR_AF': '', 'gnomAD_ASJ_AF': '', 'gnomAD_EAS_AF': '', 'gnomAD_FIN_AF': '', 'gnomAD_NFE_AF': '', 'gnomAD_OTH_AF': '', 'gnomAD_SAS_AF': '', 'vcf_pos': '151845367', 'AC': '1', 'AN': '2', 'SRC': 'Sample23,', 't_GT': '0/1', 'n_GT': '', 't_FL_AD': '0', 'n_FL_AD': '', 't_FL_ADN': '0', 'n_FL_ADN': '', 't_FL_ADP': '0', 'n_FL_ADP': '', 't_FL_DP': '91', 'n_FL_DP': '', 't_FL_DPN': '48', 'n_FL_DPN': '', 't_FL_DPP': '43', 'n_FL_DPP': '', 't_FL_RD': '91', 'n_FL_RD': '', 't_FL_RDN': '48', 'n_FL_RDN': '', 't_FL_RDP': '43', 'n_FL_RDP': '', 't_FL_VF': '0', 'n_FL_VF': '', 't_AD': '511,142', 'n_AD': '', 't_DP': '653', 'n_DP': '', 't_depth_sample': '653', 't_ref_count_sample': '511', 't_alt_count_sample': '142', 'is_fillout': 'False'},

        {'Hugo_Symbol': 'RTEL1', 'Entrez_Gene_Id': '51750', 'Center': '.', 'NCBI_Build': 'GRCh37', 'Chromosome': '20', 'Start_Position': '62321135', 'End_Position': '62321135', 'Strand': '+', 'Variant_Classification': 'Silent', 'Variant_Type': 'SNP', 'Reference_Allele': 'G', 'Tumor_Seq_Allele1': 'G', 'Tumor_Seq_Allele2': 'A', 'dbSNP_RS': 'rs746824222', 'dbSNP_Val_Status': '', 'Tumor_Sample_Barcode': 'Sample23', 'Matched_Norm_Sample_Barcode': 'Sample23-N', 'Match_Norm_Seq_Allele1': 'G', 'Match_Norm_Seq_Allele2': 'G', 'Tumor_Validation_Allele1': '', 'Tumor_Validation_Allele2': '', 'Match_Norm_Validation_Allele1': '', 'Match_Norm_Validation_Allele2': '', 'Verification_Status': '', 'Validation_Status': '', 'Mutation_Status': '', 'Sequencing_Phase': '', 'Sequence_Source': '', 'Validation_Method': '', 'Score': '', 'BAM_File': '', 'Sequencer': '', 'Tumor_Sample_UUID': '', 'Matched_Norm_Sample_UUID': '', 'HGVSc': 'c.2130G>A', 'HGVSp': 'p.=', 'HGVSp_Short': 'p.Q710=', 'Transcript_ID': 'ENST00000508582', 'Exon_Number': '24/35', 't_depth': '184', 't_ref_count': '184', 't_alt_count': '0', 'n_depth': '', 'n_ref_count': '', 'n_alt_count': '', 'all_effects': 'RTEL1,synonymous_variant,p.=,ENST00000318100,;RTEL1,synonymous_variant,p.=,ENST00000370018,NM_032957.4,NM_016434.3;RTEL1,synonymous_variant,p.=,ENST00000360203,NM_001283009.1;RTEL1,synonymous_variant,p.=,ENST00000508582,;RTEL1,synonymous_variant,p.=,ENST00000425905,;RTEL1,upstream_gene_variant,,ENST00000370003,;RTEL1-TNFRSF6B,synonymous_variant,p.=,ENST00000482936,;RTEL1-TNFRSF6B,synonymous_variant,p.=,ENST00000492259,;RTEL1-TNFRSF6B,non_coding_transcript_exon_variant,,ENST00000480273,;RTEL1-TNFRSF6B,non_coding_transcript_exon_variant,,ENST00000496281,;RTEL1,upstream_gene_variant,,ENST00000496816,;', 'Allele': 'A', 'Gene': 'ENSG00000258366', 'Feature': 'ENST00000508582', 'Feature_type': 'Transcript', 'Consequence': 'synonymous_variant', 'cDNA_position': '2476/4273', 'CDS_position': '2130/3732', 'Protein_position': '710/1243', 'Amino_acids': 'Q', 'Codons': 'caG/caA', 'Existing_variation': 'rs746824222', 'ALLELE_NUM': '1', 'DISTANCE': '', 'STRAND_VEP': '1', 'SYMBOL': 'RTEL1', 'SYMBOL_SOURCE': 'HGNC', 'HGNC_ID': '15888', 'BIOTYPE': 'protein_coding', 'CANONICAL': 'YES', 'CCDS': 'CCDS13530.3', 'ENSP': 'ENSP00000424307', 'SWISSPROT': 'Q9NZ71', 'TREMBL': '', 'UNIPARC': 'UPI00019B2219', 'RefSeq': '', 'SIFT': '', 'PolyPhen': '', 'EXON': '24/35', 'INTRON': '', 'DOMAINS': 'Superfamily_domains:SSF52540,SMART_domains:SM00491,Pfam_domain:PF13307,TIGRFAM_domain:TIGR00604,hmmpanther:PTHR11472:SF4,hmmpanther:PTHR11472', 'AF': '', 'AFR_AF': '', 'AMR_AF': '', 'ASN_AF': '', 'EAS_AF': '', 'EUR_AF': '', 'SAS_AF': '', 'AA_AF': '', 'EA_AF': '', 'CLIN_SIG': '', 'SOMATIC': '', 'PUBMED': '', 'MOTIF_NAME': '', 'MOTIF_POS': '', 'HIGH_INF_POS': '', 'MOTIF_SCORE_CHANGE': '', 'IMPACT': 'LOW', 'PICK': '1', 'VARIANT_CLASS': 'SNV', 'TSL': '', 'HGVS_OFFSET': '', 'PHENO': '', 'MINIMISED': '', 'ExAC_AF': '', 'ExAC_AF_AFR': '', 'ExAC_AF_AMR': '', 'ExAC_AF_EAS': '', 'ExAC_AF_FIN': '', 'ExAC_AF_NFE': '', 'ExAC_AF_OTH': '', 'ExAC_AF_SAS': '', 'GENE_PHENO': '1', 'FILTER': '.', 'flanking_bps': 'AGG', 'vcf_id': '.', 'vcf_qual': '.', 'ExAC_AF_Adj': '', 'ExAC_AC_AN_Adj': '', 'ExAC_AC_AN': '', 'ExAC_AC_AN_AFR': '', 'ExAC_AC_AN_AMR': '', 'ExAC_AC_AN_EAS': '', 'ExAC_AC_AN_FIN': '', 'ExAC_AC_AN_NFE': '', 'ExAC_AC_AN_OTH': '', 'ExAC_AC_AN_SAS': '', 'ExAC_FILTER': '', 'gnomAD_AF': '', 'gnomAD_AFR_AF': '', 'gnomAD_AMR_AF': '', 'gnomAD_ASJ_AF': '', 'gnomAD_EAS_AF': '', 'gnomAD_FIN_AF': '', 'gnomAD_NFE_AF': '', 'gnomAD_OTH_AF': '', 'gnomAD_SAS_AF': '', 'vcf_pos': '62321135', 'AC': '1', 'AN': '2', 'SRC': 'Sample24,', 't_GT': './.', 'n_GT': '', 't_FL_AD': '0', 'n_FL_AD': '', 't_FL_ADN': '0', 'n_FL_ADN': '', 't_FL_ADP': '0', 'n_FL_ADP': '', 't_FL_DP': '184', 'n_FL_DP': '', 't_FL_DPN': '95', 'n_FL_DPN': '', 't_FL_DPP': '89', 'n_FL_DPP': '', 't_FL_RD': '184', 'n_FL_RD': '', 't_FL_RDN': '95', 'n_FL_RDN': '', 't_FL_RDP': '89', 'n_FL_RDP': '', 't_FL_VF': '0', 'n_FL_VF': '', 't_AD': '', 'n_AD': '', 't_DP': '.', 'n_DP': '', 't_depth_sample': '', 't_ref_count_sample': '', 't_alt_count_sample': '', 'is_fillout': 'True'}
        ]

        self.assertEqual(records, expected_records)

    def test_run_fillout_workflow_small_2(self):
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


    def test_run_fillout_workflow_small_3(self):
        """
        Test case with small variant set that should have filters applied inside the pipeline

        run the fillout like this:

        sample_id       sample_type
        Sample3 clinical
        Sample5 research
        Sample2 research
        Sample1 research
        Sample4 clinical
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
        self.assertEqual(len(mutations), 121699)
        hash = md5_obj(mutations)
        expected_hash = '4b55be411af84915ab07e01bb04d0619'
        self.assertEqual(hash, expected_hash)

        comments, mutations = self.load_mutations(filtered_output_path, strip = True)
        self.assertEqual(len(mutations), 121699)
        hash = md5_obj(mutations)
        expected_hash = '4b55be411af84915ab07e01bb04d0619'
        self.assertEqual(hash, expected_hash)






    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_fillout_workflow2(self):
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
        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$be8534bcaf326de029790a832ab5b44a17a03d22',
                # 'size': 40194610,
                'path':  output_path
                }
            }
        # NOTE: for some reason, this file keeps coming out with different annotations for 'splice_acceptor_variant' or `splice_donor_variant`
        # this keeps changing the byte size and checksum so need to remove those here for now
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(output_path)

        self.assertEqual(len(mutations), 38920)

        # Need to remove these fields because they are inconsistent on the output maf file;
        for mut in mutations:
            mut.pop('all_effects')
            mut.pop('Consequence')
            mut.pop('Variant_Classification')

        hash = md5_obj(mutations)
        expected_hash = 'bc5f54f1057a7ba29f55d9d4aac92a01'
        self.assertEqual(hash, expected_hash)




    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_run_fillout_workflow3(self):
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

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.maf',
                'class': 'File',
                # 'checksum': 'sha1$2f60f58389ec65af87612c7532ad28b882fb84ba',
                # 'size': 26238820,
                'path':  output_path
                }
            }
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        self.assertCWLDictEqual(output_json, expected_output)

        comments, mutations = self.load_mutations(output_path)

        hash = md5_obj(mutations)

        self.assertEqual(len(mutations), 26404)

        # Need to remove these fields because they are inconsistent on the output maf file;
        for mut in mutations:
            mut.pop('all_effects')
            mut.pop('Consequence')
            mut.pop('Variant_Classification')

        hash = md5_obj(mutations)
        expected_hash = '4a03d128d76b72328b62a87814d89993'
        self.assertEqual(hash, expected_hash)



if __name__ == "__main__":
    unittest.main()
