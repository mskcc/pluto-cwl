#!/usr/bin/env python3

# these are the expected headers on the cBioPortal formatted data_clinical_sample.txt file
# values shown here are transposed for clarity
expected_data_clinical_sample_columns = [
['SAMPLE_ID', 'SAMPLE_ID', 'STRING', '1'],
['IGO_ID', 'IGO_ID', 'STRING', '1'],
['PATIENT_ID', 'PATIENT_ID', 'STRING', '1'],
['COLLAB_ID', 'COLLAB_ID', 'STRING', '0'],
['SAMPLE_TYPE', 'SAMPLE_TYPE', 'STRING', '1'],
['SAMPLE_CLASS', 'SAMPLE_CLASS', 'STRING', '1'],
['GENE_PANEL', 'GENE_PANEL', 'STRING', '1'],
['ONCOTREE_CODE', 'ONCOTREE_CODE', 'STRING', '1'],
['SPECIMEN_PRESERVATION_TYPE', 'SPECIMEN_PRESERVATION_TYPE', 'STRING', '1'],
['TISSUE_SITE', 'TISSUE_SITE', 'STRING', '1'],
['REQUEST_ID', 'REQUEST_ID', 'STRING', '1'],
['PROJECT_ID', 'PROJECT_ID', 'STRING', '1'],
['PIPELINE', 'PIPELINE', 'STRING', '1'],
['PIPELINE_VERSION', 'PIPELINE_VERSION', 'STRING', '1'],
['SAMPLE_COVERAGE', 'SAMPLE_COVERAGE', 'NUMBER', '1'],
['PROJECT_PI', 'PROJECT_PI', 'STRING', '1'],
['REQUEST_PI', 'REQUEST_PI', 'STRING', '1'],
['genome_doubled', 'genome_doubled', 'STRING', '0'],
['ASCN_PURITY', 'ASCN_PURITY', 'NUMBER', '1'],
['ASCN_PLOIDY', 'ASCN_PLOIDY', 'NUMBER', '1'],
['ASCN_VERSION', 'ASCN_VERSION', 'STRING', '0'],
['ASCN_WGD', 'ASCN_WGD', 'STRING', '1'],
['CMO_MSI_SCORE', 'CMO_MSI_SCORE', 'NUMBER', '0'],
['CMO_MSI_STATUS', 'CMO_MSI_STATUS', 'STRING', '0'],
['CMO_TMB_SCORE', 'CMO_TMB_SCORE', 'NUMBER', '1'],
]
# TODO: add CMO_ASSAY_COVERAGE


# transpose the values again so they match what is read in the file
expected_data_clinical_sample_comments = list(map(list, zip(*expected_data_clinical_sample_columns)))

# looks like this;
# expected_comments = [
# ['SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD', 'CMO_TMB_SCORE', 'CMO_MSI_SCORE', 'CMO_MSI_STATUS'],
# ['SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD', 'CMO_TMB_SCORE', 'CMO_MSI_SCORE', 'CMO_MSI_STATUS'],
# ['STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'NUMBER', 'STRING', 'STRING', 'STRING', 'NUMBER', 'NUMBER', 'STRING', 'STRING', 'NUMBER', 'NUMBER', 'STRING'],
# ['1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1', '1', '0', '0']
# ]






# some example maf rows for test cases
maf_row1 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
}
maf_row2 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
}
maf_row3 = { # extra row with no match in facets
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865788",
"End_Position": "99865789",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo" # dummy value that would only be in portal data_mutations_extended.txt output maf file
}
maf_row4 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample2-T", # different sample here !!
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo"
}
maf_row5 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample2-T", # different sample here !!
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo"
}
facets_row1 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"ASCN.TOTAL_COPY_NUMBER": "1" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
}
facets_row2 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"ASCN.TOTAL_COPY_NUMBER": "2" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
}
facets_row3 = { # extra row with no match in maf
"Hugo_Symbol" : "PNISR2",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"ASCN.TOTAL_COPY_NUMBER": "2" # dummy value that would only be in facets Tumor1.Normal1_hisens.ccf.maf output maf file
}
facets_row4 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample2-T", # different sample !!
"Matched_Norm_Sample_Barcode": "Sample2-N",
"ASCN.TOTAL_COPY_NUMBER": "5"
}
facets_row5 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample2-T", # different sample !!
"Matched_Norm_Sample_Barcode": "Sample2-N",
"ASCN.TOTAL_COPY_NUMBER": "6"
}

demo_comments = [
['# comment 1'],
['# comment 2']
]

expected_row1 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo",
"ASCN.CLONAL": "1"
}
expected_row2 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo",
"ASCN.CLONAL": "2"
}
expected_row3 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865788",
"End_Position": "99865789",
"Tumor_Sample_Barcode": "Sample1-T",
"Matched_Norm_Sample_Barcode": "Sample1-N",
"portal_val": "foo",
"ASCN.CLONAL": "."
}
expected_row4 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample2-T",
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo",
"ASCN.CLONAL": "5"
}
expected_row4_2 = {
"Hugo_Symbol" : "FGF3",
"Entrez_Gene_Id" : "2248",
"Chromosome" : "11",
"Start_Position" : "69625447",
"End_Position": "69625448",
"Tumor_Sample_Barcode": "Sample2-T",
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo",
"ASCN.CLONAL": "." # no matching sample, default empty value
}
expected_row5 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample2-T",
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo",
"ASCN.CLONAL": "6"
}
expected_row5_2 = {
"Hugo_Symbol" : "PNISR",
"Entrez_Gene_Id" : "25957",
"Chromosome" : "6",
"Start_Position" : "99865784",
"End_Position": "99865785",
"Tumor_Sample_Barcode": "Sample2-T",
"Matched_Norm_Sample_Barcode": "Sample2-N",
"portal_val": "foo",
"ASCN.CLONAL": "." # no matching sample, default empty value
}
