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
