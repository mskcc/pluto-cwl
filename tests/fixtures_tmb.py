#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common location for some Python objects used in some test cases
"""
comments = [
['# comment 1'],
['# comment 2']
]

row1 = {
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
}
row1_2 = {
't_af': '0.51',
't_depth': '551',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
}

row2 = {
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'splice_region_variant,synonymous_variant' # exclude due to synonymous_variant
}

row3 = { # this one should pass filter
't_af': '0.50',
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}
row3_2 = { # this one should pass filter
't_af': '0.52',
't_depth': '552',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}
row3_3 = {
't_af': '0.01', # exclude due to low AF
't_depth': '552',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}


row4 = {
't_af': '0.01', # exclude due to low AF
't_depth': '550',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}

row5 = {
't_af': '0.51',
't_depth': '90', # exclude due to low coverage
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'missense_variant'
}

row6 = { # this one should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'EGFR',
'Start_Position': '1',
'Consequence': 'splice_region_variant'
}

row7 = { # this one should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295340', # good value; is_TERT_promoter = True
'Consequence': 'splice_region_variant'
}

row8 = { # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295339', # good value; is_TERT_promoter = True
'Consequence': 'splice_region_variant'
}

row9 = { # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295341', # bad value; is_TERT_promoter = False
'Consequence': 'splice_region_variant' # include anyway because its not synonymous_variant
}

row10 ={ # this should pass filter
't_af': '0.45',
't_depth': '590',
'Hugo_Symbol': 'TERT',
'Start_Position': '1295339', # good value; is_TERT_promoter = True
'Consequence': 'synonymous_variant' # include even though its synonymous_variant
}


# PASS: 6
# FAIL: 4
rows1 = [
    row1, # exclude due to synonymous_variant
    row2, # exclude due to synonymous_variant
    row3, # this one should pass filter
    row4, # exclude due to low AF
    row5, # this one should pass filter
    row6, # this one should pass filter
    row7,  # this one should pass filter
    row8, # this should pass filter
    row9, # this should pass filter
    row10 # this should pass filter
]


# PASS: 5
# FAIL: 5
rows2 = [
    row1_2, # exclude due to synonymous_variant
    row2, # exclude due to synonymous_variant
    row3_2, # this one should pass filter
    row3_3, # exclude due to low AF
    row4, # exclude due to low AF
    row5, # this one should pass filter
    row6, # this one should pass filter
    row8, # this should pass filter
    row9, # this should pass filter
    row10 # this should pass filter
]
