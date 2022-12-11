#!/usr/bin/env python
"""
"""
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REF_DIR = os.path.join(os.path.dirname(THIS_DIR), "ref") # ../ref

FIXTURES_DIR = os.environ.get('FIXTURES_DIR', '/juno/work/ci/helix_filters_01/fixtures')
FACETS_SNPS_VCF = os.environ.get('FACETS_SNPS_FILE', '/juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf')
KNOWN_FUSIONS_FILE = os.path.join(REF_DIR, "known_fusions_at_mskcc.txt")
IMPACT_FILE=os.environ.get('IMPACT_file', '/work/ci/helix_filters_01/reference_data/gene_lists/all_IMPACT_genes.tsv')
CONPAIR_MARKERS_BED = os.environ.get("CONPAIR_MARKERS_BED", "/juno/work/ci/helix_filters_01/reference_data/concordance/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.bed")
CONPAIR_MARKERS_TXT = os.environ.get("CONPAIR_MARKERS_TXT", "/juno/work/ci/helix_filters_01/reference_data/concordance/markers/IMPACT468/FP_tiling_genotypes_for_Conpair.txt")

ARGOS_VERSION_STRING = os.environ.get('ARGOS_VERSION_STRING', '2.x') # TODO: deprecate this
IS_IMPACT = os.environ.get('IS_IMPACT', "True") # TODO: deprecate this
PORTAL_FILE = os.environ.get('PORTAL_FILE', 'data_mutations_extended.txt') # TODO: deprecate this
PORTAL_CNA_FILE = os.environ.get('PORTAL_CNA_FILE', 'data_CNA.txt') # TODO: deprecate this

REF_FASTA = os.environ.get('REF_FASTA', '/juno/work/ci/resources/genomes/GRCh37/fasta/b37.fasta')
MICROSATELLITES_LIST = os.environ.get("MICROSATELLITES_LIST", "/work/ci/resources/request_files/msisensor/microsatellites.list")
# $ md5sum /work/ci/resources/request_files/msisensor/microsatellites.list
# dc982a3bfe1e33b201b99a8ebf3acd61  /work/ci/resources/request_files/msisensor/microsatellites.list
# $ wc -l /work/ci/resources/request_files/msisensor/microsatellites.list
# 33422661 /work/ci/resources/request_files/msisensor/microsatellites.list



DATA_SETS = {
    "Proj_08390_G": { # full sample Argos output
        "DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G"),
        "MAF_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "maf"),
        "BAM_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "bam"),
        # "SNP_PILEUP_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "snp_pileup"),
        "FACETS_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "facets"),
        "FACETS_SUITE_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "facets-suite"),
        "INPUTS_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "inputs"),
        "QC_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "qc"),
        "targets_list": "/juno/work/ci/resources/roslin_resources/targets/HemePACT_v4/b37/HemePACT_v4_b37_targets.ilist",
        "analyst_file": "Proj_08390_G.muts.maf", # TODO: deprecate this
        "analysis_gene_cna_file": "Proj_08390_G.gene.cna.txt", # TODO: deprecate this
        "MAF_FILTER_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "maf_filter"),
        "SNP_PILEUP_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "snp-pileup"),
        'REF_FASTA': REF_FASTA,
        'microsatellites_file': MICROSATELLITES_LIST,
        "MSI_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "msi"),
        "TMB_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "tmb"),
    },
    "Proj_1": { # same as Proj_08390_G but both filenames and file contents have been scrubbed; results in different file md5's
        "MAF_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "maf"),
        "BAM_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "bam"),
        "FACETS_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "facets"),
        "FACETS_SUITE_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "facets-suite"),
        "QC_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "qc"),
        "INPUTS_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "inputs"),
        'REF_FASTA': REF_FASTA,
        "targets_list": "/juno/work/ci/resources/roslin_resources/targets/HemePACT_v4/b37/HemePACT_v4_b37_targets.ilist",
        "MSI_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "msi"),
        "TMB_DIR": os.path.join(FIXTURES_DIR, "Proj_1", "tmb"),
    },
    "demo":{ # small subset of samples on a full project
        "DIR": os.path.join(FIXTURES_DIR, "demo"),
        "MAF_DIR": os.path.join(FIXTURES_DIR, "demo", "maf"),
        "BAM_DIR": os.path.join(FIXTURES_DIR, "demo", "bam"),
        "QC_DIR": os.path.join(FIXTURES_DIR, "demo", "qc"),
        "INPUTS_DIR": os.path.join(FIXTURES_DIR, "demo", "inputs"),
        "SNP_PILEUP_DIR": os.path.join(FIXTURES_DIR, "demo", "snp-pileup"),
        "FACETS_DIR": os.path.join(FIXTURES_DIR, "demo", "facets"),
        "targets_list": "/juno/work/ci/resources/roslin_resources/targets/HemePACT_v4/b37/HemePACT_v4_b37_targets.ilist",
        'microsatellites_file': os.path.join(FIXTURES_DIR, "demo", "microsatellites", 'microsatellites.head500000.list'),
        # $ md5sum microsatellites.head500000.list
        # aa0126e6a916ec82a2837989458918b3  microsatellites.head500000.list
        'REF_FASTA': REF_FASTA
    },
    # dataset selected for use with fillout since it has pooled normals
    "07618_AG": {
        "DIR": os.path.join(FIXTURES_DIR, "07618_AG"),
        "BAM_DIR": os.path.join(FIXTURES_DIR, "07618_AG", "bam"),
        "MAF_DIR": os.path.join(FIXTURES_DIR, "07618_AG", "maf")
    },
    "Fillout01": {
        "DIR": os.path.join(FIXTURES_DIR, "Fillout01"),
        "BAM_DIR": os.path.join(FIXTURES_DIR, "Fillout01", "bam"),
        "MAF_DIR": os.path.join(FIXTURES_DIR, "Fillout01", "maf"),
        "VCF_DIR": os.path.join(FIXTURES_DIR, "Fillout01", "vcf"),
        "OUTPUT_DIR": os.path.join(FIXTURES_DIR, "Fillout01", "output")
    },
    "Conpair_1": {
        "BAM_DIR": os.path.join(FIXTURES_DIR, "Conpair_1", "bam"),
        "LIKELIHOODS": os.path.join(FIXTURES_DIR, "Conpair_1", "likelihoods"),
    }
}
