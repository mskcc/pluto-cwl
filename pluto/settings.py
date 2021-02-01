"""
Put settings to use for the tests in here for easier access
"""
import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CWL_DIR = os.path.join(os.path.dirname(THIS_DIR), "cwl")
REF_DIR = os.path.join(os.path.dirname(THIS_DIR), "ref")
EXAMPLES_DIR = os.path.join(os.path.dirname(THIS_DIR), "examples")

# common args to be included in all cwltool invocations
CWL_ARGS = [
    "--preserve-environment", "PATH",
    "--preserve-environment", "SINGULARITY_CACHEDIR",
    "--singularity"
]
TOIL_ARGS = [
    '--singularity', # run with Singularity instead of Docker
    '--batchSystem', 'lsf', '--disableCaching', 'True', # use LSF on the HPC to submit jobs
    '--disable-user-provenance', '--disable-host-provenance',
    '--preserve-entire-environment', # need to propagate the env vars for Singularity, etc., into the HPC jobs
    '--retryCount', '1',
    '--maxLocalJobs', '500', # run up to 500 jobs at once; not actually "local", this includes HPC jobs
    '--statePollingWait', '10', # check available jobs every 10 seconds instead of after every job is submitted
    '--clean', 'onSuccess'
]

# location on the filesystem for static fixtures
FIXTURES_DIR = os.environ.get('FIXTURES_DIR', '/juno/work/ci/helix_filters_01/fixtures')
FACETS_SNPS_VCF = os.environ.get('FACETS_SNPS_FILE', '/juno/work/ci/resources/genomes/GRCh37/facets_snps/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf')
KNOWN_FUSIONS_FILE = os.path.join(REF_DIR, "known_fusions_at_mskcc.txt")
IMPACT_FILE=os.environ.get('IMPACT_file', '/work/ci/helix_filters_01/reference_data/gene_lists/all_IMPACT_genes.tsv')

ARGOS_VERSION_STRING = os.environ.get('ARGOS_VERSION_STRING', '2.x') # TODO: deprecate this
IS_IMPACT = os.environ.get('IS_IMPACT', "True") # TODO: deprecate this
PORTAL_FILE = os.environ.get('PORTAL_FILE', 'data_mutations_extended.txt') # TODO: deprecate this
PORTAL_CNA_FILE = os.environ.get('PORTAL_CNA_FILE', 'data_CNA.txt') # TODO: deprecate this


DATA_SETS = {
    "Proj_08390_G": {
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
        "SNP_PILEUP_DIR": os.path.join(FIXTURES_DIR, "Proj_08390_G", "snp-pileup")
    }
}
