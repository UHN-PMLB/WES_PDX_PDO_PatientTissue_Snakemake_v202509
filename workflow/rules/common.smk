import glob
import pandas as pd
from snakemake.utils import validate
#from snakemake.remote import FTP
#ftp = FTP.RemoteProvider()
configfile: "config/config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

def get_sample_type(wildcards):
    """Get sample type (PDX, PDO, patient_tissue) from unit sheet."""
    return samples.loc[wildcards.sample, "sample_type"]


def get_fq1(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_paired_end(wildcards.sample):
        u = samples.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
    else:
        u = samples.loc[ (wildcards.sample), ["fq1"] ].dropna()
    return [ f"{u.fq1}" ]

def get_fq2(wildcards):
    """Get raw FASTQ files from unit sheet."""
    if is_paired_end(wildcards.sample):
        u = samples.loc[ (wildcards.sample), ["fq1", "fq2"] ].dropna()
        fq2 = [ f"{u.fq2}" ]
    else:
        fq2 = ""
    return fq2

def get_rsem_output_all_units(level="genes"):
    return [
        f"results/rsem/{sample.sample_name}.{level}.results"
        for sample in samples.itertuples()
    ]

def is_paired_end(sample):
    sample_units = samples.loc[sample]
    fq2_null = sample_units.isnull()["fq2"]
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired
    
def get_star_output_all_units(wildcards, fi="coord", use='all'):
    if fi == "coord":
        outfile = "Aligned.sortedByCoord.out.bam"
    elif fi == 'transcriptome':
        outfile = "Aligned.toTranscriptome.out.bam"
    else:
	    outfile = "ReadsPerGene.out.tab"
    res = []
    if use == 'all':
        for sample in samples.itertuples():
            if is_paired_end(sample.sample_name):
                lib = "pe"
            else:
                lib = "se"
            res.append(
                "results/star/{}/{}/{}".format(
                    lib, sample.sample_name, outfile
                )
            )
    else:
        if is_paired_end(wildcards.sample):
            lib = 'pe'
        else:
            lib = 'se'
        res.append("results/star/{}/{}/{}".format(
                lib, wildcards.sample, outfile
                )
        )
    return res

def get_star_output(wildcards, fi="coord", bai=False):
    if fi == "coord":
        outfile = "Aligned.sortedByCoord.out.bam"
    elif fi == 'transcriptome':
        outfile = "Aligned.toTranscriptome.out.bam"
    else:
        outfile = "ReadsPerGene.out.tab"
    if bai:
        outfile = outfile + ".bai"
        
    if is_paired_end(wildcards.sample):
        lib = 'pe'
    else:
        lib = 'se'
    
    res = []
    res.append("results/star/{}/{}/{}".format(
        lib, wildcards.sample, outfile
        )
    )
    return res
    

