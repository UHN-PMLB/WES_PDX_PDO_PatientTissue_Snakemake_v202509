rule vep_annotate_unfiltered:
    input:
        vcf = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz"
    output:
        vcf_annotated = "results/mutect2/{strand}/{sample}.unfiltered.vep.vcf.gz"
    params:
        cache = "/cluster/tools/data/commondata/ensembl/vep/113",
        fasta = config["ref_index"]["genome"]
    threads: 8
    shell:
        """
        module load vep/113

        vep \
          -i {input.vcf} \
          -o {output.vcf_annotated} \
          --vcf \
          --assembly GRCh38 \
          --offline \
          --cache \
          --dir_cache {params.cache} \
          --fasta {params.fasta} \
          --symbol \
          --canonical \
          --everything \
          --fork {threads} \
          --force_overwrite
        """

rule rescue_driver_variants:
    input:
        vcf_annotated = "results/mutect2/{strand}/{sample}.unfiltered.vep.vcf.gz",
        drivers = "ref/cosmic_cgc/cgc_drivers_tier1_v102_GRCh38.txt"
    output:
        rescue_vcf = "results/mutect2/{strand}/{sample}.rescue.candidates.vcf"
    threads: 4
    shell:
        """
        python workflow/scripts/rescue_driver_variants.py \
            {input.vcf_annotated} {input.drivers} {output.rescue_vcf}
        """


rule vcf2maf_rescue:
    input:
        vcf = "results/mutect2/{strand}/{sample}.rescue.candidates.vcf",
        ref = config["ref_index"]["genome"]
    output:
        maf = "results/maf/{strand}/{sample}.rescue.raw.maf"
    params:
        vep_cache = "/cluster/tools/data/commondata/ensembl/vep/113"
    threads: 4
    shell:
        """
        module load vcf2maf/1.6.22
        module load samtools/1.20
        module load vep/113

        # Explicitly specify VEP executable directory
        VEP_PATH=/cluster/tools/software/rocky9/vep/113

        vcf2maf.pl \
          --input-vcf {input.vcf} \
          --output-maf {output.maf} \
          --tumor-id {wildcards.sample} \
          --normal-id NONE \
          --ref-fasta {input.ref} \
          --ncbi-build GRCh38 \
          --species homo_sapiens \
          --vep-path $VEP_PATH \
          --vep-data {params.vep_cache}
        """

rule filter_rescue_maf:
    input:
        maf = "results/maf/{strand}/{sample}.rescue.raw.maf"
    output:
        maf = "results/maf/{strand}/{sample}.rescue.filtered.maf"
    shell:
        """
        python workflow/scripts/filter_rescue_maf.py {input.maf} {output.maf}
        """
