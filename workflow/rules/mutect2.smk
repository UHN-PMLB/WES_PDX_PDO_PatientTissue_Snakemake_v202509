rule gatk_mutect2:
    input:
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai"
    output:
        vcf = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz.tbi"
    params:
        ref = config["ref_index"]["genome"],
        pon = config["ref"]["pon"],
        gnomad = config["ref"]["gnomad"]
    threads: 16
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx32G" Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            --native-pair-hmm-threads {threads} \
            --panel-of-normals {params.pon} \
            --germline-resource {params.gnomad} \
            -O {output.vcf}
        """

rule gatk_filtermutectcalls:
    input:
        vcf = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz.tbi",
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai"
    output:
        vcf = "results/mutect2/{strand}/{sample}.filtered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.filtered.vcf.gz.tbi"
    params:
        ref = config["ref_index"]["genome"]
    threads: 8
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx16G" FilterMutectCalls \
            -R {params.ref} \
            -V {input.vcf} \
            -O {output.vcf}
        """

rule vcf2maf:
    input:
        vcf = "results/mutect2/{strand}/{sample}.filtered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.filtered.vcf.gz.tbi",
        ref = config["ref_index"]["genome"]
    output:
        maf = "results/maf/{strand}/{sample}.maf"
    threads: 4
    shell:
        """
        module load vcf2maf/1.6.22
        module load samtools/1.20
        module load vep/113

        mkdir -p results/maf/{wildcards.strand}

        # Step 1: keep only PASS variants
        bcftools view -f PASS {input.vcf} -Oz -o results/maf/{wildcards.strand}/{wildcards.sample}.pass.vcf.gz
        tabix -p vcf results/maf/{wildcards.strand}/{wildcards.sample}.pass.vcf.gz

        # Step 2: convert to MAF using vcf2maf (with ExAC filter)
        vcf2maf.pl \
            --input-vcf results/maf/{wildcards.strand}/{wildcards.sample}.pass.vcf.gz \
            --output-maf {output.maf} \
            --tumor-id {wildcards.sample} \
            --normal-id NONE \
            --ref-fasta {input.ref} \
            --ncbi-build GRCh38 \
            --species homo_sapiens \
            --vep-forks {threads} \
            --vep-path $(dirname $(which vep)) \
            --vep-data /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/VEP_cache \
            --filter-vcf /cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz
        """

rule merge_mafs:
    input:
        expand("results/maf/{strand}/{sample}.maf",
               strand=["pe","se"],  # adjust if only one type
               sample=SAMPLES)      # your sample list from config or snakefile
    output:
        maf = "results/maf/cohort_merged.maf"
    shell:
        """
        module load perl

        # Concatenate header from the first file + all non-header lines
        head -n 1 {input[0]} > {output.maf}
        for f in {input}; do
            tail -n +2 $f >> {output.maf}
        done
        """
