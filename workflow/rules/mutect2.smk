rule mutect2:
    input:
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai"
    output:
        vcf = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz.tbi"
    params:
        ref = config["ref_index"]["genome"],
        pon = config["ref"]["pon"],
        gnomad = config["ref"]["gnomad"],
        intervals = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.bed"
    threads: 16
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx32G" Mutect2 \
            -R {params.ref} \
            -I {input.bam} \
            -L {params.intervals} \
            --native-pair-hmm-threads {threads} \
            --panel-of-normals {params.pon} \
            --germline-resource {params.gnomad} \
            -O {output.vcf}
        """


rule filtermutectcalls:
    input:
        vcf = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz",
        vcf_index = "results/mutect2/{strand}/{sample}.unfiltered.vcf.gz.tbi",
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai",
        contamination = "results/mutect2/{strand}/{sample}.contamination.table"
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
            --contamination-table {input.contamination} \
            -O {output.vcf}
        """


rule filter_pass_variants:
    input:
        vcf = "results/mutect2/{strand}/{sample}.filtered.vcf.gz"
    output:
        pass_vcf = "results/mutect2/{strand}/{sample}.pass_rare.vcf"
    params:
        gnomad = config["ref"]["gnomad"]
    threads: 2
    shell:
        r"""
        module load samtools/1.20

        bcftools annotate -a {params.gnomad} -c INFO/AF {input.vcf} \
          | bcftools view -i 'FILTER="PASS" && (INFO/AF<0.01 || INFO/AF=".")' \
          -Ov -o {output.pass_vcf}
        """


###############################################
#     Standard strict (filtered) variants     #
###############################################
rule vcf2maf_strict:
    input:
        vcf = "results/mutect2/{strand}/{sample}.pass_rare.vcf",
        ref = config["ref_index"]["genome"]
    output:
        maf = "results/maf/{strand}/{sample}.strict.maf"
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


###############################################
#   Merge strict + rescue variants per sample #
###############################################
rule merge_maf_final:
    input:
        strict = "results/maf/{strand}/{sample}.strict.maf",
        rescue = "results/maf/{strand}/{sample}.rescue.filtered.maf"
    output:
        maf = "results/maf/{strand}/{sample}.final.maf"
    shell:
        """
        mkdir -p $(dirname {output.maf})

        head -n 1 {input.strict} > {output.maf}
        tail -n +2 {input.strict} >> {output.maf}
        tail -n +2 {input.rescue} >> {output.maf}

        # Deduplicate by Chrom:Start:Alt
        awk -F'\\t' 'NR==1{{print; next}} NR>1{{key=$2":"$3":"$7":"$8; if(!(key in seen)){{print; seen[key]=1}}}}' {output.maf} > {output.maf}.tmp && mv {output.maf}.tmp {output.maf}
        """


###############################################
#   Merge all samples into cohort MAF         #
###############################################
rule merge_mafs:
    input:
        expand("results/maf/pe/{sample}.final.maf", sample=samples["sample_name"])
    output:
        maf = "results/maf/cohort_merged_final.maf"
    shell:
        """
        head -n 1 {input[0]} > {output.maf}
        for f in {input}; do
            tail -n +2 $f >> {output.maf}
        done
        """
