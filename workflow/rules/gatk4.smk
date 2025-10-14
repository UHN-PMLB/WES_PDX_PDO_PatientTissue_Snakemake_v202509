rule gatk_markduplicates:
    input:
        bam = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bam = "results/gatk/{strand}/{sample}.dedup.bam",
        bai = "results/gatk/{strand}/{sample}.dedup.bam.bai",
        metrics = "results/gatk/{strand}/{sample}.dedup.metrics.txt"
    threads: 8
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx32G" MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            ----create-output-bam-index true \
            --tmp-dir tmp/ \
            --spark-master local[{threads}]
        """


rule gatk_bqsr:
    input:
        bam = "results/gatk/{strand}/{sample}.dedup.bam",
        bai = "results/gatk/{strand}/{sample}.dedup.bai"
    output:
        table = "results/gatk/{strand}/{sample}.recal.table"
    params:
        ref = config["ref_index"]["genome"],
        known_sites = expand("--known-sites {site}", site=config["ref"]["known_sites"]) # list of dbSNP + Mills/1000G VCFs
    threads: 8
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx32G" BaseRecalibrator \
            -R {params.ref} \
            -I {input.bam} \
            {params.known_sites} \
            -O {output.table}
        """


rule gatk_applybqsr:
    input:
        bam = "results/gatk/{strand}/{sample}.dedup.bam",
        bai = "results/gatk/{strand}/{sample}.dedup.bai",
        table = "results/gatk/{strand}/{sample}.recal.table"
    output:
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai"
    params:
        ref = config["ref_index"]["genome"]
    threads: 8
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx32G" ApplyBQSR \
            -R {params.ref} \
            -I {input.bam} \
            --bqsr-recal-file {input.table} \
            -O {output.bam}
        """
