rule gatk_markduplicates:
    input:
        bam = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam",
        bai = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        bam = temp("results/gatk/{strand}/{sample}.dedup.bam"),
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
            --create-output-bam-index true \
            --tmp-dir tmp/ \
            --spark-master local[{threads}]
        """

rule gatk_bqsr:
    input:
        bam = "results/gatk/{strand}/{sample}.dedup.bam",
        bai = "results/gatk/{strand}/{sample}.dedup.bam.bai"
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
        bai = "results/gatk/{strand}/{sample}.dedup.bam.bai",
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

rule gatk_getpileupsummaries:
    input:
        bam = "results/gatk/{strand}/{sample}.recal.bam",
        bai = "results/gatk/{strand}/{sample}.recal.bai",
        vcf = "/cluster/tools/data/genomes/human/GRCh38/iGenomes/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    output:
        pileups = "results/mutect2/{strand}/{sample}.pileups.table"
    params:
        ref = config["ref_index"]["genome"]
    threads: 4
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx56G" GetPileupSummaries \
            -I {input.bam} \
            -V {input.vcf} \
            -L {input.vcf} \
            -R {params.ref} \
            -O {output.pileups}
        """

rule gatk_calculatecontamination:
    input:
        pileups = "results/mutect2/{strand}/{sample}.pileups.table"
    output:
        contamination = "results/mutect2/{strand}/{sample}.contamination.table",
        segments = "results/mutect2/{strand}/{sample}.segments.table"
    threads: 2
    shell:
        """
        module load gatk/4.6.0.0

        gatk --java-options "-Xmx4G" CalculateContamination \
            -I {input.pileups} \
            -O {output.contamination} \
            --tumor-segmentation {output.segments}
        """
