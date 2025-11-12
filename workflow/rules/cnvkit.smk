rule cnvkit_prep_targets:
    input:
        bed = config["cnv"]["target_bed"],
        genome_len = config["ref_index"]["genome_length"],
    output:
        targets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.bed"
    shell:
        """
        module load bedtools/2.27.1

        mkdir -p $(dirname {output.targets})

        # Step 1: sort intervals by chromosome and coordinate
        # Step 2: merge overlapping or adjacent exons
        # Step 3: pad merged intervals by 100 bp (slop)
        # Step 4: merge again in case padding created overlaps
        bedtools sort -i {input.bed} \
        | bedtools merge -i - \
        | bedtools slop -i - -g {input.genome_len} -b 100 \
        | bedtools merge -i - > {output.targets}
        """

rule cnvkit_target_annotated:
    input:
        targets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.bed",
        refFlat = config["ref_index"]["refFlat"]
    output:
        annotated_targets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.annotated.bed"
    shell:
        """
        module load CNVkit/0.9.9
        mkdir -p $(dirname {output.annotated_targets})

        cnvkit.py target {input.targets} \
            --annotate {input.refFlat} \
            -o {output.annotated_targets}
        """

rule cnvkit_access:
    input:
        ref = config["ref_index"]["genome"]   # path to hg38.fa used for alignment
    output:
        access = "ref/genomes/human/WES_target/cnvkit/hg38.access.bed"
    shell:
        """
        module load CNVkit/0.9.9

        mkdir -p $(dirname {output.access})

        cnvkit.py access {input.ref} -o {output.access}
        """

rule cnvkit_antitargets:
    input:
        targets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.annotated.bed",
        access = "ref/genomes/human/WES_target/cnvkit/hg38.access.bed"
    output:
        antitargets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.antitargets.bed"
    shell:
        """
        module load CNVkit/0.9.9

        mkdir -p $(dirname {output.antitargets})

        cnvkit.py antitarget {input.targets} -g {input.access} -o {output.antitargets}
        """

rule cnvkit_flat_reference:
    input:
        fasta = config["ref_index"]["genome"],
        targets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.merged.pad100.annotated.bed",
        antitargets = "ref/genomes/human/WES_target/cnvkit/SureSelectXT_HS_V8.antitargets.bed"
    output:
        reference = "ref/genomes/human/WES_target/cnvkit/flat_reference.cnn"
    shell:
        """
        module load CNVkit/0.9.9
        mkdir -p $(dirname {output.reference})

        cnvkit.py reference \
            --fasta {input.fasta} \
            --targets {input.targets} \
            --antitargets {input.antitargets} \
            -o {output.reference}
        """

rule cnvkit_batch:
    input:
        bam = "results/gatk/pe/{sample}.recal.bam",
        reference = "ref/genomes/human/WES_target/cnvkit/flat_reference.cnn"
    output:
        cnr = "results/cnvkit/{sample}/{sample}.cnr"
    params:
        outdir = "results/cnvkit/{sample}"
    threads: 8
    shell:
        """
        module load CNVkit/0.9.9
        mkdir -p $(dirname {output.cnr})

        cnvkit.py batch {input.bam} \
            --reference {input.reference} \
            --output-dir {params.outdir} \
            --method hybrid \
            --drop-low-coverage \
            --processes {threads}
        """

rule cnvkit_segment:
    input:
        cnr = "results/cnvkit/{sample}/{sample}.cnr"
    output:
        cns = "results/cnvkit/{sample}/{sample}.cns"
    threads: 4
    shell:
        """
        module load CNVkit/0.9.9
	module load R/4.4.1

        cnvkit.py segment {input.cnr} \
            --method cbs \
            --processes {threads} \
            -o {output.cns}
        """

rule cnvkit_call:
    input:
        cns = "results/cnvkit/{sample}/{sample}.cns",
        vcf = "results/mutect2/pe/{sample}.pass_rare.vcf"   # path to SNV file
    output:
        call = "results/cnvkit/{sample}/{sample}.call.cns"
    params:
        ploidy = 2,
        purity = 0.7
    shell:
        """
        module load CNVkit/0.9.9
        cnvkit.py call {input.cns} \
            -v {input.vcf} \
            --center median \
            --ploidy {params.ploidy} \
            --purity {params.purity} \
            -o {output.call}
        """

rule cnvkit_gene_cnv:
    input:
        call = "results/cnvkit/{sample}/{sample}.call.cns",
        cnr = "results/cnvkit/{sample}/{sample}.cnr"
    output:
        gene_cnv = "results/cnvkit/{sample}/{sample}.gene_cnv.tsv"
    shell:
        """
        module load CNVkit/0.9.9
        cnvkit.py genemetrics {input.cnr} -s {input.call} \
            -o {output.gene_cnv}
        """

rule cnvkit_export_seg:
    input:
        call = "results/cnvkit/{sample}/{sample}.call.cns"
    output:
        seg = "results/cnvkit/{sample}/{sample}.seg"
    shell:
        """
        module load CNVkit/0.9.9
        cnvkit.py export seg {input.call} -o {output.seg}
        """

rule cnvkit_scatter_after_call:
    input:
        cnr = "results/cnvkit/{sample}/{sample}.cnr",
        call = "results/cnvkit/{sample}/{sample}.call.cns",
        vcf = "results/mutect2/pe/{sample}.filtered.vcf.gz"
    output:
        scatter = "results/cnvkit/{sample}/{sample}_scatter.pdf",
        diagram = "results/cnvkit/{sample}/{sample}_diagram.pdf"
    shell:
        """
        module load CNVkit/0.9.9

        # scatter: CN log2 ratios + corrected segments + BAF
        cnvkit.py scatter {input.cnr} -s {input.call} -v {input.vcf} -o {output.scatter}

        # diagram: genome-wide CN gains/losses overview
        cnvkit.py diagram {input.call} -o {output.diagram}
        """

