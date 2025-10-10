rule fastp:
   input:
      f1 = get_fq1,
      f2 = get_fq2,
   output:
      trim_f1 = "results/fastp/{sample}-trim.1.fq.gz",
      trim_f2 = "results/fastp/{sample}-trim.2.fq.gz",
      json = "results/fastp/{sample}.fastp.json",
      html = "results/fastp/{sample}.fastp.html",
   threads: 4
   shell:
       """
       module load fastp/0.23.1

       fastp -i {input.f1} -I {input.f2} -o {output.trim_f1} -O {output.trim_f2} -h {output.html} -j {output.json}
       """


rule deconvolutexengsort:
   input:
      f1 =  "results/fastp/{sample}-trim.1.fq.gz",
      f2 =  "results/fastp/{sample}-trim.2.fq.gz",
   output:
      graftf1 =  "results/xengsort/{sample}-graft.1.fq.gz",
      graftf2 =  "results/xengsort/{sample}-graft.2.fq.gz",
   params:
      xengsortidx = config["ref"]["xengsortidx"],
      xengsortcontainer = config['env']['xengsort'],
      sampleid = "{sample}",
      sample_type = lambda wildcards: get_sample_type(wildcards),
   threads: 4
   shell:
    r"""
    module load apptainer/1.0.2
    module load pigz/2.6
    mkdir -p results/xengsort
    mkdir -p tmp

    if [[ "{params.sample_type}" == "PDX" ]]; then

      zcat {input.f1} > tmp/{params.sampleid}_R1.fastq
      zcat {input.f2} > tmp/{params.sampleid}_R2.fastq

      apptainer run {params.xengsortcontainer} \
        xengsort classify \
        --index {params.xengsortidx} \
        --fastq tmp/{params.sampleid}_R1.fastq \
        --pairs tmp/{params.sampleid}_R2.fastq \
        --prefix results/xengsort/{params.sampleid} \
        --compression none \
        -T {threads} \
        --progress

      rm tmp/{params.sampleid}_R1.fastq tmp/{params.sampleid}_R2.fastq

      for suffix in graft neither both ambiguous host; do
        pigz -{threads} results/xengsort/{params.sampleid}-${{suffix}}.1.fq
        pigz -{threads} results/xengsort/{params.sampleid}-${{suffix}}.2.fq
      done

    else
      ln -sf $(realpath {input.f1}) {output.graftf1}
      ln -sf $(realpath {input.f2}) {output.graftf2}
    fi
    """


rule bwa_mem2_index:
    input:
        ref = config["ref_index"]["genome"]
    output:
        idx = expand("{ref}.{ext}", 
                     ref=config["ref_index"]["genome"], 
                     ext=["0123","amb","ann","bwt.2bit.64","pac"])
    params:
        bwamem2container = config['env']['bwa_mem2']
    threads: 4
    shell:
        """
        module load apptainer/1.0.2

        apptainer run {params.bwamem2container} \
        bwa-mem2 index {input.ref}
        """

rule bwa_mem2_align_pe:
    input:
        fq1 = "results/xengsort/{sample}-graft.1.fq.gz",
        fq2 = "results/xengsort/{sample}-graft.2.fq.gz",
        idx = expand("{ref}.{ext}", 
                     ref=config["ref_index"]["genome"], 
                     ext=["0123","amb","ann","bwt.2bit.64","pac"])
    output:
        bam = "results/bwa/pe/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        ref = config["ref_index"]["genome"],
        bwamem2container = config['env']['bwa_mem2'],
    threads: 16
    shell:
        """
        module load apptainer/1.0.2
        module load samtools/1.20

        mkdir -p results/bwa/pe/{wildcards.sample}

        apptainer run {params.bwamem2container} \
        bwa-mem2 mem -t {threads} \
		      -R "@RG\tID:{wildcards.sample}\tLB:Exome\tSM:{wildcards.sample}\tPL:ILLUMINA" \
		      {params.ref} {input.fq1} {input.fq2} \
            | samtools sort -@ {threads} -o {output.bam}
        """

rule bwa_mem2_align_se:
    input:
        fq1 = "results/xengsort/{sample}-graft.1.fq.gz",
        idx = expand("{ref}.{ext}", 
                     ref=config["ref_index"]["genome"], 
                     ext=["0123","amb","ann","bwt.2bit.64","pac"])
    output:
        bam = "results/bwa/se/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        ref = config["ref_index"]["genome"],
        bwamem2container = config['env']['bwa_mem2'],
    threads: 16
    shell:
        """
        module load apptainer/1.0.2
        module load samtools/1.20

        mkdir -p results/bwa/se/{wildcards.sample}

        apptainer run {params.bwamem2container} \
        bwa-mem2 mem -t {threads} \
		      -R "@RG\tID:{wildcards.sample}\tLB:Exome\tSM:{wildcards.sample}\tPL:ILLUMINA" \
		      {params.ref} {input.fq1} \
            | samtools sort -@ {threads} -o {output.bam}
        """

rule index_coord:
    input:
        bam = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam",
    output:
        bai = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam.bai",
    shell:
        """
        module load samtools/1.20
        samtools index {input.bam} {output.bai}
        """

rule picard_align_matrix:
    input:
        bam = "results/bwa/{strand}/{sample}/Aligned.sortedByCoord.out.bam",
    params:
        ref = config['ref_index']['genome'],
    output:
        align_matrix = "results/picard/{strand}/{sample}_picard_align_matrix.txt",
    shell:
        """
        module load picard/2.10.9

        java -jar $picard_dir/picard.jar CollectAlignmentSummaryMetrics \
            R={params.ref} \
            I={input.bam} \
            O={output.align_matrix}
        """
