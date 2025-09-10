rule prepare_reference:
  input:
    reference_genome=config['ref_index']['genome'],
  output:
    seq="ref/genomes/human/GRCh38/RSEM/reference.seq",
    grp="ref/genomes/human/GRCh38/RSEM/reference.grp",
    ti="ref/genomes/human/GRCh38/RSEM/reference.ti",
  params:
    extra="--gtf {}".format(config["star"]["gtf"]),
  log:
    "logs/rsem/prepare-reference.log",
  shell:
    '''
    module load rsem/1.3.0
    mkdir -p ref/genomes/human/GRCh38/RSEM
    python workflow/scripts/prepare-rsem-reference.py {input.reference_genome} {output.seq} '{params.extra}' {log}
    '''

rule calculate_expression:
  input:
    bam=lambda wc: get_star_output_all_units(wc, fi='transcriptome', use='single'),
    reference="ref/genomes/human/GRCh38/RSEM/reference.seq",
  output:
    genes_results="results/rsem/{sample}.genes.results",
    isoforms_results="results/rsem/{sample}.isoforms.results",
    tbam=temp("results/rsem/{sample}.transcript.bam"),
    gbam=temp("results/rsem/{sample}.genome.bam"),
  params:
    outprefix="results/rsem/{sample}",
    paired_end=lambda w: "--paired-end" if is_paired_end(w.sample) else "",
    extra="-bam --estimate-rspd --output-genome-bam --time --forward-prob 0 --seed 42",
  log:
    "logs/rsem/calculate_expression/{sample}.log",
  shell:
    "ref=$(echo {input.reference} | sed 's/\\..*//'); "
    "module load rsem/1.3.0; "
    "rsem-calculate-expression "
    "--estimate-rspd --output-genome-bam --time --forward-prob 0 --seed 42 "
    "--bam "
    "{params.paired_end} "
    "{input.bam} "
    "$ref "
    "{params.outprefix} "
    "> {log} 2>&1"

rule rsem_generate_genes_tpm_matrix:
  input:
    get_rsem_output_all_units(level="genes"),
  output:
    temp("results/counts/all_gene_tpm_raw.tmp"),
  params:
    extra="TPM",
  log:
    "logs/rsem/generate_gene_tpm_matrix.log",
  shell:
    """
    module load rsem/1.3.0
    perl workflow/scripts/rsem-generate-data-matrix-modified.pl {params.extra} {input} > {output} 2>{log}
    """

rule format_gene_tpm_matrix:
  input:
    "results/counts/all_gene_tpm_raw.tmp",
  output:
    tpm="results/counts/all_gene_tpm_raw.tsv",
  shell:
    '''
    sed 's/\"//g' {input} |
    sed 's/results\/rsem\///g' |
    sed 's/-merged.genes.results//g' |
    sed '1 s/^/gene/' > {output.tpm}
    '''

rule rsem_generate_isoforms_tpm_matrix:
  input:
    get_rsem_output_all_units(level="isoforms"),
  output:
    temp("results/counts/all_isoform_tpm_raw.tmp"),
  params:
    extra="TPM",
  log:
    "logs/rsem/generate_isoform_tpm_matrix.log",
  shell:
    """
    module load rsem/1.3.0
    perl workflow/scripts/rsem-generate-data-matrix-modified.pl {params.extra} {input} > {output} 2>{log}
    """

rule format_isoform_tpm_matrix:
  input:
    "results/counts/all_isoform_tpm_raw.tmp",
  output:
    tpm="results/counts/all_isoform_tpm_raw.tsv",
  shell:
    '''
    sed 's/\"//g' {input} |
    sed 's/results\/rsem\///g' |
    sed 's/-merged.genes.results//g' |
    sed '1 s/^/gene/' > {output.tpm}
    '''

rule rsem_generate_genes_count_matrix:
  input:
    get_rsem_output_all_units(level="genes"),
  output:
    temp("results/counts/all_gene_count_raw.tmp"),
  log:
    "logs/rsem/generate_gene_count_matrix.log",
  shell:
    """
    module load rsem/1.3.0
    rsem-generate-data-matrix {input} > {output} 2>{log}
    """

rule format_gene_count_matrix:
  input:
    "results/counts/all_gene_count_raw.tmp",
  output:
    "results/counts/all_gene_count_raw.tsv",
  shell:
    "sed 's/\"//g' {input} |  "
    "sed 's/results\/rsem\///g' | "
    "sed 's/-merged.genes.results//g' | "
    "sed '1 s/^/gene/' | "
    "sed -e 's/\.[0-9]*//g' > {output}"

rule rsem_generate_isoforms_count_matrix:
  input:
    get_rsem_output_all_units(level="isoforms"),
  output:
    temp("results/counts/all_isoform_count_raw.tmp"),
  log:
    "logs/rsem/generate_isoform_count_matrix.log",
  shell:
    """
    module load rsem/1.3.0
    rsem-generate-data-matrix {input} > {output} 2>{log}
    """

rule format_isoform_count_matrix:
  input:
    "results/counts/all_isoform_count_raw.tmp",
  output:
    "results/counts/all_isoform_count_raw.tsv",
  shell:
    "sed 's/\"//g' {input} |  "
    "sed 's/results\/rsem\///g' | "
    "sed 's/-merged.genes.results//g' | "
    "sed '1 s/^/gene/' | "
    "sed -e 's/\.[0-9]*//g' > {output}"

