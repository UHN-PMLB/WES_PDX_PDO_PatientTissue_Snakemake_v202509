rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    wrapper:
        "0.64.0/bio/star/index"
