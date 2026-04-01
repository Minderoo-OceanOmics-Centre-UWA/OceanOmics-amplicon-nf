process BLAST_BLASTN {
    tag "$prefix"
    label 'process_medium'

    container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

    input:
    tuple val(prefix), path(fasta)
    path db

    output:
    tuple val(prefix), path("*blastn_results.txt"), emit: txt
    tuple val(prefix), path("*default_format.txt"), emit: default_format
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${prefix}"
    ////"-outfmt '6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp'",
    """
    DB=`echo \$(ls *.ndb | sed 's/\\.ndb\$//')`

    blastn \\
        -num_threads $task.cpus \\
        -db "\$DB" \\
        -query $fasta \\
        $args \\
        -outfmt 11 \\
        -out ${prefix}_${fasta}_blastn_results.asn

    blast_formatter -archive ${prefix}_${fasta}_blastn_results.asn -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" > ${prefix}_${fasta}_blastn_results.txt
    blast_formatter -archive ${prefix}_${fasta}_blastn_results.asn > ${prefix}_${fasta}_blastn_results_default_format.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
