process FILT_WITH_MASTERLIST {
    tag "$prefix"
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'

    input:
    tuple val(prefix), path(fasta)
    path(masterlist)

    output:
    tuple val(prefix), path("non_masterlist_seqs_*"), emit: fasta
    tuple val(prefix), path("removed_seqs_*")      , emit: removed_seqs

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat ${masterlist} | grep -v '>' > seqs.txt
    seqkit grep -v -s -f seqs.txt ${fasta} > non_masterlist_seqs_${fasta}
    seqkit grep -s -f seqs.txt ${fasta} > tmp_${fasta}
    cat tmp_${fasta} | grep -v '>' > tmp_seqs.txt
    seqkit grep -s -f tmp_seqs.txt ${masterlist} > removed_seqs_${fasta}
    """
}
