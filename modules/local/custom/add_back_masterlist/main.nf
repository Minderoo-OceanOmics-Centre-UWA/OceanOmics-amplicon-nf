process ADD_BACK_MASTERLIST {
    tag "$prefix"
    label 'process_low'

    container 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'

    input:
    tuple val(prefix), path(otu_table), path(lca_table), path(nbc_table)
    tuple val(prefix2), path(removed_seqs)

    output:
    tuple val(prefix), path("*"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    get ASVs in lca_table and otu_table. If an ASV in otu isn't in lca, add it to lca. pull taxa from removed_seqs matching the ASV_sequence. Pull nbc info from nbc_table. Pull count info from otu table
    """
}
