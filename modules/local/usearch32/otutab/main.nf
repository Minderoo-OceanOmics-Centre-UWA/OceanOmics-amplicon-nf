process USEARCH32_OTUTAB {
    tag "$reads"
    label 'process_medium'
    container 'sunqiangkun/usearch:v1'

    input:
    path reads
    tuple val(prefix), path(db)

    output:
    tuple val(prefix), path("zotu_table.tsv"), emit: zotu_tsv
    tuple val(prefix), path("lca_input.tsv") , emit: lca_input_tsv
    path "zmap.txt"                          , emit: zmap
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${reads} > input.fa
    usearch -otutab input.fa -zotus ${db} -otutabout zotu_table.tsv -mapout zmap.txt ${args}
    rm input.fa

    # Delete any quote marks from the table and add #ID at the start for lca format
    cat zotu_table.tsv | tr -d '"' | sed '1s/^#OTU ID/#ID/' > lca_input.tsv

    # Also add ZOTU at the start of the zotu table
    cat zotu_table.tsv | tr -d '"' | sed '1s/^#OTU ID/ZOTU/' > tmp.tsv

    # ADD ZOTU_sequence column to zotu_table.tsv
    awk 'BEGIN { FS="\\t"; OFS="\\t" } FNR==NR { if (/^>/) { sub(">", "", \$0); sample=\$0 } else { seq[sample]=\$0 } next } FNR==1 { print \$0, "ZOTU_sequence" } FNR>1 { print \$0, seq[\$1] }' ${db} tmp.tsv > zotu_table.tsv

    rm tmp.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        usearch: \$(usearch --version 2>&1 | head -n 1 | sed 's/usearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
