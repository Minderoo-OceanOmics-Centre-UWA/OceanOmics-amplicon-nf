process RENAME {
    tag ""
    label 'process_medium'
    container 'adbennett/mmv:v0.1'

    input:
    tuple val(prefix), path(reads)
    path sample_rename

    output:
    tuple val(prefix), path("Assigned/*.fq.gz"), emit: reads
    tuple val(prefix), path("Unknown/*.fq.gz") , emit: unknown
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Here we reference the rename pattern file, which contains the information that maps the
    # index ID pairings with the sample IDs
    mmv < ${sample_rename} -g

    # move the unnamed and unknowns into separate folders
    mkdir -p Assigned Unknown
    mv *unknown*.fq.gz Unknown

    # Get sample names
    samples=\$(awk '{split(\$2, a, "\\\\.[#]*1\\\\.fq\\\\.gz"); print a[1]}' "${sample_rename}")

    # Moved assigned files to Assigned Dir
    for sample in \$samples; do
        if [[ \${sample} == *"_forward"* ]]; then
            sam=\$(echo \${sample} | sed 's/_forward//g')
            zcat \${sam}_forward.1.fq.gz \${sam}_reverse.2.fq.gz > Assigned/\${sam}.R1.fq
            zcat \${sam}_forward.2.fq.gz \${sam}_reverse.1.fq.gz > Assigned/\${sam}.R2.fq
            rm \${sam}_forward*fq.gz
            rm \${sam}_reverse*fq.gz
        fi
    done

    gzip Assigned/*
    

    # Any leftover .fq.gz files are unknown samples
    if compgen -G "*.fq.gz"; then
        mv *.fq.gz Unknown
    fi

    VERSION=\$(apt-cache show mmv | grep ^Version: | sed 's/Version: //')
    echo '"${task.process}":
        mmv: \$VERSION' > versions.yml
    """
}
