# OceanOmics-amplicon-nf: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and the output directory will contain these subdirectories:

- [01-cutadapt](#cutadapt) - Demultiplexed reads
- [01-fastp](#fastp) - Fastp trimmed reads
- [01-fastqc](#fastqc) - Read QC
- [01-seqkit_stats](#seqkit_stats) - Read stats
- [01-seqtk](#seqtk) - Trimmed reads after seqtk
- [02-dada2](#dada2) - Output files created during ASV workflow
- [02-usearch](#usearch) - Output files created during ZOTU workflow when using usearch mode
- [02-vsearch](#vsearch) - Output files created during ZOTU workflow when using vsearch mode
- [03-filtered_asvs_from_masterlist](#filtered_asvs_from_masterlist) - Fasta files created when using a masterlist.
- [03-lulu](#lulu) - LULU output files and databases
- [04-blast](#blast) - Blastn results
- [04-ocomnbc](#ocomnbc) - Output of OceanOmics Naive Bayes Classifier
- [05-lca](#lca) - Lowest Common Ancestor output
- [06-aquamap](#aquamap) - Aquamap probabilities
- [06-files_with_masterlist_seqs_added_back] - taxonomy files with masterlist sequences added back in.
- [06-phyloseq](#phyloseq) - Phyloseq objects
- [07-faire](#faire) - FAIRe metadata with taxonomy information added.
- [07-multiqc](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [07-proportional_filter](#proportional_filter) - Filtered Phyloseq objects and proportional filtering stats files.
- [07-pipeline info](#pipeline-info) - Report metrics generated during the workflow execution

### 01-cutadapt

<details markdown="1">
<summary>Output files</summary>

- `01-cutadapt/`
  - `all-primers-trimmed`
    - `*fq.gz`: Fq files after primer trimming.
  - `assigned`
    - `*fq.gz`: Fq files after demultiplexing, but before primer trimming.
  - `unknown`
    - `Unknown`
      - `*fq.gz`: Reads that couldn't be assigned to samples.

</details>

### 01-fastp

<details markdown="1">
<summary>Output files</summary>

- `01-fastp/`
  - `*fq.gz`: Fastp trimed reads.

</details>

### 01-fastqc

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### 01-seqkit_stats

<details markdown="1">
<summary>Output files</summary>

- `01-seqkit_stats/`
  - `assigned_seqkit_stats.txt`: Seqkit stats file produced after demultiplexing.
  - `final_seqkit_stats.txt/`: Seqkit stats file produced after primer trimming and filtering.
  - `*prefilter_seqkit_stats.txt/`: Seqkit stats file produced before primer trimming and filtering. Similar to assigned_seqkit_stats.txt, but this is produced even if --skip_demux was used.
  - `raw_seqkit_stats.txt`:Seqkit stats file produced before demultiplexing.
  - `unknown_seqkit_stats.txt`: Seqkit stats file produced after demultiplexing, but only on reads that couldn't be assigned to samples.

</details>

### 01-seqtk

<details markdown="1">
<summary>Output files</summary>

- `01-seqtk`
  - `*fq.gz`: Seqtk trimmed reads.

</details>

### 02-dada2

<details markdown="1">
<summary>Output files</summary>

- `02-dada2/`
  - `plots/`
    - `*png`: QC plots produced during the dada2 workflow.
    - `*asv_track_reads.txt`: Stats file tracking the read retention through different stages of DADA2.
  - `*asv_final_table.tsv`: ASV count table transposed.
  - `*asv_table.csv`: ASV count table.
  - `*asv.fa`: ASV fasta file.
  - `*lca_input.tsv`: ASV table reformated for default LCA script.

</details>

### 02-usearch

<details markdown="1">
<summary>Output files</summary>

- `02-usearch/`
  - `relabeled/`
    - `*fq.gz`: Fastq files with the reads relabled.
  - `filtered/`
    - `*fq.gz`: Fastq files after adapterRemoval.
  - `*zotus.fa`: Final zOTUs.
  - `*Unq.fa.gz`: Dereplicated reads.
  - `*zotu_table.tsv`: zOTU count table.
  - `*lca_input.tsv`: zOTU count table reformated.

</details>

### 02-vsearch

<details markdown="1">
<summary>Output files</summary>

- `02-vsearch/`
  - `relabeled/`
    - `*fq.gz`: Fastq files with the reads relabled.
  - `filtered/`
    - `*fq.gz`: Fastq files after adapterRemoval.
  - `*zotus.fa`: Final zOTUs.
  - `*centroids.fa`: Denoised reads.
  - `*Unq.fa.gz`: Dereplicated reads.
  - `*zotu_table.tsv`: zOTU count table.
  - `*lca_input.tsv`: zOTU count table reformated.

</details>

### 03-filtered_asvs_from_masterlist

<details markdown="1">
<summary>Output files</summary>

- `03-filtered_asvs_from_masterlist/`
  - `*non_masterlist_seqs_asv.fa`: ASVs found in this data but not in masterlist.
  - `*removed_seqs_asv.fa`: ASVs found in this data and in masterlist.

</details>

### 03-lulu

<details markdown="1">
<summary>Output files</summary>

- `03-lulu/`
  - `*asv_db/`
    - `asv*`: Blast database created from ASVs.
  - `*asv_curated_table.tab`: ASV count table after LULU curation.
  - `*asv_lulu_map.tab`: Stat file produced by LULU.
  - `*asv_match_list.txt`: ASV matchlist created by BLASTing ASVs with each other. This was used by LULU.
  - `*curated_asv.fa`: ASV fasta file after LULU curation.

</details>

### 04-blast

<details markdown="1">
<summary>Output files</summary>

- `04-blast/`
  - `*blastn_results.txt`: BLAST results file presented as a table.
  - `*default_format_blastn_results.txt`: BLAST results file in the default format.

</details>

### 04-ocomnbc

<details markdown="1">
<summary>Output files</summary>

- `04-ocomnbc/`
  - `*ocom_nbc_output.tsv`: The Naive Bayes classifier output file..

</details>

### 05-lca

<details markdown="1">
<summary>Output files</summary>

- `05-lca/`
  - `*lca_with_fishbase_output.tsv`: LCA output file that contains info on what happened during the LCA step.
  - `*taxa_final.tsv`: taxa table in the FAIRe format with LCA results.
  - `*taxa_raw.tsv`: taxa table in the FAIRe format with BLAST results (up to 10 hits per ASV).

</details>

### 06-aquamap

<details markdown="1">
<summary>Output files</summary>

- `06-aquamap/`
  - `*csv`: Aquamaps probability file that has the probability (from 0 - 1) that a species might be found at the sampling locations.

</details>

### 06-files_with_masterlist_seqs_added_back

<details markdown="1">
<summary>Output files</summary>

- `06-files_with_masterlist_seqs_added_back/`
  - `*lca_with_fishbase_output.tsv`: LCA output file with masterlist sequences added back.
  - `*ocom_nbc_output.tsv`: NBC output file with masterlist sequences added back.
  - `*taxa_final.tsv`: Taxonomy file with masterlist sequences added back.

</details>

### 06-phyloseq

<details markdown="1">
<summary>Output files</summary>

- `06-phyloseq/`
  - `*final_taxa.tsv`: Final taxa file before being added to the phyloseq object.
  - `*flagged_phyloseq.rds`: The Phyloseq object.

</details>

### 07-faire

<details markdown="1">
<summary>Output files</summary>

- `07-faire/`
  - `*final_faire_metadata.xlsx`: FAIRe formated data.

</details>

### 07-multiqc

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### 07-proportional_filter

<details markdown="1">
<summary>Output files</summary>

- `07-proportional_filter/`
  - `*faire_taxa_filtered.tsv`: The taxa file in the FAIRe format after filtering.
  - `*OTU_filtered.tsv`: The count table after filtering.
  - `*phyloseq_filtered.rds`: The phyloseq object after filtering.
  - `*phyloseq_taxa_filtered.tsv`: The taxa table after filtering.
  - `*proportional_stats.txt`: Stats file that explains what happened during the filter step.

</details>

### 07-pipeline info

<details markdown="1">
<summary>Output files</summary>

- `07-pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
