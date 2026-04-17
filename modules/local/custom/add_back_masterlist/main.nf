process ADD_BACK_MASTERLIST {
    tag "$prefix"
    label 'process_low'

    //container 'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0'
    //container 'quay.io/biocontainers/pandas:1.5.2'
    container 'quay.io/fhcrc-microbiome/biopython-pandas:3dd2f94'

    input:
    tuple val(prefix) , path(otu_table), path(lca_table), path(nbc_table), path(taxa_final)
    tuple val(prefix2), path(removed_seqs)

    output:
    tuple val(prefix), path("with_masterlist_seqs_*final_table*"), path("with_masterlist_seqs_*lca*"), path("with_masterlist_seqs_*nbc*"), emit: tsv
    tuple val(prefix), path("with_masterlist_seqs_*taxa_final*"), emit: taxa_final

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    from Bio import SeqIO

    upper_prefix = "${prefix}".split("_")[0].upper()

    lca_table = pd.read_csv("${lca_table}", sep = "\\t") # domain phylum class order family genus species OTU asv_blast_length numberOfUnq_BlastHits %ID queryCoverage species_in_LCA sources OcOm_2101_17_1 OcOm_2101_17_2 OcOm_2101_17_3 OcOm_2101_17_4 OcOm_2101_17_5 OcOm_2101_18_1 OcOm_2101_18_2 OcOm_2101_18_3 OcOm_2101_18_4 OcOm_2101_18_5
    otu_table = pd.read_csv("${otu_table}", sep = "\\t") # ASV OcOm_2101_17_1 OcOm_2101_17_2 OcOm_2101_17_3 OcOm_2101_17_4 OcOm_2101_17_5 OcOm_2101_18_1 OcOm_2101_18_2 OcOm_2101_18_3 OcOm_2101_18_4 OcOm_2101_18_5 ASV_sequence
    nbc_table = pd.read_csv("${nbc_table}", sep = "\\t") # Gene Label Genus.prediction Genus.score Species.prediction Species.score
    taxa_final = pd.read_csv("${taxa_final}", sep = "\\t") # 'seq_id', 'dna_sequence', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'specificEpithet', 'scientificName', 'scientificNameAuthorship', 'taxonRank', 'taxonID', 'taxonID_db', 'verbatimIdentification', 'accession_id', 'accession_id_ref_db', 'percent_match', 'percent_query_cover', 'confidence_score', 'identificationRemarks', 'length'

    #removed_seqs = read_fasta("${removed_seqs}") # >Eukaryota;Chordata;Actinopteri;Perciformes/Serranoidei;Anthiadidae;Pseudanthias;dropped
    samples = otu_table.columns.tolist()[1:-1]

    removed_dict = {}
    for record in SeqIO.parse("${removed_seqs}", "fasta"):
        removed_dict[str(record.seq)] = record.id

    lca_asvs = lca_table["OTU"].tolist()
    otu_asvs = otu_table[upper_prefix].tolist()

    for asv in otu_asvs:
        if asv not in lca_asvs:
            curr_row = otu_table[otu_table[upper_prefix] == asv]

            ASV_sequence = curr_row[upper_prefix + "_sequence"].item()

            if ASV_sequence in list(removed_dict.keys()):
                taxa = removed_dict[ASV_sequence]
                taxa = taxa.lstrip(">")
                split_taxa = taxa.split(";")
                domain = split_taxa[0]
                phylum = split_taxa[1]
                clss = split_taxa[2]
                order = split_taxa[3]
                family = split_taxa[4]
                genus = split_taxa[5]
                species = split_taxa[6]
            else:
                domain = "NA"
                phylum = "NA"
                clss = "NA"
                order = "NA"
                family = "NA"
                genus = "NA"
                species = "NA"

            new_row = {
                "domain": domain, 
                "phylum":phylum, 
                "class": clss, 
                "order": order, 
                "family": family, 
                "genus": genus, 
                "species": species, 
                "OTU": asv, 
                "length": "NA", 
                "numberOfUnq_BlastHits": "NA", 
                "%ID": "NA",
                "queryCoverage": "NA",
                "species_in_LCA": "NA",
                "sources": "NA" 
            }

            # new_row (sample columns) fill with otu_table
            for sam in samples:
                read_count = curr_row[sam].item()
                new_row[sam] = read_count

            lca_table = pd.concat([lca_table, pd.DataFrame([new_row])], ignore_index=True)

            try:
                epithet = species.split(" ")[1]
            except IndexError:
                epithet = ""

            new_row = {
                "seq_id": asv,
                "dna_sequence": ASV_sequence,
                "domain": domain, 
                "phylum":phylum, 
                "class": clss, 
                "order": order, 
                "family": family, 
                "genus": genus, 
                "specificEpithet": epithet,
                "scientificName": species,
                "scientificNameAuthorship": "NA",
                "taxonRank": "NA",
                "taxonID": "NA",
                "taxonID_db": "NA",
                "verbatimIdentification": "NA",
                "accession_id": "NA",
                "accession_id_ref_db": "NA",
                "percent_match": "NA",
                "percent_query_cover": "NA",
                "confidence_score": "NA",
                "identificationRemarks": "NA",
                "length": "NA"
            }

            taxa_final = pd.concat([taxa_final, pd.DataFrame([new_row])], ignore_index=True)

    lca_table.to_csv("with_masterlist_seqs_${lca_table}", sep = "\\t", index = False)
    otu_table.to_csv("with_masterlist_seqs_${otu_table}", sep = "\\t", index = False)
    nbc_table.to_csv("with_masterlist_seqs_${nbc_table}", sep = "\\t", index = False)
    taxa_final.to_csv("with_masterlist_seqs_${taxa_final}", sep = "\\t", index = False)
    """
}
