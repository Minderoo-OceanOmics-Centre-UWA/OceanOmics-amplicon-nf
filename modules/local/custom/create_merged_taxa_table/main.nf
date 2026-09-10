process CREATE_MERGED_TAXA_TABLE {
    tag ""
    label 'process_low'
    container 'rocker/tidyverse:4.3.0'

    input:
    path faire_files

    output:
    path "*.xlsx", emit: xlsx

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    # CHECK IF THERE IS A FILE WITH curateddb AND A FILE WITH nt
    if () {
        # LOOP THROUGH ALL FILES THAT HAVE CURATEDDB
        for () {
            # EXTRACT ASSAY FROM FILE NAME AND ASV/ZOTU
            # GET nt VERSION OF FILE
            # DO MAIN CODE
            # OUTPUT TAXA FILE

            cur_path <- "16SFishD/07-faire/OcOm_2518_16SFishD_asv_curateddb_final_faire_metadata.xlsx"
            nt_path  <- "16SFishD/07-faire/OcOm_2518_16SFishD_asv_nt_final_faire_metadata.xlsx"

            cur_taxa <- read_excel(cur_path, sheet = "taxaFinal", skip = 2)
            nt_taxa  <- read_excel(nt_path, sheet = "taxaFinal", skip = 2)

            cur_pids        <- c()
            cur_qcovs       <- c()
            cur_asses       <- c()
            nt_pids         <- c()
            nt_qcovs        <- c()
            nt_asses        <- c()
            cur_ids         <- c()
            nt_ids          <- c()
            greater_pids    <- c()
            greater_qcovs   <- c()
            greater_asses   <- c()
            greater_sources <- c()
            for (asv in cur_taxa\$seq_id) {
                cur_pid  <- cur_taxa[cur_taxa\$seq_id == asv, "percent_match"][[1]]
                cur_qcov <- cur_taxa[cur_taxa\$seq_id == asv, "percent_query_cover"][[1]]
                nt_pid   <- nt_taxa[nt_taxa\$seq_id == asv, "percent_match"][[1]]
                nt_qcov  <- nt_taxa[nt_taxa\$seq_id == asv, "percent_query_cover"][[1]]
                nt_rank  <- nt_taxa[nt_taxa\$seq_id == asv, "taxonRank"][[1]]
                cur_id   <- cur_taxa[cur_taxa\$seq_id == asv, "accession_id"][[1]]
                nt_id    <- nt_taxa[nt_taxa\$seq_id == asv, "accession_id"][[1]]
                if (nt_rank == "not applicable" | nt_rank == "not applicable: domain level dropped") {
                    nt_ass <- "not applicable"
                } else if (nt_rank == "species") {
                    nt_ass <- nt_taxa[nt_taxa\$seq_id == asv, "scientificName"][[1]]
                } else {
                    nt_ass <- nt_taxa[nt_taxa\$seq_id == asv, nt_rank][[1]]
                }
                cur_rank <- cur_taxa[cur_taxa\$seq_id == asv, "taxonRank"][[1]]
                if (cur_rank == "not applicable" | cur_rank == "not applicable: domain level dropped") {
                    cur_ass <- "not applicable"
                } else if (cur_rank == "species") {
                    cur_ass <- cur_taxa[cur_taxa\$seq_id == asv, "scientificName"][[1]]
                } else {
                    cur_ass <- cur_taxa[cur_taxa\$seq_id == asv, cur_rank][[1]]
                }
                
                cur_pids  <- c(cur_pids, cur_pid)
                cur_qcovs <- c(cur_qcovs, cur_qcov)
                cur_asses <- c(cur_asses, cur_ass)
                nt_pids   <- c(nt_pids, nt_pid)
                nt_qcovs  <- c(nt_qcovs, nt_qcov)
                nt_asses  <- c(nt_asses, nt_ass)
                cur_ids   <- c(cur_ids, cur_id)
                nt_ids    <- c(nt_ids, nt_id)
                
                if (! is.na(nt_pid) & ! is.na(cur_pid) & nt_pid != "not applicable") {
                    if (cur_pid == "not applicable") {
                        greater_pids    <- c(greater_pids, nt_pid)
                        greater_qcovs   <- c(greater_qcovs, nt_qcov)
                        greater_asses   <- c(greater_asses, nt_ass)
                        greater_sources <- c(greater_sources, "nt")
                    } else if (as.numeric(nt_pid) > as.numeric(cur_pid)) {
                        greater_pids    <- c(greater_pids, nt_pid)
                        greater_qcovs   <- c(greater_qcovs, nt_qcov)
                        greater_asses   <- c(greater_asses, nt_ass)
                        greater_sources <- c(greater_sources, "nt")
                    } else {
                        greater_pids    <- c(greater_pids, cur_pid)
                        greater_qcovs   <- c(greater_qcovs, cur_qcov)
                        greater_asses   <- c(greater_asses, cur_ass)
                        greater_sources <- c(greater_sources, "curateddb")
                    }
                } else {
                    greater_pids    <- c(greater_pids, cur_pid)
                    greater_qcovs   <- c(greater_qcovs, cur_qcov)
                    greater_asses   <- c(greater_asses, cur_ass)
                    greater_sources <- c(greater_sources, "curateddb")
                }
            }
            cur_taxa\$curateddb_assignment          <- cur_asses
            cur_taxa\$curateddb_id                  <- cur_ids
            cur_taxa\$curateddb_percent_match       <- cur_pids
            cur_taxa\$curateddb_percent_query_cover <- cur_qcovs
            cur_taxa\$nt_assignment                 <- nt_asses
            cur_taxa\$nt_id                         <- nt_ids
            cur_taxa\$nt_percent_match              <- nt_pids
            cur_taxa\$nt_percent_query_cover        <- nt_qcovs
            cur_taxa\$best_assignment               <- greater_asses
            cur_taxa\$best_assignment_database      <- greater_sources
            cur_taxa\$best_percent_match            <- greater_pids
            cur_taxa\$best_percent_query_cover      <- greater_qcovs

            write_csv(cur_taxa, file = paste0("taxaFinal_", mode, "_", assay, ".csv"))

            ogs   <- c()
            names <- c()
            for (asv in cur_taxa\$seq_id) {
                if (! is.na(cur_taxa[cur_taxa\$seq_id == asv, "accession_id_ref_db"])) {
                    if (cur_taxa[cur_taxa\$seq_id == asv, "accession_id_ref_db"] == "OG") {
                        acc_id <- cur_taxa[cur_taxa\$seq_id == asv, "accession_id"][[1]]
                        og     <- unlist(strsplit(acc_id, split = ".", fixed = TRUE))[1]
                        name   <- cur_taxa[cur_taxa\$seq_id == asv, "scientificName"][[1]]
                        if (! og %in% ogs & name != "dropped") {
                            ogs   <- c(ogs, og)
                            names <- c(names, name)
                        }
                    }
                }
            }
            df = data.frame(ogs, names)
        }
    }
    """
}
