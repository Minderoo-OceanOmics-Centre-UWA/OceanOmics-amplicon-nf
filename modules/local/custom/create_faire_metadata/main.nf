process CREATE_FAIRE_METADATA {
    tag "$prefix"
    label 'process_low'
    //container 'docker.io/pedrofeijao/pandas-openpyxl:v1.0'
    container 'docker.io/pawelqs/tidyverse_jsonlite_openxlsx:v1'

    input:
    tuple val(prefix), path(taxa_raw), path(taxa_final), path(otu_raw), path(otu_final), path(full_taxa_table)
    path(metadata)

    output:
    path "*xlsx"       , emit: xlsx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "\"${prefix}\""
    def taxa_raw = "\"${taxa_raw}\""
    def taxa_final = "\"${taxa_final}\""
    def otu_raw = "\"${otu_raw}\""
    def otu_final = "\"${otu_final}\""
    def metadata = "\"${metadata}\""
    def full_taxa_table = "\"${full_taxa_table}\""
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(openxlsx)

    # Faster TSV reading
    taxa_raw        <- read.table(${taxa_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    taxa_final      <- read.table(${taxa_final}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    otu_raw         <- read.table(${otu_raw}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    otu_final       <- read.table(${otu_final}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")
    full_taxa_table <- read.table(${full_taxa_table}, sep = "\\t", header = TRUE, stringsAsFactors = FALSE, quote="", comment.char="")

    # Detect ID column
    upper_prefix <- if ("ASV" %in% names(full_taxa_table)) "ASV" else "ZOTU"

    id_col <- upper_prefix

    # Missing IDs
    missing_ids <- setdiff(full_taxa_table[[id_col]], taxa_raw\$seq_id)

    # Only process missing rows
    if (length(missing_ids) > 0) {

        missing_rows <- full_taxa_table[
            full_taxa_table[[id_col]] %in% missing_ids,
        ]
        

        # Build new rows vectorized
        new_rows <- data.table(
            seq_id = missing_rows[[id_col]],
            dna_sequence = missing_rows[[paste0(upper_prefix, "_sequence")]],
            domain = "not applicable",
            phylum = "not applicable",
            class = "not applicable",
            order = "not applicable",
            family = "not applicable",
            genus = "not applicable",
            specificEpithet = "not applicable",
            scientificName = "not applicable",
            scientificNameAuthorship = "not applicable",
            taxonRank = "not applicable",
            taxonID = "not applicable",
            taxonID_db = "not applicable",
            verbatimIdentification = "not applicable",
            accession_id = "not applicable",
            accession_id_ref_db = "not applicable",
            percent_match = "not applicable",
            percent_query_cover = "not applicable",
            confidence_score = "not applicable",
            identificationRemarks = "not applicable"
        )

        # Optional columns
        length_col <- paste0(${prefix}, "_length")

        if (length_col %in% names(missing_rows)) {
            new_rows[[length_col]] <- missing_rows[[length_col]]
        }

        if ("length" %in% names(missing_rows)) {
            new_rows[["length"]] <- missing_rows[["length"]]
        }

        # Add unusual_size only for taxa_final
        new_rows_final <- copy(new_rows)

        if ("unusual_size" %in% names(missing_rows)) {
            new_rows_final[["unusual_size"]] <- missing_rows[["unusual_size"]]
        }

        # Single append
        taxa_raw   <- rbindlist(list(taxa_raw, new_rows), fill = TRUE)
        taxa_final <- rbindlist(list(taxa_final, new_rows_final), fill = TRUE)
    }

    # Workbook
    wb <- loadWorkbook(${metadata})

    write_sheet <- function(wb, sheet, data) {

        if (!(sheet %in% names(wb))) {
            addWorksheet(wb, sheet)
            writeData(wb, sheet, data)
        } else {

            # Faster than reading sheet back in
            existing_rows <- nrow(readWorkbook(wb, sheet = sheet, colNames = TRUE))

            writeData(
                wb,
                sheet = sheet,
                x = data,
                startRow = existing_rows + 2,
                colNames = FALSE
            )
        }
    }

    write_sheet(wb, "taxaRaw", taxa_raw)
    write_sheet(wb, "taxaFinal", taxa_final)
    write_sheet(wb, "otuRaw", otu_raw)
    write_sheet(wb, "otuFinal", otu_final)

    saveWorkbook(
        wb,
        paste0(${prefix}, "_final_faire_metadata.xlsx"),
        overwrite = TRUE
    )

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    openxlsx: ", packageVersion("openxlsx"))), "versions.yml")
    """
}
