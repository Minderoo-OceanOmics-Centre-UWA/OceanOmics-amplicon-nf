process AQUAMAPS {
    tag "$samplesheet"
    label 'process_medium'
    container 'quay.io/biocontainers/pandas:0.23.4--py36hf8a1672_0'

    input:
    tuple val(prefix), path(final_taxa)
    path(aquamaps_db)
    path(samplesheet)

    output:
    tuple val(prefix), path("*csv"), emit: prob_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env python3

    import sqlite3
    import pandas as pd

    def round_to_25_or_75(x):
        base = int(x)
        frac = abs(x - base)

        if frac < 0.5:
            if base >= 0:
                return base + 0.25

            else:
                return base - 0.25

        else:
            if base >= 0:
                return base + 0.75

            else:
                return base - 0.75

    def round_to_closest_50(x):
        base = int(x)
        frac = abs(x - base)

        if frac >= 0.75:
            if base >= 0:
                return base + 1.25

            else:
                return base - 1.25

        elif frac <= 0.25:
            if base >= 0:
                return base - 0.25

            else:
                return base + 0.25

        else:
            return None

    try:
        conn = sqlite3.connect("${aquamaps_db}")
        cursor = conn.cursor()
        samp_df = pd.read_csv("${samplesheet}", sep = ",") # samp_name decimalLatitude decimalLongitude
        taxa_df = pd.read_csv("${final_taxa}", sep = "\\t").fillna("NA") # species
        hits = list(set(taxa_df['species']))
        hits.remove("dropped")
        hits.remove("NA")

        prob_df = pd.DataFrame(columns = list(samp_df['samp_name']),  index = hits)

        for i in hits:
            i_split = i.split(" ")

            if len(i_split) == 2:
                genus = i_split[0]
                species = i_split[1]

                cursor.execute("SELECT SpeciesID FROM speciesoccursum_r WHERE Genus = '" + genus + "' AND Species = '" + species + "'")
                species_id = cursor.fetchone()

                if species_id is not None:
                    species_id = species_id[0]

                    for y in samp_df.index:
                        samp_name = samp_df.at[y, "samp_name"]
                        latitude = samp_df.at[y, "decimalLatitude"]
                        longitude = samp_df.at[y, "decimalLongitude"]

                        if not pd.isna(latitude) and not pd.isna(longitude):
                            rounded_lat = round_to_25_or_75(latitude)
                            rounded_long = round_to_25_or_75(longitude)

                            cursor.execute("SELECT Probability FROM hcaf_species_native WHERE SpeciesID = '" + species_id + "' AND CenterLat = " + str(rounded_lat) + " AND CenterLong = " + str(rounded_long))
                            prob = cursor.fetchone()

                            if prob is not None:
                                prob_df.at[i, samp_name] = prob[0]

                            else:
                                rounded_lat = round_to_closest_50(latitude)
                                rounded_long = round_to_closest_50(longitude)

                                if rounded_lat is not None and rounded_long is not None:
                                    cursor.execute("SELECT Probability FROM hcaf_species_native WHERE SpeciesID = '" + species_id + "' AND CenterLat = " + str(rounded_lat) + " AND CenterLong = " + str(rounded_long))
                                    prob = cursor.fetchone()

                                    if prob is not None:
                                        prob_df.at[i, samp_name] = prob[0]
                                    else:
                                        prob_df.at[i, samp_name] = "species is missing data for this location"
                                else:
                                    prob_df.at[i, samp_name] = "species is missing data for this location"
                        else:
                            prob_df.at[i, samp_name] = "latitude and/or longitude is missing"
                
                else:
                    prob_df.loc[i] = "species not in Aquamaps"

        prob_df.to_csv("aquamaps_" + "${prefix}" + ".csv")
        conn.close()

    except Exception as e:
        print(e)
        conn.close()
    """
}
