/*
 * Demultiplex with Cutadapt
 */

include { ALT_CUTADAPT                   } from '../../modules/local/alt_cutadapt/main.nf'
include { CREATE_DEMUX_DEPENDENCIES      } from '../../modules/local/custom/createdemuxdependencies/main.nf'
include { RENAME                         } from '../../modules/local/custom/rename/main.nf'
include { SEQKIT_STATS as \
            ASSIGNED_STATS;
            SEQKIT_STATS as \
            UNKNOWN_STATS;
            SEQKIT_STATS as \
            RAW_STATS                    } from '../../modules/local/seqkit_stats/main.nf'
include { CONCAT                         } from '../../modules/local/custom/concat/main.nf'

workflow ALT_CUTADAPT_WORKFLOW {
    take:
    ch_input    // channel: [ samplesheet.csv ]
    ch_ulimit   // channel: [ ulimit ]
    ch_reads

    main:
    ch_versions = Channel.empty()
    //ch_reads.view()

    // MODULE: Demultiplex
    ALT_CUTADAPT (
        ch_reads,
        params.ulimit
    )
    ch_versions = ch_versions.mix(ALT_CUTADAPT.out.versions)

    // MODULE: Rename the Cutadapt output files
    //RENAME (
    //    ALT_CUTADAPT.out.reads,
    //    CREATE_DEMUX_DEPENDENCIES.out.sample_rename
    //)
    //ch_versions = ch_versions.mix(RENAME.out.versions)

    // MODULE: Check stats of reads that couldn't be assigned to samples after demultiplexing
    //UNKNOWN_STATS (
    //    RENAME.out.unknown,
    //    "unknown"
    //)
    //ch_versions = ch_versions.mix(UNKNOWN_STATS.out.versions)

    //ch_demuxed_reads = CONCAT.out.reads
    ch_demuxed_reads = ALT_CUTADAPT.out.reads


     // MODULE: Check stats of reads assigned to samples after demultiplexing
    //ASSIGNED_STATS (
    //    ch_demuxed_reads,
    //    "assigned"
    //)
    //ch_versions = ch_versions.mix(ASSIGNED_STATS.out.versions)

    emit:
    reads           = ch_demuxed_reads                    // channel: [ val(meta), reads ]
    //assigned_stats  = ASSIGNED_STATS.out.stats            // channel: [ val(meta), stats ]
    //unknown_stats   = UNKNOWN_STATS.out.stats             // channel: [ val(meta), stats ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
