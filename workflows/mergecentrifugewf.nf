include { mergeCentrifuge } from '../modules/merge_centrifuge'

workflow MERGECENTRIFUGE {
    take:
        ch_centrifuge_reports

    main:
        mergeCentrifuge(ch_centrifuge_reports)

    emit:
        ch_merged_centrifuge = mergeCentrifuge.out
}
