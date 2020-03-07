# Scatter the DeepVariant workflow on several sample BAMs, then use GLnexus to
# merge the resulting gVCFs.
#
#               +-----------------------------------------------------------+
#               |                                                           |
#               |  DeepVariant_GLnexus.wdl                                  |
#               |                                                           |
#               |       +--------------------------+                        |
#               |       |                          |   sample gVCF          |
#               |   +--->     DeepVariant.wdl      |----+                   |
#               |   |   |                          |    |                   |
#               |   |   +--------------------------+    |    +-----------+  |
#               |   |                                   +---->           |  |
# sample BAMs-------+---> ...                      ...  ----->  GLnexus  +----> project VCF
#               |   |                                   +---->           |  |
#               |   |   +--------------------------+    |    +-----------+  |
#               |   |   |                          |    |                   |
#               |   +--->     DeepVariant.wdl      |----+                   |
#               |       |                          |   sample gVCF          |
#               |       +--------------------------+                        |
#               |                                                           |
#               +-----------------------------------------------------------+
version 1.0

import "DeepVariant.wdl" as dv
import "GLnexus.wdl" as glx

workflow DeepVariant_GLnexus {
    input {
        Array[File]+ bam
        Array[File]+? bai

        # Genomic range(s) to call. Provide at most one of range or ranges_bed.
        # If neither is provided, calls the whole reference genome.
        String? range       # e.g. chr12:111760000-111820000
        File? ranges_bed

        # reference genome & index
        File ref_fasta
        File? ref_fasta_idx

        # DeepVariant settings
        String model_type = "wgs"
        Int shards_per_sample = 32
        Int? gvcf_gq_binsize

        # GLnexus settings
        String glnexus_config = "DeepVariantWGS"

        # pVCF output name
        String output_name
    }

    File? file_none_

    if (!defined(ref_fasta_idx)) {
        # index ref_fasta if needed
        call dv.samtools_faidx {
            input:
                fasta = ref_fasta
        }
    }
    File ref_fasta_idx2 = select_first([ref_fasta_idx, samtools_faidx.fai])

    scatter (i in range(length(bam))) {
        call dv.DeepVariant {
            input:
                ref_fasta = ref_fasta,
                ref_fasta_idx = ref_fasta_idx2,
                range = range,
                ranges_bed = ranges_bed,
                model_type = model_type,
                shards = shards_per_sample,
                gvcf_gq_binsize = gvcf_gq_binsize,
                bam = bam[i],
                bai = if defined(bai) then bai[i] else file_none_
        }
    }

    call glx.GLnexus {
        input:
            gvcf = DeepVariant.gvcf_gz,
            range = range,
            ranges_bed = ranges_bed,
            config = glnexus_config,
            output_name = output_name
    }

    output {
        Array[File]+ gvcf_gz = DeepVariant.gvcf_gz
        File pvcf_gz = GLnexus.pvcf_gz
    }
}

