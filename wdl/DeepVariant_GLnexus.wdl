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

import "DeepVariant.wdl" as swf

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
        Int? gvcf_gq_binsize

        # GLnexus settings
        String glnexus_config = "DeepVariantWGS"

        # pVCF output name
        String output_name
    }

    File? _file_none

    if (!defined(ref_fasta_idx)) {
        # index ref_fasta if needed
        call swf.samtools_faidx {
            input:
                fasta = ref_fasta
        }
    }
    File ref_fasta_idx2 = select_first([ref_fasta_idx, samtools_faidx.fai])

    scatter (i in range(length(bam))) {
        call swf.DeepVariant as dv { input:
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx2,
            range = range,
            ranges_bed = ranges_bed,
            model_type = model_type,
            gvcf_gq_binsize = gvcf_gq_binsize,
            bam = bam[i],
            bai = if defined(bai) then bai[i] else _file_none
        }
    }

    call GLnexus { input:
        gvcf = dv.gvcf_gz,
        range = range,
        ranges_bed = ranges_bed,
        config = glnexus_config,
        output_name = output_name
    }

    output {
        Array[File] vcf_gz = dv.vcf_gz
        Array[File] gvcf_gz = dv.gvcf_gz
        File pvcf_gz = GLnexus.pvcf_gz
    }
}

task GLnexus {
    input {
        Array[File]+ gvcf
        String? range
        File? ranges_bed
        String config
        String output_name
    }

    command <<<
        set -euxo pipefail
        bed_arg=""
        if [ -n "~{ranges_bed}" ]; then
            bed_arg="--bed ~{ranges_bed}"
        elif [ -n "~{range}" ]; then
            echo "~{range}" | tr :- '\t' | tr -d , > range.bed
            bed_arg="--bed range.bed"
        fi
        glnexus_cli --config "~{config}" --list $bed_arg "~{write_lines(gvcf)}" | \
            bcftools view - | bgzip -c > "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "quay.io/mlin/glnexus:v1.2.2"
        disks: "local-disk 64 HDD"
        cpu: 16
    }

    output {
        File pvcf_gz = "${output_name}.vcf.gz"
    }
}

